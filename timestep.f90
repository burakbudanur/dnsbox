#include "macros.h"
module timestep
    ! Time stepping routines
    use numbers
    use openmpi
    use io
    use parameters
    use fieldio
    use fftw
    use vfield
    use rhs

    integer(i4) :: ncorr_last = 1
    real(dp)    :: error_last = 0, courant

    integer(i4)    :: steps_ch
    character(255) :: steps_file = 'steps.gp'
    logical :: steps_written = .false.
    
    contains

!==============================================================================

    subroutine timestep_precorr(vel_vfieldxx, vel_vfieldk, fvel_vfieldk)
        
        real(dp), intent(inout) :: vel_vfieldxx(:, :, :, :)
        complex(dpc), intent(inout) :: vel_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: fvel_vfieldk(:, :, :, :)

        real(dp)     :: timestep_prefieldxx(nyy, nzz_perproc, nxx, 3)
        complex(dpc), dimension(nx_perproc, ny_half, nz, 3) :: &
            timestep_prefieldk, timestep_corfieldk, timestep_nonlinterm_next, &
            timestep_nonlinterm_prev

        ! semi-implicit predictor-corrector method
        _indices
        integer(i4) :: c
        real(dp)    :: invdt, norm, error 

        ! initial rhs
        call rhs_nonlin_term(vel_vfieldxx, vel_vfieldk, fvel_vfieldk)
        timestep_nonlinterm_prev = fvel_vfieldk

        ! add the linear term to the rhs for others's use
        _loop_spec_begin
            fvel_vfieldk(ix,iy,iz,1:3) = fvel_vfieldk(ix,iy,iz,1:3) &
                                    + (- laplacian(ix,iy,iz) / Re) & 
                                    * vel_vfieldk(ix,iy,iz,1:3)
        _loop_spec_end
        
        invdt = 1.0_dp/dt
        
        ! Prediction: u(n+1)_1 = ((1/dt + (1 - implicitness) L) u(n) + N(n)) /
        !                        (1/dt - implicitness L)

        _loop_spec_begin
            timestep_prefieldk(ix,iy,iz,1:3) = &
                ((invdt + (1 - implicitness) * (- laplacian(ix,iy,iz) / Re)) &
                * vel_vfieldk(ix,iy,iz,1:3) + timestep_nonlinterm_prev(ix,iy,iz,1:3)) & !nonlintrm
                / (invdt - implicitness * (- laplacian(ix,iy,iz) / Re))
        _loop_spec_end
        
        ! Corrector iterations
        do c = 1, ncorr

            ! Compute the prediction norm:
            call vfield_norm(timestep_prefieldk, norm, .true.)
            
            call fftw_vk2x(timestep_prefieldk, timestep_prefieldxx)
            call rhs_nonlin_term(timestep_prefieldxx, timestep_prefieldk, timestep_nonlinterm_next)

            ! Now we have N(n+1)^c in state(:,:,:,1:3)
            _loop_spec_begin
                timestep_corfieldk(ix,iy,iz,1:3) = implicitness &
                * (timestep_nonlinterm_next(ix,iy,iz,1:3) - timestep_nonlinterm_prev(ix,iy,iz,1:3)) &
                / (invdt - implicitness * (-laplacian(ix,iy,iz) / Re))
                ! Update the prediction
                timestep_prefieldk(ix,iy,iz,1:3) = timestep_prefieldk(ix,iy,iz,1:3) &
                                                 + timestep_corfieldk(ix,iy,iz,1:3)
            _loop_spec_end

            ! Compute the corrector norm:
            call vfield_norm(timestep_corfieldk, error, .true.)

            ! update
            timestep_nonlinterm_prev(:,:,:,1:3) = timestep_nonlinterm_next(:,:,:,1:3)

            if (error/norm < steptol) then
                ! accept step
                vel_vfieldk(:,:,:,1:3) = timestep_prefieldk(:,:,:,1:3)

                ! apply pressure, galilean invariance and symmetry projections
                call vfield_pressure(vel_vfieldk)
                call vfield_galinv(vel_vfieldk)
                call symmopps_project_all(vel_vfieldk,vel_vfieldk)
                call vfield_solvediv(vel_vfieldk)      
    
                ! update the physical space version
                call fftw_vk2x(vel_vfieldk, vel_vfieldxx)

                ! log the step
                ncorr_last = c
                error_last = error/norm
                exit
            end if
            
        end do
        
        if (error/norm > steptol) then 
            write(out, *) 'timestep: Timestep did not converge.'
            write(out, *) 'timestep: Error = ', error_last
            write(out, *) 'timestep: Terminating at t = ', time
            flush(out)
            error stop            
        end if
        
    end subroutine timestep_precorr

!==============================================================================

    subroutine timestep_set_dt
        
        dt = (courant_target / courant) * dt
        
        if (dtmax > 0 .and. dt > dtmax) then
            dt = dtmax
            courant = courant * dtmax / dt
        else
            courant = courant_target
        end if
        
    end subroutine timestep_set_dt

!==============================================================================

    subroutine timestep_courant(vfieldxx)
        real(dp), intent(in) :: vfieldxx(:, :, :, :)
        real(dp) :: sfieldxx(nyy, nzz_perproc, nxx)
        real(dp) :: my_courant

        sfieldxx(:, :, :) = abs(vfieldxx(:,:,:,1)) * dt / dx &
                                + abs(vfieldxx(:,:,:,2)) * dt / dy &
                                + abs(vfieldxx(:,:,:,3)) * dt / dz
        my_courant = maxval(sfieldxx(:, :, :))
        if (.not. adaptive_dt) then 
            call MPI_REDUCE(&
            my_courant,courant,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,mpi_err)
        else
            call MPI_ALLREDUCE(&
            my_courant,courant,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpi_err)
        end if

    end subroutine timestep_courant

!==============================================================================

    subroutine timestep_write
        
        ! outputting statistics
        
        if (my_id==0) then
            
           ! outputting all this in the stat file

            inquire(file=TRIM(steps_file), exist=there, opened=there2)
            if (.not.there) then
            open(newunit=steps_ch,file=TRIM(steps_file),form='formatted')
                write(steps_ch,"(A2,"//i4_len//","//"4"//sp_len//","//i4_len//")") &
                    "# ", "itime", "time", "dt", "courant", "err_corr", "ncorr"
            end if
            if(there.and..not.there2) then
            open(newunit=steps_ch,file=TRIM(steps_file),position='append')
            end if
            write(steps_ch,"(A2,"//i4_f//","//"4"//sp_f//","//i4_f//")")&
                "  ", itime, time, dt, courant, error_last, ncorr_last

           steps_written = .true.

        end if

    end subroutine timestep_write

!==============================================================================

end module timestep

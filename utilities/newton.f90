#include "../macros.h"
module mnewton
    use numbers
    use openmpi
    use io
    use run
    use solver
    use symmpos

    complex(dpc), allocatable, dimension(:, :, :, :) :: &
    part_x_vfieldk, part_z_vfieldk

    contains

    include "NewtonHook.f90"
    include "GMRESm.f90"

    subroutine getrhs(x, y)

        ! function to be minimised

        real(dp), intent(in)  :: x(nnewt)
        real(dp), intent(out) :: y(nnewt)
        real(dp) :: y_(nnewt)
        integer(i4) :: ims, ims_, i_delta_t

        real(dp) :: newton_shift_x, newton_shift_y

        do ims = 0, ms -1
            ims_ = modulo(ims+1,ms)
            
            if (my_id == 0 .and. ims == 0) then

                if (find_shift_x) then 
                    newton_shift_x = x(i_find_period + i_find_shift_x) * scale_dex 
                end if
                
                if (find_shift_z) then 
                    newton_shift_z = &
                        x(i_find_period + i_find_shift_x + 1) * scale_dez 
                end if

            end if

            call solver_steporbit(ndtss(ims+1), ims, x)
            
            if (ims == ms - 1 .and. find_shift_x) then 
                call MPI_BCAST(newton_shift_x, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
                call symmops_shiftx(- newton_shift_x, vel_vfieldk_now, vel_vfieldk_now)
            end if

            if (ims == ms - 1 .and. find_shift_z) then 
                call MPI_BCAST(newton_shift_z, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
                call symmops_shiftz(- newton_shift_z, vel_vfieldk_now, vel_vfieldk_now)
            end if

            if (ims == ms - 1 .and. (find_shift_x .or. find_shift_z) ) then 
                call fftw_vk2x(vel_vfieldk_now, vel_vfieldxx_now)
            end if

            if (ims == ms -1) call solver_relative_symmetries_apply(vel_vfieldxx_now, vel_vfieldk_now)
            call solver_vectorize(vel_vfieldk_now, ims_, y_)
        end do
        y = y_ - x            ! diff
        
        if (my_id == 0) then

            if (find_period) then
                
                do ims=0, ms-1

                    i_delta_t = ims * nnewt_pershoot + 1

                    if (ims /= 0) then 
                        i_delta_t = i_delta_t + i_find_shift_x + i_find_shift_z
                    end if
                    
                    y(i_delta_t) = 0 ! constraints, rhs=0

                end do

            end if

            if (find_shift_x) then

                y(i_find_period + i_find_shift_x) = 0
            
            end if

            if (find_shift_z) then 

                y(i_find_period + i_find_shift_x + i_find_shift_z) = 0

            end if

        end if
        
    end subroutine getrhs

!==============================================================================

    subroutine multJ(x, y)
        
        !  Jacobian of function + lhs of constraints on update
        
        real(dp), intent(in)     :: x(nnewt)
        real(dp), intent(out)    :: y(nnewt)   
        real(dp) :: eps, s(nnewt), d, dt_new, dt_old
        integer(i4)  :: ims, i_delta_t, k
                        ! (F(x0+eps.x)-F(x0))/eps
        
        write(out, *) 'Jacobian call '
        
        eps = sqrt(solver_dotprod_ms(x,x))
        
        write(out, *) 'eps = ', eps
        if(abs(eps) < small)  then
            write(out,*) 'multJ: eps=0 (1)'
            flush(out)
            stop
        end if  
        eps = epsJ * sqrt(solver_dotprod_ms(new_x,new_x)) / eps
        write(out, *) 'eps_ = ', eps
        if(abs(eps) < small)  then 
            write(out,*) 'multJ: eps=0 (2)'
            flush(out)
            stop
        end if
        y = new_x + eps*x
        ! This is for debugging
        d = solver_dotprod_ms(y, y)
        write(out, *) 'Jacobian: normy2 = ', d
        
        call getrhs(y, s)
        call getrhs(new_x, new_fx)
        y = (s - new_fx) / eps

        do ims = 0, ms -1
            if (my_id == 0 .and. find_period) then 
                
                i_delta_t = ims * nnewt_pershoot + 1

                if (ims /= 0) then 
                    i_delta_t = i_delta_t + i_find_shift_x + i_find_shift_z
                end if

                dt_new = new_x(i_delta_t) * scaleT / ndtss(ims + 1)

            end if
            
            if (find_period) then
                call MPI_BCAST(dt_new, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
                write(out, *) 'MultJ: dt_new = ', dt_new

                dt_old = dt
                dt = dt_new
            end if

            ! contstraint, 
            ! no update in trajectory direction

            call solver_steporbit(1,ims,new_x)

            call solver_vectorize(vel_vfieldk_now, ims, s)
            s(ims*nnewt_pershot+1:(ims+1)*nnewt_pershot) = &
                (s(ims*nnewt_pershot+1:(ims+1)*nnewt_pershot) - new_x(ims*nnewt_pershot+1:(ims+1)*nnewt_pershot)) / dt
            d = solver_dotprod(ims,s,x)

            if (find_period) dt = dt_old
            
            if (my_id == 0 .and. i_find_period > 0) then
                if (find_period) y(i_delta_t) = d
            end if

            ! contstraint, 
            ! no update in x-direction
            
            if (ims == 0) then

                call solver_tensorize(vel_vfieldk_now, ims, new_x)
                do k = 1, 3
                    diffops_partx(vel_vfieldk_now(:, :, :, k), part_x_vfieldk(:, :, :, k))
                    diffops_partz(vel_vfieldk_now(:, :, :, k), part_z_vfieldk(:, :, :, k))
                end do

                solver_vectorize(part_x_vfieldk, ims, s)
                d = solver_dotprod(ims, s, x)

                if (my_id == 0 .and. find_shift_x) then
                    y(i_find_period + 1) = d
                end if

                solver_vectorize(part_z_vfieldk, ims, s)
                d = solver_dotprod(ims, s, x)

                if (my_id == 0 .and. find_shift_z) then
                    y(i_find_period + i_find_shift_x + 1) = d
                end if

            end if

        end do
        
    end subroutine multJ

!==============================================================================

    subroutine saveorbit
        
        real(dp) :: norm_x
        integer(i4) :: KILL_SWITCH_PERIOD = 0,un, ims, i_delta_t
        character*1 :: ims_str
        complex(dpc) :: vfieldk(nx_perproc, ny_half, nz, 3)

        ndts = sum(ndtss)
        norm_x = sqrt(solver_dotprod_ms(new_x,new_x))
        
        if (my_id == 0) then
            open(newunit=un,status='unknown',access='append',file='newton.dat')
            if(new_nits==0) then
                write(un,"(A2,"//"4"//i4_f//")") "# ", ndts, mgmres, nnewt, ms
                write(un, "(A2,"//"2"//i4_len//","//"5"//sp_len//")") "# ", "nits", "gits", &
                                            "rel_err", "tol_ratio", "del", "tol", "norm_x"
            end if
            write(un, "(A2,"//"2"//i4_f//","//"5"//sp_f//")") "  ", new_nits, new_gits, &
                                new_tol / norm_x, new_tol / tol, new_del, new_tol, norm_x
            close(un)
        end if
        
        if (my_id == 0 .and. (i_find_period + i_find_shift_x + i_find_shift_z) > 0) then
            if (find_period) then
                period = 0
                do ims = 0, ms -1

                    i_delta_t = ims * nnewt_pershoot + 1
                    if (ims /= 0) then 
                        i_delta_t = i_delta_t + i_find_shift_x + i_find_shift_z
                    end if        

                    period = period + new_x(i_delta_t) * scaleT
                    ! Kill switch for negative period guesses
                    if (new_x(i_delta_t) * scaleT < 0) then
                        KILL_SWITCH_PERIOD = 1
                    end if
                end do
            end if

            open(newunit=un,status='unknown',access='append',file='guesses.dat')
            if (ms > 1) then
                write(un,"(A2,"//dp_f//","//i4_f//")") "# ", period, ndts
                do ims = 0, ms -1
                    write(un, "("//i4_f//","//dp_f//","//i4_f//")") &
                            new_nits, new_x(ims*nnewt_pershot+1) * scaleT, ndtss(ims+1)
                end do
            else
                if(new_nits==0) write(un,"(A2,"//i4_f//","//i4_f//","//i4_f//")") "# ", ndts
                write(un, "("//i4_f//","//dp_f//","//dp_f//")") new_nits, period
            end if
            close(un)

            if (find_shift_x) then 
                shift_x = new_x(i_find_period + 1) * scale_dex
            end if

            if (find_shift_x .and. find_shift_z) then
                shift_z = new_x(i_find_period + 2) * scale_dez 
            else if (find_shift_z) then
                shift_z = new_x(i_find_period + 1) * scale_dez
            end if

            if (find_shift_x .or. find_shift_z) then

                open(newunit=un,status='unknown',access='append',file='shifts.dat')
                write(un,"(A2,"//dp_f//","//dp_f//")") "# ", shift_x, shift_z
                close(un)
            
            end if

        end if

        ! Broadcast the KILL_SWITCH status
        call MPI_BCAST(KILL_SWITCH_PERIOD, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)

        if (KILL_SWITCH_PERIOD == 1) then
            write(out, *) "Period guess is negative, stopping."
            call newton_signal_not_converged
            call run_exit
        end if
        
        ! quit run if already converged and do not start again
        if (new_nits == 0 .and. new_tol / tol < 1) then
            call newton_done
            call run_exit
        end if

        ! Save the state file:
        write(file_ext, "(i6.6)") new_nits
        do ims = 0, ms - 1
            write(ims_str, "(i1.1)") ims
            if (ms > 1) then
                fname = 'newton.'//file_ext//'-'//ims_str
            else
                fname = 'newton.'//file_ext
            end if
            call solver_tensorize(vfieldk, ims, new_x)
            call fieldio_write(vfieldk)
        end do

    end subroutine saveorbit

!==============================================================================

    subroutine newton_done
        
        integer(i4) :: un
        if (my_id == 0) then
            open(newunit=un,file='NEWTON_DONE',position='append')
            write(un,*)
            close(un)
        end if
    end subroutine newton_done

!==============================================================================

    subroutine newton_signal_converged
        
        integer(i4) :: un
        if (my_id == 0) then
            open(newunit=un,file='NEWTON_CONVERGED',position='append')
            write(un,*)
            close(un)
        end if
    end subroutine newton_signal_converged

!==============================================================================

    subroutine newton_signal_not_converged
        
        integer(i4) :: un
        if (my_id == 0) then
            open(newunit=un,file='NEWTON_NOT_CONVERGED',position='append')
            write(un,*)
            close(un)
        end if
    end subroutine newton_signal_not_converged

end module mnewton

!==============================================================================

program newton
    !*************************************************************************
    !  Example with one extra parameter T.
    !  Can put parameter and constraint = 0 when not required.
    !- - - - - - - - - - - - - - -
    !  Newton vector:
    !    x(1)   = T / scaleT
    !    x(2:) = state vector x
    !
    !  Extra constraint:
    !    (F(x)-x). dx/dt = 0 .
    !       no update along direction of trajectory.
    !
    !  Jacobian approximation:
    !    dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps
    !
    !*************************************************************************
    
    ! modules:
    use numbers
    use openmpi
    use io
    use parameters
    use fftw
    use symmops
    use fieldio
    use vfield
    use timestep
    use run
    use solver
    use mnewton
    
    real(dp)           :: d
    logical            :: forbidNewton
    
    integer(i4) :: ims, info, i_delta_t
    logical     :: fexist
    character*1 :: ims_str

    call run_init
    IC = 0
    adaptive_dt = .false.
    
    inquire(file='NEWTON_NOT_CONVERGED', exist=forbidNewton)
    if (forbidNewton) call run_exit

    inquire(file='NEWTON_DONE', exist=forbidNewton)
    if (forbidNewton) call run_exit

    call solver_data_read
    call solver_set_problem_size
    
    ! allocate vectors:
    allocate(new_x(nnewt))
    allocate(new_fx(nnewt))
    new_x = 0

    allocate(part_x_vfieldk(nx_perproc, ny_half, nz, 3))
    allocate(part_z_vfieldk(nx_perproc, ny_half, nz, 3))
    
    do ims = 1, ms
        if (ndtss(ims) < 0) ndtss(ims) = nint(periods(ims)/dt)
    end do
    
    if (my_id == 0 .and. nscalars > 0) then
        if (find_period) then
            do ims = 0, ms - 1

                i_delta_t = ims * nnewt_pershoot + 1

                if (ims /= 0) then 
                    i_delta_t = i_delta_t + i_find_shift_x + i_find_shift_z
                end if
    
                new_x(i_delta_t) = periods(ims + 1)
            end do
        end if
        
        if (find_shift_x) then
            new_x(i_find_period + 1) = shift_x
        end if
        
        if (find_shift_z) then 
            new_x(i_find_period + i_find_shift_x + 1) = shift_z
        end if

    end if
    
    ! Load the state    
    ! Initial time
    if (ms > 1) then
        write(file_ext, "(i6.6)") 0
        fname = 'state.'//file_ext//'-0'
        inquire(file=fname, exist=fexist)
        if (fexist) then
            ! read from disk
            do ims = 0, ms - 1
                write(ims_str, "(i1.1)") ims
                fname = 'state.'//file_ext//'-'//ims_str
                call fieldio_read(vel_vfieldk_now)
                call solver_vectorize(vel_vfieldk_now, ims, new_x)
            end do
        else
            ! construct from time zero
            write(file_ext, "(i6.6)") 0
            fname = 'state.'//file_ext
            call fieldio_read(vel_vfieldk_now)
            call solver_vectorize(vel_vfieldk_now, 0, new_x)
            scaleT = 1.0_dp
            do ims = 1, ms - 1
                call solver_steporbit(ndtss(ims), ims-1, new_x)
                call solver_vectorize(vel_vfieldk_now, ims, new_x)
            end do
        end if
    else
        write(file_ext, "(i6.6)") 0
        fname = 'state.'//file_ext
        call fieldio_read(vel_vfieldk_now)
        call solver_vectorize(vel_vfieldk_now, 0, new_x)
    end if
    time = 0
    itime = 0

    ! Set the scales
    d = solver_dotprod(0,new_x,new_x)
    if(abs(d) < small) d=1.0_dp
    scaleT = periods(1) / sqrt(d)
    scale_dex = shift_x / sqrt(d)
    scale_dez = shift_z / sqrt(d)

    if (my_id == 0 .and. nscalars > 0) then
        if (find_period) then
            do ims  = 0, ms - 1

                i_delta_t = ims * nnewt_pershoot + 1

                if (ims /= 0) then 
                    i_delta_t = i_delta_t + i_find_shift_x + i_find_shift_z
                end if

                new_x(i_delta_t) = new_x(i_delta_t) / scaleT
            end do
        end if

        if (find_shift_x) then 
            new_x(i_find_period + 1) = new_x(i_find_period + 1) / scale_dex 
        end if

        if (find_shift_z) then 
            new_x(i_find_period + i_find_shift_x + 1) = &
                new_x(i_find_period + i_find_shift_x + 1) / scale_dez 
        end if

    end if

    tol  = rel_err * sqrt(d)
    del  = del     * sqrt(d)
    mndl = mndl    * sqrt(d)
    mxdl = mxdl    * sqrt(d)
    
    info = 1
    call newtonhook(getrhs, multJ, saveorbit, solver_dotprod_ms, &
                    mgmres, nnewt, gtol, tol, del, mndl, mxdl, nits, info, out)

    if (info == 0) then
        call newton_signal_converged
    else
        call newton_signal_not_converged
    end if
     
    call run_exit   
    
end program newton

#include "macros.h"
#include "mt19937_64.f90"
module vfield
    use numbers
    use openmpi
    use io
    use parameters
    use fieldio
    use fftw
    use diffops
    use symmops
    
    logical     :: seed_init = .false.
    integer(i8) :: seed
    complex(dpc), allocatable, dimension(:, :, :, :) :: laminar_vfieldk ! Laminar solution
    
    contains 

!==============================================================================

    subroutine vfield_init

        allocate(laminar_vfieldk(nx_perproc, ny_half, nz, 3))
        call vfield_laminar(laminar_vfieldk)

    end subroutine vfield_init

!==============================================================================

    subroutine vfield_inprod(vfieldk1, vfieldk2, res, allreduce)
        ! res = <vfieldk1, vfieldk2>, returns only the real part
        ! of the result
        complex(dpc), intent(in) :: vfieldk1(:, :, :, :)
        complex(dpc), intent(in) :: vfieldk2(:, :, :, :)
        real(dp), intent(out) :: res
        logical, intent(in) :: allreduce
        real(dp)          :: res1
        complex(dpc)      :: dummy
        _indices
        integer(i4) :: n

        res1 = 0
        do n=1,3
            _loop_spec_begin
                dummy = conjg(vfieldk1(ix,iy,iz,n)) * vfieldk2(ix,iy,iz,n)
                ! correct for double counting
                if (iy == 1) dummy = dummy * 0.5_dp
                res1 = res1 + dummy%re
            _loop_spec_end
        end do

        if (.not. allreduce) then
            call MPI_REDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, 0, &
            MPI_COMM_WORLD, mpi_err)
        else
            call MPI_ALLREDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, &
            MPI_COMM_WORLD, mpi_err)
        end if
        
    end subroutine vfield_inprod

!==============================================================================

    subroutine vfield_enstrophy(vfieldk, res, allreduce)

        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        real(dp), intent(out) :: res
        logical, intent(in) :: allreduce
        real(dp)          :: res1
        complex(dpc)      :: dummy
        _indices
        integer(i4) :: n

        res1 = 0
        do n=1,3
            _loop_spec_begin
                dummy = laplacian(ix,iy,iz) * &
                            conjg(vfieldk(ix,iy,iz,n)) * vfieldk(ix,iy,iz,n)
                ! correct for double counting
                if (iy == 1) dummy = dummy * 0.5_dp
                res1 = res1 + dummy%re
            _loop_spec_end
        end do

        if (.not. allreduce) then
            call MPI_REDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, 0, &
            MPI_COMM_WORLD, mpi_err)
        else
            call MPI_ALLREDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, &
            MPI_COMM_WORLD, mpi_err)
        end if
        
    end subroutine vfield_enstrophy

!==============================================================================

    subroutine vfield_norm2_horizontal(vfieldk, res, allreduce)

        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        real(dp), intent(out) :: res
        logical, intent(in) :: allreduce
        real(dp)          :: res1
        complex(dpc)      :: dummy
        _indices

        res1 = 0
        _loop_spec_begin
            dummy = conjg(vfieldk(ix,iy,iz,1)) * vfieldk(ix,iy,iz,1) + &
                    conjg(vfieldk(ix,iy,iz,3)) * vfieldk(ix,iy,iz,3)
            ! correct for double counting
            if (iy == 1) dummy = dummy * 0.5_dp
            res1 = res1 + dummy%re
        _loop_spec_end

        if (.not. allreduce) then
            call MPI_REDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, 0, &
            MPI_COMM_WORLD, mpi_err)
        else
            call MPI_ALLREDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, &
            MPI_COMM_WORLD, mpi_err)
        end if
        
    end subroutine vfield_norm2_horizontal

!==============================================================================

    subroutine vfield_norm2(vfieldk, res, allreduce)
        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        logical, intent(in) :: allreduce
        real(dp), intent(out) :: res

        call vfield_inprod(vfieldk, vfieldk, res, allreduce)
        
    end subroutine vfield_norm2

!==============================================================================

    subroutine vfield_norm(vfieldk, res, allreduce)
        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        logical, intent(in) :: allreduce
        real(dp), intent(out) :: res
        
        call vfield_norm2(vfieldk, res, allreduce)
        res = sqrt(res)
        
    end subroutine vfield_norm

!==============================================================================

    subroutine vfield_random(vfieldk, symmetric)
        use mt19937_64
        complex(dpc), intent(out) :: vfieldk(:, :, :, :)
        logical, intent(in) :: symmetric
        real(dp) :: random_vfieldxx(nyy, nzz_perproc, nxx, 3)
        
        _indices
        _indicess

        integer(i4) :: un, istat, n
        real(dp)    :: dummy, e_tot        

        ! If the parameters.in contains a seed number, use it
        if (random_seed /= -1) then
            seed = random_seed
        else if (my_id==0) then
            ! otherwise read from entropy
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
               seed_init = .false.
            else
                write(out,*) "m_runs: /dev/urandom not found, stopping."
                flush(out)
                error stop
            end if
        end if

        ! Broadcast the seed
        call MPI_BCAST(seed,1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpi_err)
        
        ! initialize the seed
        if (.not. seed_init) then
            call init_genrand64(seed)
            seed_init = .true.
        end if

        ! bring the processors to their own places in the random sequence
        do n = 1,my_id*nyy*nzz_perproc*nxx*3
            dummy = genrand64_real1()
        end do
        
        ! fill the arrays with [-1, 1]
        do n=1,3;
            _loop_phys_begin
                random_vfieldxx(iyy,izz,ixx,n) = (2 * genrand64_real1() - 1)
            _loop_phys_end
        end do;        

        call fftw_vx2k(random_vfieldxx, vfieldk)

        ! scale the arrays with e^{-sk^2}
        do n=1,3;
            _loop_spec_begin
                vfieldk(ix,iy,iz,n) = vfieldk(ix,iy,iz,n) &
                                   * exp(-random_smooth * laplacian(ix,iy,iz))
            _loop_spec_end
        end do;

        call vfield_pressure(vfieldk)
        call vfield_galinv(vfieldk)

        if (symmetric) call symmopps_project_all(vfieldk,vfieldk)
        ! Normalize to a specific energy
        call vfield_norm2(vfieldk, e_tot, .true.)
        vfieldk(:,:,:,1:3) = sqrt(random_energy * ekin_lam / e_tot) * vfieldk(:,:,:,1:3)
        
        call vfield_solvediv(vfieldk)

        ! once more, doesn't hurt
        call vfield_pressure(vfieldk)
        call vfield_galinv(vfieldk)
        if (symmetric) call symmopps_project_all(vfieldk,vfieldk)
        call vfield_solvediv(vfieldk)

    end subroutine vfield_random
    
!==============================================================================

    subroutine vfield_laminar(vfieldk)
        complex(dpc), intent(out) :: vfieldk(:, :, :, :)
        
        vfieldk(:,:,:,1:3) = 0
        if (ix_zero /= -1 .and. iy_force /= -1) then
            if (tilting) then
                if (forcing == 1) then
                    vfieldk(ix_zero,iy_force,1,1)%im = -cos(tilt_angle * PI / 180.0_dp) * 0.5_dp
                    vfieldk(ix_zero,iy_force,1,3)%im = -sin(tilt_angle * PI / 180.0_dp) * 0.5_dp
                elseif (forcing == 2) then
                    vfieldk(ix_zero,iy_force,1,1)%re = cos(tilt_angle * PI / 180.0_dp) * 0.5_dp
                    vfieldk(ix_zero,iy_force,1,3)%re = sin(tilt_angle * PI / 180.0_dp) * 0.5_dp
                endif
            else
                if (forcing == 1) then
                    vfieldk(ix_zero,iy_force,1,1)%im = -0.5_dp
                elseif (forcing == 2) then
                    vfieldk(ix_zero,iy_force,1,1)%re = 0.5_dp
                endif
            end if
        endif
    end subroutine vfield_laminar

!==============================================================================

    subroutine vfield_shapiro(t, vfieldk)
        real(dp), intent(in) :: t
        complex(dpc), intent(out) :: vfieldk(:, :, :, :)
        real(dp) :: shapiro_vfieldxx(nyy, nzz_perproc, nxx, 3)
        real(dp) :: decay
        _indicess

        decay = exp(-shapiro_lambda**2 * t / Re)

        _loop_phys_begin
            shapiro_vfieldxx(iyy,izz,ixx,1) = -(decay / (shapiro_alpha**2 + shapiro_beta**2)) * &
                (shapiro_lambda * shapiro_beta * &
                    dcos(shapiro_alpha * vfield_coordinatexx(iyy,izz,ixx,1)) * &
                    dsin(shapiro_beta * vfield_coordinatexx(iyy,izz,ixx,2)) * &
                    dsin(shapiro_gamma * vfield_coordinatexx(iyy,izz,ixx,3)) &
                 + shapiro_gamma * shapiro_alpha * &
                    dsin(shapiro_alpha * vfield_coordinatexx(iyy,izz,ixx,1)) * &
                    dcos(shapiro_beta * vfield_coordinatexx(iyy,izz,ixx,2)) * &
                    dcos(shapiro_gamma * vfield_coordinatexx(iyy,izz,ixx,3)) &
                )

            shapiro_vfieldxx(iyy,izz,ixx,2) = (decay / (shapiro_alpha**2 + shapiro_beta**2)) * &
                (shapiro_lambda * shapiro_alpha * &
                    dsin(shapiro_alpha * vfield_coordinatexx(iyy,izz,ixx,1)) * &
                    dcos(shapiro_beta * vfield_coordinatexx(iyy,izz,ixx,2)) * &
                    dsin(shapiro_gamma * vfield_coordinatexx(iyy,izz,ixx,3)) &
                    - shapiro_gamma * shapiro_beta * &
                    dcos(shapiro_alpha * vfield_coordinatexx(iyy,izz,ixx,1)) * &
                    dsin(shapiro_beta * vfield_coordinatexx(iyy,izz,ixx,2)) * &
                    dcos(shapiro_gamma * vfield_coordinatexx(iyy,izz,ixx,3)) &
                )

            shapiro_vfieldxx(iyy,izz,ixx,3) = decay * &
                    dcos(shapiro_alpha * vfield_coordinatexx(iyy,izz,ixx,1)) * &
                    dcos(shapiro_beta * vfield_coordinatexx(iyy,izz,ixx,2)) * &
                    dsin(shapiro_gamma * vfield_coordinatexx(iyy,izz,ixx,3))
                
        _loop_phys_end

        call fftw_vx2k(shapiro_vfieldxx, vfieldk)
        

    end subroutine vfield_shapiro

!==============================================================================

    subroutine vfield_pressure(vfieldk)
        ! Apply pressure solver on vfieldk(:,:,:,1:3) 
        complex(dpc), intent(inout) :: vfieldk(:, :, :, :)
        _indices
        integer(i4)  :: d
        complex(dpc) :: dpressure
        
        _loop_spec_begin
            dpressure = 0
            do d=1,3
                dpressure = dpressure + vfieldk(ix,iy,iz,d) * &
                                        vfield_coordinatek(ix,iy,iz,d) * &
                                        inverse_laplacian(ix,iy,iz)
            end do
            
            do d=1,3
                vfieldk(ix,iy,iz,d) = vfieldk(ix,iy,iz,d) - dpressure * &
                                                vfield_coordinatek(ix,iy,iz,d)
            end do
        _loop_spec_end
                
    end subroutine vfield_pressure

!==============================================================================

    subroutine vfield_galinv(vfieldk)
        ! Set (0,0,0)-mode to 0 (Galilean invariance)
        complex(dpc), intent(inout) :: vfieldk(:, :, :, :)

        ! index ix_zero has the zero mode in x
        if(ix_zero /= -1) vfieldk(ix_zero,1,1,1:3) = 0
    end subroutine vfield_galinv

!==============================================================================

    subroutine vfield_solvediv(vfieldk)
        complex(dpc), intent(inout) :: vfieldk(:, :, :, :)

        _indices

        ! uses
        !   u k_x + v k_y + w k_z = 0
        ! to solve for one velocity component in terms of the other two

        if (nx - 1 >= ny_half .and. nx - 1 >= nz - 1) then
            _loop_spec_begin
                if (ix_zero /= -1 .and. ix==ix_zero) cycle
                vfieldk(ix,iy,iz,1) = -(vfieldk(ix,iy,iz,2)*ky(iy) &
                                            + vfieldk(ix,iy,iz,3)*kz(iz)) / kx(ix)
            _loop_spec_end
        elseif (ny_half >= nx - 1 .and. ny_half >= nz - 1) then
            _loop_spec_begin
                if (iy==1) cycle
                vfieldk(ix,iy,iz,2) = -(vfieldk(ix,iy,iz,1)*kx(ix) &
                                            + vfieldk(ix,iy,iz,3)*kz(iz)) / ky(iy)
            _loop_spec_end
        else
            _loop_spec_begin
                if (iz==1) cycle
                vfieldk(ix,iy,iz,3) = -(vfieldk(ix,iy,iz,1)*kx(ix) &
                                                + vfieldk(ix,iy,iz,2)*ky(iy)) / kz(iz)
            _loop_spec_end
        end if

    end subroutine vfield_solvediv

!==============================================================================

    subroutine vfield_truncate(k_c, vfieldk)
        ! low-pass filter setting fluctuations with spatial frequencies higher 
        ! than k_c to 0
        real(dp), intent(in) :: k_c 
        complex(dpc), intent(inout) :: vfieldk(:, :, :, :)
        
        real(dp) :: normk 
        
        _indices

        _loop_spec_begin

            normk = sqrt(kx(ix) ** 2 + ky(iy) ** 2 + kz(iz) ** 2)

            if (normk > k_c) then 
                
                vfieldk(ix, iy, iz, 1:3) = 0

            end if

        _loop_spec_end
    
    end subroutine vfield_truncate

!==============================================================================
    real(dp) function vfield_coordinatek(ix,iy,iz,d)
        integer(i4), intent(in) :: ix, iy, iz, d
        select case (d)
            case(1)
                vfield_coordinatek = kx(ix)
            case(2)
                vfield_coordinatek = ky(iy)
            case(3)
                vfield_coordinatek = kz(iz)
        end select
    end function vfield_coordinatek

!==============================================================================

    real(dp) function vfield_coordinatexx(iyy,izz,ixx,d)
        integer(i4), intent(in) :: ixx, iyy, izz, d
        
        select case (d)
            case(1)
                vfield_coordinatexx = dx * (ixx - 1)
            case(2)
                vfield_coordinatexx = dy * (iyy - 1)
            case(3)
                vfield_coordinatexx = dz * (nzz_perproc * my_id + izz - 1)
        end select
    end function vfield_coordinatexx

!==============================================================================

end module vfield

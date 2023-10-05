#include "macros.h"
module rhs
    use numbers
    use openmpi
    use io
    use parameters
    use fieldio 
    use fftw
    use vfield

    contains

!==============================================================================

    subroutine rhs_nonlin_term(vel_vfieldxx, vel_vfieldk, fvel_vfieldk)

        real(dp), intent(in)  :: vel_vfieldxx(:, :, :, :)
        complex(dpc), intent(in) :: vel_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: fvel_vfieldk(:, :, :, :)

        complex(dpc) :: rhs_vfieldk(nx_perproc, ny_half, nz, 6), advect(3), div
        real(dp) :: rhs_vfieldxx(nyy, nzz_perproc, nxx, 5), u_out_u(6), trace

        integer(i4) :: n, i, j

        ! LES fields
        real(dp), allocatable :: &
            rhs_nu_T_fieldxx(:, :, :), &
            rhs_rofs_tensorxx(:, :, :, :)
        complex(dpc), allocatable :: &
            rhs_model_vfieldk(:, :, :, :), &
            rhs_div_model_vfieldk(:, :, :, :)

        _indices
        _indicess

        if (LES) then
            if (.not. allocated(rhs_nu_T_fieldxx)) &
                allocate(rhs_nu_T_fieldxx(nyy, nzz_perproc, nxx))
            if (.not. allocated(rhs_rofs_tensorxx)) &
                allocate(rhs_rofs_tensorxx(nyy, nzz_perproc, nxx, 6))
            if (.not. allocated(rhs_model_vfieldk)) &
                allocate(rhs_model_vfieldk(nx_perproc, ny_half, nz, 3))
            if (.not. allocated(rhs_div_model_vfieldk)) &
                allocate(rhs_div_model_vfieldk(nx_perproc, ny_half, nz, 3))

        end if
        
        ! get 6 velocity products:
        _loop_phys_begin

            do n = 1, 6
                u_out_u(n) = vel_vfieldxx(iyy,izz,ixx,isym(n)) &
                                * vel_vfieldxx(iyy,izz,ixx,jsym(n))
            end do

            ! Basdevant 1983, subtract the trace
            trace = u_out_u(nsym(1,1)) + u_out_u(nsym(2,2)) &
                                                + u_out_u(nsym(3,3))
            u_out_u(nsym(1,1)) = u_out_u(nsym(1,1)) - trace / 3.0_dp
            u_out_u(nsym(2,2)) = u_out_u(nsym(2,2)) - trace / 3.0_dp
            ! No need to do
            !   u_out_u(nsym(3,3)) = u_out_u(nsym(3,3)) - trace / 3.0_dp
            ! as it doesn't get used.
            
            ! No need to write to rhs_vfieldxx(:,:,:,6) as it is not used.
            do n = 1, 5
                rhs_vfieldxx(iyy,izz,ixx,n) = u_out_u(n)
            end do
        _loop_phys_end

        do n = 1,5
            call fftw_sx2k(rhs_vfieldxx(:,:,:,n), rhs_vfieldk(:,:,:,n))
        end do

        ! Get one element on the diagonal from the tracelessness
        rhs_vfieldk(:,:,:,nsym(3,3)) = &
            -(rhs_vfieldk(:,:,:,nsym(1,1)) + rhs_vfieldk(:,:,:,nsym(2,2)))
        
        ! Now rhs_vfieldk(:,:,:,:, 1:6) holds the products uu, uv, uw, ...
        ! in Fourier space.

        if (LES) then

            call rhs_les(vel_vfieldk, rhs_nu_T_fieldxx, rhs_rofs_tensorxx)
            
            do j = 1, 3
                do i = 1, 3
                    rhs_vfieldxx(:, :, :, i) = &
                            -2.0_dp * rhs_nu_T_fieldxx * rhs_rofs_tensorxx(:, :, :, nsym(i,j))
                end do

                call fftw_vx2k(rhs_vfieldxx, rhs_model_vfieldk)
                call diffops_div(rhs_model_vfieldk, rhs_div_model_vfieldk(:, :, :, j))

            end do
        end if    

        _loop_spec_begin                       
            ! Nonlinear term:

            ! advect(j) = - F{u_i d_i u_j}
            do j = 1, 3
                advect(j) = 0
                do i = 1 ,3
                    advect(j) = advect(j) - nabla(ix, iy, iz, i) &
                                * rhs_vfieldk(ix, iy, iz, nsym(i,j))
                end do

                if(LES) advect(j) = advect(j) - rhs_div_model_vfieldk(ix, iy, iz, j)
            end do

            ! Pressure terms
            div = 0
            do n = 1, 3
                div = div + vfield_coordinatek(ix,iy,iz,n) * advect(n)
            end do
            
            ! Extra terms to the pressure gradient due to non-divergence-free forcing terms
            if (rayleigh_friction) div = div + vfield_coordinatek(ix, iy, iz, 2) * vel_vfieldk(ix, iy, iz, 2) * sigma_R

            do n = 1, 3
                fvel_vfieldk(ix,iy,iz,n) = advect(n) &
                                    - vfield_coordinatek(ix,iy,iz,n) * div &
                                        * inverse_laplacian(ix,iy,iz)
            end do

        _loop_spec_end
        
        ! Add forcing
        if (forcing /= 0 .and. ix_zero /= -1) then
            if (tilting) then
                if (forcing == 1) then ! sine
                    fvel_vfieldk(ix_zero,iy_force,1,1) = fvel_vfieldk(ix_zero,iy_force,1,1) &
                            - imag_1 * cos(tilt_angle * PI / 180.0_dp) * 0.5_dp * amp / (4.0_dp*Re)
                    fvel_vfieldk(ix_zero,iy_force,1,3) = fvel_vfieldk(ix_zero,iy_force,1,3) &
                            - imag_1 * sin(tilt_angle * PI / 180.0_dp) * 0.5_dp * amp / (4.0_dp*Re)
                elseif (forcing == 2) then ! cosine
                    fvel_vfieldk(ix_zero,iy_force,1,1) = fvel_vfieldk(ix_zero,iy_force,1,1) &
                            + cos(tilt_angle * PI / 180.0_dp) * 0.5_dp * amp / (4.0_dp*Re)
                    fvel_vfieldk(ix_zero,iy_force,1,3) = fvel_vfieldk(ix_zero,iy_force,1,3) &
                            + sin(tilt_angle * PI / 180.0_dp) * 0.5_dp * amp / (4.0_dp*Re)
                end if
            else
                if (forcing == 1) then ! sine
                    fvel_vfieldk(ix_zero,iy_force,1,1) = fvel_vfieldk(ix_zero,iy_force,1,1) &
                            - imag_1 * 0.5_dp * amp / (4.0_dp*Re)
                elseif (forcing == 2) then ! cosine
                    fvel_vfieldk(ix_zero,iy_force,1,1) = fvel_vfieldk(ix_zero,iy_force,1,1) &
                            + 0.5_dp * amp / (4.0_dp*Re)
                end if
            end if            
        end if

        if (rayleigh_friction) then
            fvel_vfieldk(:, :, :, 1) = fvel_vfieldk(:, :, :, 1) &
                - sigma_R * (vel_vfieldk(:, :, :, 1) - laminar_vfieldk(:, :, :, 1)) 
                
            fvel_vfieldk(:, :, :, 3) = fvel_vfieldk(:, :, :, 3) &
                - sigma_R * (vel_vfieldk(:, :, :, 3) - laminar_vfieldk(:, :, :, 3))             
        end if

    end subroutine rhs_nonlin_term

!==============================================================================

    subroutine rhs_all(vel_vfieldxx, vel_vfieldk, fvel_vfieldk)
        real(dp), intent(in)  :: vel_vfieldxx(:, :, :, :)
        complex(dpc), intent(in)  :: vel_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: fvel_vfieldk(:, :, :, :)

        _indices

        call rhs_nonlin_term(vel_vfieldxx, vel_vfieldk, fvel_vfieldk)
        _loop_spec_begin
            fvel_vfieldk(ix,iy,iz,1:3) = fvel_vfieldk(ix,iy,iz,1:3) &
                                    + (- laplacian(ix,iy,iz) / Re) & 
                                    * vel_vfieldk(ix,iy,iz,1:3)
        _loop_spec_end

    end subroutine rhs_all

!==============================================================================

    subroutine rhs_les(vel_vfieldk, nu_T_fieldxx, rofs_tensorxx)
        ! compute eddy viscosity and rate of strain tensor required for LES
        
        complex(dpc), intent(in) :: vel_vfieldk(:, :, :, :)
        real(dp), intent(out) :: nu_T_fieldxx(:, :, :) ! eddy viscosity
        real(dp), intent(out) :: rofs_tensorxx(:, :, :, :) ! rate-of-strain tensor

        complex(dpc) :: rofs_tensork(nx_perproc, ny, nz, 5) 
        complex(dpc) :: rhs_temp_fieldk(nx_perproc, ny, nz)

        integer(i4) :: n
        real(dp) :: multiplier

        ! diagonal elementes of the rate of strain tensor
        call diffops_partx(vel_vfieldk(:, :, :, 1), rofs_tensork(:, :, :, nsym(1,1)))
        call diffops_party(vel_vfieldk(:, :, :, 2), rofs_tensork(:, :, :, nsym(2,2)))
        ! no need to compute (3, 3) (tracelessness)

        ! off-diagonal elements

        ! (1,2)
        call diffops_partx(vel_vfieldk(:, :, :, 2), rofs_tensork(:, :, :, nsym(1,2)))
        call diffops_party(vel_vfieldk(:, :, :, 1), rhs_temp_fieldk(:, :, :))
        rofs_tensork(:, :, :, nsym(1,2)) = &
                        (rofs_tensork(:, :, :, nsym(1,2)) + rhs_temp_fieldk) * 0.5_dp

        ! (1,3)
        call diffops_partx(vel_vfieldk(:, :, :, 3), rofs_tensork(:, :, :, nsym(1,3)))
        call diffops_partz(vel_vfieldk(:, :, :, 1), rhs_temp_fieldk(:, :, :))
        rofs_tensork(:, :, :, nsym(1,3)) = &
            (rofs_tensork(:, :, :, nsym(1,3)) + rhs_temp_fieldk) * 0.5_dp

        ! (2,3)
        call diffops_party(vel_vfieldk(:, :, :, 3), rofs_tensork(:, :, :, nsym(2,3)))
        call diffops_partz(vel_vfieldk(:, :, :, 2), rhs_temp_fieldk(:, :, :))
        rofs_tensork(:, :, :, nsym(2,3)) = &
            (rofs_tensork(:, :, :, nsym(2,3)) + rhs_temp_fieldk) * 0.5_dp

        do n = 1, 5
            call fftw_sk2x(rofs_tensork(:, :, :, n), rofs_tensorxx(:, :, :, n))
        end do

        ! get the (3 ,3) element from tracelessness
        rofs_tensorxx(:, :, :, nsym(3,3)) = &
                    -(rofs_tensorxx(:, :, :, nsym(1,1)) &
                        + rofs_tensorxx(:, :, :, nsym(2,2)))
        
        nu_T_fieldxx = 0

        do n = 1, 6
            ! correct for non-diagonal entries
            multiplier = 1.0_dp
            if (isym(n) /= jsym(n)) multiplier = 2.0_dp

            nu_T_fieldxx(:, :, :) = nu_T_fieldxx(:, :, :) &
                                    + 2.0_dp * multiplier &
                                        * rofs_tensorxx(:, :, :, n) ** 2

        end do


        nu_T_fieldxx = (smag_const * Delta_LES) ** 2 * sqrt(nu_T_fieldxx)

    end subroutine rhs_les

end module rhs

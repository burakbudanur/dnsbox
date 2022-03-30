#include "macros.h"

module symmops
    use numbers
    use openmpi
    use io
    use parameters
    use fftw

    contains 

!==============================================================================

    subroutine symmops_mirrorx(in_vfieldxx, out_vfieldxx)
        real(dp), intent(inout)  :: in_vfieldxx(:, :, :, :)
        real(dp), intent(out) :: out_vfieldxx(:, :, :, :)

        integer(i4) :: ixx
        real(dp) :: swap(nyy, nzz_perproc)

        out_vfieldxx(:,:,:,1:3) = in_vfieldxx(:,:,:,1:3)
                
        ! Reflect in x
        do ixx = 1, nxx / 2 + 1
            
            if (ixx == 1 .or. (mod(nxx, 2) == 0 .and. ixx == nxx / 2 + 1)) then
            
                ! sigma_x u(x) = - u(- x)
                out_vfieldxx(:, :, ixx, 1) = - out_vfieldxx(:, :, ixx, 1)
                        
            else             
            
                ! sigma_x u(x) = - u(- x)
                swap(:, :) = out_vfieldxx(:, :, ixx, 1)
                out_vfieldxx(:, :, ixx, 1) = - out_vfieldxx(:, :, nxx - (ixx - 2), 1)
                out_vfieldxx(:, :, nxx - (ixx - 2), 1) = -swap(:, :)

                ! sigma_x [v,w](x) = [v,w](- x)
                swap(:, :) = out_vfieldxx(:, :, ixx, 2)
                out_vfieldxx(:, :, ixx, 2) = out_vfieldxx(:, :, nxx - (ixx - 2), 2)
                out_vfieldxx(:, :, nxx - (ixx - 2), 2) = swap(:, :)

                swap(:, :) = out_vfieldxx(:, :, ixx, 3)
                out_vfieldxx(:, :, ixx, 3) = out_vfieldxx(:, :, nxx - (ixx - 2), 3)
                out_vfieldxx(:, :, nxx - (ixx - 2), 3) = swap(:, :)
            
            end if
            
        end do
        
    end subroutine symmops_mirrorx

!==============================================================================

    subroutine symmops_mirrory(in_vfieldxx, out_vfieldxx)

        real(dp), intent(inout) :: in_vfieldxx(:, :, :, :)
        real(dp), intent(out) :: out_vfieldxx(:, :, :, :)
        
        integer(i4) :: iyy
        real(dp) :: swap(nzz_perproc, nxx)

        out_vfieldxx(:,:,:,1:3) = in_vfieldxx(:,:,:,1:3)
        
       ! Reflect in y
        do iyy = 1, nyy / 2 + 1
            
            if(iyy == 1 .or. (mod(nyy, 2) == 0 .and. iyy == nyy / 2 + 1)) then
            
                ! sigma_y v(y) = - v(- y)
                out_vfieldxx(iyy, :, :, 2) = -out_vfieldxx(iyy, :, :, 2)
                        
            else             
            
                ! sigma_y v(y) = - v(- y)
                swap(:, :) = out_vfieldxx(iyy, :, :, 2)
                out_vfieldxx(iyy, :, :, 2) = -out_vfieldxx(nyy - (iyy - 2), :, :, 2)
                out_vfieldxx(nyy - (iyy - 2), :, :, 2) = -swap(:, :)

                ! sigma_y [u,w](y) = [u,w](- y)
                swap(:, :) = out_vfieldxx(iyy, :, :, 1)
                out_vfieldxx(iyy, :, :, 1) = out_vfieldxx(nyy - (iyy - 2), :, :, 1)
                out_vfieldxx(nyy - (iyy - 2), :, :, 1) = swap(:, :)

                swap(:, :) = out_vfieldxx(iyy, :, :, 3)
                out_vfieldxx(iyy, :, :, 3) = out_vfieldxx(nyy - (iyy - 2), :, :, 3)
                out_vfieldxx(nyy - (iyy - 2), :, :, 3) = swap(:, :)
            
            end if
            
        end do            
        
    end subroutine symmops_mirrory
    
!==============================================================================
    
    subroutine symmops_mirrorz(in_vfieldk, out_vfieldk)

        complex(dpc), intent(inout) :: in_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: out_vfieldk(:, :, :, :)
        complex(dpc) :: swap(nx_perproc,ny_half)
        
        integer(i4) :: iz
        
        out_vfieldk(:,:,:,1:3) = in_vfieldk(:,:,:,1:3)

        do iz = 1, nz / 2 ! iz = nz / 2 + 1 is zeroed

            if (iz == 1) then
                out_vfieldk(:, :, iz, 3) = -out_vfieldk(:, :, iz, 3)
                out_vfieldk(:, :, iz, 1:2) = out_vfieldk(:, :, iz, 1:2)

            else

                ! j-th index corresponds to (j - 1)-th positive wave number
                ! corresponding negative frequency is at N - (j - 2)
                swap(:, :) = out_vfieldk(:, :, iz, 3)
                out_vfieldk(:, :, iz, 3) = -out_vfieldk(:, :, nz - iz + 2, 3)
                out_vfieldk(:, :, nz - iz + 2, 3) = -swap(:, :)

                swap(:, :) = out_vfieldk(:, :, iz, 1)
                out_vfieldk(:, :, iz, 1) = out_vfieldk(:, :, nz - iz + 2, 1)
                out_vfieldk(:, :, nz - iz + 2, 1) = swap(:, :)

                swap(:, :) = out_vfieldk(:, :, iz, 2)
                out_vfieldk(:, :, iz, 2) = out_vfieldk(:, :, nz - iz + 2, 2)
                out_vfieldk(:, :, nz - iz + 2, 2) = swap(:, :)

            end if
            
        end do
            
    end subroutine symmops_mirrorz

!==============================================================================

    subroutine symmops_shiftx(sx, in_vfieldk, out_vfieldk)
        
        complex(dpc), intent(inout) :: in_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: out_vfieldk(:, :, :, :)
        complex(dpc) :: temp_vfieldk(ny_half, nz, 3)

        real(dp), intent(in) :: sx
        integer(i4) :: ix
        
        out_vfieldk(:,:,:,1:3) = in_vfieldk(:,:,:,1:3)
        
        do ix = 1, nx_perproc
            if(ix_max /= -1 .and. ix == ix_max) cycle

            temp_vfieldk(:, :, 1:3) = in_vfieldk(ix, :, :, 1:3)
            
            ! Real part:
            out_vfieldk(ix, :, :, 1:3)%re = temp_vfieldk(:, :, 1:3)%re * cos(kx(ix) * sx) &
                                        + temp_vfieldk(:, :, 1:3)%im * sin(kx(ix) * sx)
            
            ! Imaginary part: 
            out_vfieldk(ix, :, :, 1:3)%im = temp_vfieldk(:, :, 1:3)%im * cos(kx(ix) * sx) &
                                        - temp_vfieldk(:, :, 1:3)%re * sin(kx(ix) * sx)
            
        end do
            
    end subroutine symmops_shiftx

!==============================================================================

    subroutine symmops_halfshiftx(in_vfieldk, out_vfieldk)
        
        complex(dpc), intent(inout) :: in_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: out_vfieldk(:, :, :, :)

        integer(i4) :: ix
        
        out_vfieldk(:,:,:,1:3) = in_vfieldk(:,:,:,1:3)
        
        do ix = 1, nx_perproc
            if(ix_max /= -1 .and. ix == ix_max) cycle

            if (mod(qx(ix), 2) /= 0) then
            
                out_vfieldk(ix, :, :, 1:3) = -out_vfieldk(ix, :, :, 1:3)
            
            end if
        end do
        
    end subroutine symmops_halfshiftx

!==============================================================================

    subroutine symmops_halfshifty(in_vfieldxx, out_vfieldxx)
        
        real(dp), intent(inout) :: in_vfieldxx(:, :, :, :)
        real(dp), intent(out) :: out_vfieldxx(:, :, :, :)        

        integer(i4) :: n
        real(dp) :: swap(nyy/2, nzz_perproc, nxx)

        if (mod(nyy, 2) /= 0) then
            write(out,*) "symmops_halfshifty: can't work with odd nyy."
            flush(out)
            error stop
        end if

        out_vfieldxx(:,:,:,1:3) = in_vfieldxx(:,:,:,1:3)
        
        do n=1,3
            swap(:, :, :) = out_vfieldxx(nyy/2 + 1 :nyy, :, :, n)
            out_vfieldxx(nyy/2 + 1 :nyy, :, :, n) = out_vfieldxx(1:nyy/2,:,:,n)
            out_vfieldxx(1:nyy/2,:,:,n) = swap(:, :, :)
        end do

    end subroutine symmops_halfshifty

!==============================================================================

    subroutine symmops_shiftz(sz, in_vfieldk, out_vfieldk)
        
        complex(dpc), intent(inout) :: in_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: out_vfieldk(:, :, :, :)
        complex(dpc) :: temp_vfieldk(nx_perproc, ny_half, 3)
        
        real(dp), intent(in) :: sz
        integer(i4) :: iz
        
        out_vfieldk(:, :, :, 1:3) = in_vfieldk(:, :, :, 1:3) 
        
       do iz = 1, nz
            if(iz == iz_max) cycle

            temp_vfieldk(:, :, 1:3) = in_vfieldk(:, :, iz, 1:3)
           
            ! Real part:
            out_vfieldk(:, :, iz, 1:3)%re = temp_vfieldk(:, :, 1:3)%re * cos(kz(iz) * sz) &
                                        + temp_vfieldk(:, :, 1:3)%im * sin(kz(iz) * sz)
            
            ! Imaginary part: 
            out_vfieldk(:, :, iz, 1:3)%im = temp_vfieldk(:, :, 1:3)%im * cos(kz(iz) * sz) &
                                        - temp_vfieldk(:, :, 1:3)%re * sin(kz(iz) * sz)
            
        end do
            
    end subroutine symmops_shiftz

!==============================================================================

    subroutine symmops_halfshiftz(in_vfieldk, out_vfieldk)

        complex(dpc), intent(inout) :: in_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: out_vfieldk(:, :, :, :)
        
        ! Shift fields in z-direction by Lz/2
        integer(i4) :: iz
        
        out_vfieldk(:, :, :, 1:3) = in_vfieldk(:, :, :, 1:3) 

        ! in Fourier space, j indexes z-coordinate
        do iz = 1, nz
            if(iz == iz_max) cycle;

            if (mod(qz(iz), 2) /= 0) then
            
                out_vfieldk(:, :, iz, 1:3) = -out_vfieldk(:, :, iz, 1:3)
            
            end if

        end do
            
    end subroutine symmops_halfshiftz

!==============================================================================

    subroutine symmops_Sx(in_vfieldxx, out_vfieldxx)

        real(dp), intent(inout) :: in_vfieldxx(:, :, :, :)
        real(dp), intent(out) :: out_vfieldxx(:, :, :, :)
        
        call symmops_halfshifty(in_vfieldxx, out_vfieldxx)
        call symmops_mirrorx(out_vfieldxx, out_vfieldxx)
        
    end subroutine symmops_Sx

!==============================================================================

    subroutine symmops_Sy(in_vfieldxx, out_vfieldxx)

        real(dp), intent(inout) :: in_vfieldxx(:, :, :, :)
        real(dp), intent(out) :: out_vfieldxx(:, :, :, :)
        
        call symmops_halfshifty(in_vfieldxx, out_vfieldxx)
        call symmops_mirrory(out_vfieldxx, out_vfieldxx)
        
    end subroutine symmops_Sy

!==============================================================================

    subroutine symmops_Rxy(in_vfieldxx, out_vfieldxx)

        real(dp), intent(inout) :: in_vfieldxx(:, :, :, :)
        real(dp), intent(out) :: out_vfieldxx(:, :, :, :)

        call symmops_mirrory(in_vfieldxx, out_vfieldxx)
        call symmops_mirrorx(out_vfieldxx, out_vfieldxx)
        
    end subroutine symmops_Rxy

!==============================================================================
    
    subroutine symmops_kproject(ksymmetry, vfieldk)

        abstract interface
            subroutine symmops_ksymmetry(in_vfieldk, out_vfieldk)
                import :: dpc
                complex(dpc), intent(inout) :: in_vfieldk(:, :, :, :)
                complex(dpc), intent(out)   :: out_vfieldk(:, :, :, :)
            end subroutine
        end interface

        procedure(symmops_ksymmetry) :: ksymmetry
        complex(dpc), intent(inout) :: vfieldk(:, :, :, :)
        complex(dpc) :: symmopps_vfieldk(nx_perproc, ny_half, nz, 3)

        call ksymmetry(vfieldk(:,:,:,1:3), symmopps_vfieldk(:,:,:,1:3))
        vfieldk(:,:,:,1:3) = (vfieldk(:,:,:,1:3) + symmopps_vfieldk(:,:,:,1:3)) / 2.0_dp
        
    end subroutine symmops_kproject

!==============================================================================

    subroutine symmops_xxproject(xxsymmetry, vfieldxx)

        abstract interface
            subroutine symmops_xxsymmetry(in_vfieldxx, out_vfieldxx)
                import :: dp
                real(dp), intent(inout) :: in_vfieldxx(:, :, :, :)
                real(dp), intent(out)   :: out_vfieldxx(:, :, :, :)
            end subroutine
        end interface

        procedure(symmops_xxsymmetry) :: xxsymmetry
        real(dp), intent(inout) :: vfieldxx(:, :, :, :)
        real(dp) :: symmopps_vfieldxx(nyy, nzz_perproc, nxx, 3)

        call xxsymmetry(vfieldxx(:,:,:,1:3), symmopps_vfieldxx(:,:,:,1:3))
        vfieldxx(:,:,:,1:3) = (vfieldxx(:,:,:,1:3) + symmopps_vfieldxx(:,:,:,1:3)) / 2.0_dp
        
    end subroutine symmops_xxproject

!==============================================================================

    subroutine symmopps_project_all(in_vfieldk, out_vfieldk)
        complex(dpc), intent(inout) :: in_vfieldk(:, :, :, :)
        complex(dpc), intent(out)   :: out_vfieldk(:, :, :, :)
        real(dp) :: out_vfieldxx(nyy, nzz_perproc, nxx, 3)

        out_vfieldk(:,:,:,1:3) = in_vfieldk(:,:,:,1:3)

        if (Rz) call symmops_kproject(symmops_mirrorz, out_vfieldk)

        if (Rxy .or. Ry .or. Sx .or. Sy) then

            call fftw_vk2x(out_vfieldk, out_vfieldxx)

            if (Rxy) call symmops_xxproject(symmops_Rxy, out_vfieldxx)

            if (Ry) call symmops_xxproject(symmops_mirrory, out_vfieldxx)

            if (Sx) call symmops_xxproject(symmops_Sx, out_vfieldxx)

            if (Sy) call symmops_xxproject(symmops_Sy, out_vfieldxx)

            call fftw_vx2k(out_vfieldxx, out_vfieldk)

        endif

    end subroutine symmopps_project_all

!==============================================================================

end module symmops

#include "macros.h"

module test_symmops
    use numbers
    use openmpi
    use io
    use parameters
    use fftw
    use symmops

    contains 

!==============================================================================

    subroutine test_symmops_mirrorx(in_vfieldk, out_vfieldk)

        complex(dpc), intent(inout) :: in_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: out_vfieldk(:, :, :, :)
        complex(dpc) :: swap(ny_half,nz)
        
        integer(i4) :: ix
        
        out_vfieldk(:,:,:,1:3) = in_vfieldk(:,:,:,1:3)

        do ix = 1, nx / 2 ! ix = nx / 2 + 1 is zeroed

            if (ix == 1) then
                out_vfieldk(ix, :, :, 1) = -out_vfieldk(ix, :, :, 1)
                out_vfieldk(ix, :, :, 2:3) = out_vfieldk(ix, :, :, 2:3)

            else
                swap(:, :) = out_vfieldk(ix, :, :, 1)
                out_vfieldk(ix, :, :, 1) = -out_vfieldk(nx - ix + 2, :, :, 1)
                out_vfieldk(nx - ix + 2, :, :, 1) = -swap(:,:)

                swap(:, :) = out_vfieldk(ix, :, :, 2)
                out_vfieldk(ix, :, :, 2) = out_vfieldk(nx - ix + 2, :, :, 2)
                out_vfieldk(nx - ix + 2, :, :, 2) = swap(:,:)

                swap(:, :) = out_vfieldk(ix, :, :, 3)
                out_vfieldk(ix, :, :, 3) = out_vfieldk(nx - ix + 2, :, :, 3)
                out_vfieldk(nx - ix + 2, :, :, 3) = swap(:,:)

            end if
            
        end do
            
    end subroutine test_symmops_mirrorx

!==============================================================================

    subroutine test_symmops_mirrory(in_vfieldk, out_vfieldk)

        complex(dpc), intent(in)  :: in_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: out_vfieldk(:, :, :, :)
        
        _indices

        out_vfieldk(:, :, :, 1:3) = in_vfieldk(:, :, :, 1:3) 
        
        do iy = 1, ny_half

            if (iy == 1) then
                ! w(0, kx, ky) -> -w(0, kx, ky)
                out_vfieldk(:, iy, :, 2) = -in_vfieldk(:, iy, :, 2)
                ! [u,v](0, kx, ky) -> [u,v](0, kx, ky)
                out_vfieldk(:, iy, :, (/1, 3/)) = in_vfieldk(:, iy, :, (/1, 3/))

            else

                do ix=1, nx; do iz=1,nz
                    if(iz == iz_max) cycle; if(ix_max /= -1 .and. ix == ix_max) cycle;
                    if (ix==1 .and. iz==1) then
                        ! Re w(kz, 0, 0) -> -Re w(kz, 0, 0)
                        out_vfieldk(ix,iy,iz, 2)%re = -in_vfieldk(ix,iy,iz, 2)%re
                        ! Re [u,v](kz, 0, 0) -> Re [u,v](kz, 0, 0)
                        out_vfieldk(ix,iy,iz, (/1, 3/))%re = in_vfieldk(ix,iy,iz, (/1, 3/))%re
                        ! Im w(kz, 0, 0) -> Im w(kz, 0, 0)
                        out_vfieldk(ix,iy,iz, 2)%im = in_vfieldk(ix,iy,iz, 2)%im
                        out_vfieldk(ix,iy,iz, (/1, 3/))%im = -in_vfieldk(ix,iy,iz, (/1, 3/))%im
                    else if (ix==1 .and. iz/=1) then
                        out_vfieldk(ix,iy,iz, 2)%re = -in_vfieldk(ix,iy,nz-iz+2, 2)%re
                        out_vfieldk(ix,iy,iz, (/1, 3/))%re = in_vfieldk(ix,iy,nz-iz+2, (/1, 3/))%re
                        out_vfieldk(ix,iy,iz, 2)%im = in_vfieldk(ix,iy,nz-iz+2, 2)%im
                        out_vfieldk(ix,iy,iz, (/1, 3/))%im = -in_vfieldk(ix,iy,nz-iz+2, (/1, 3/))%im
                    else if (ix/=1 .and. iz==1) then
                        out_vfieldk(ix,iy,iz, 2)%re = -in_vfieldk(nx-ix+2,iy,iz, 2)%re
                        out_vfieldk(ix,iy,iz, (/1, 3/))%re = in_vfieldk(nx-ix+2,iy,iz, (/1, 3/))%re
                        out_vfieldk(ix,iy,iz, 2)%im = in_vfieldk(nx-ix+2,iy,iz, 2)%im
                        out_vfieldk(ix,iy,iz, (/1, 3/))%im = -in_vfieldk(nx-ix+2,iy,iz, (/1, 3/))%im
                    else

                        ! Real parts
                        out_vfieldk(ix,iy,iz, 2)%re = -in_vfieldk(nx-ix+2,iy,nz-iz+2, 2)%re
                        out_vfieldk(ix,iy,iz, (/1, 3/))%re = in_vfieldk(nx-ix+2,iy,nz-iz+2, (/1, 3/))%re
                        ! Imaginary parts
                        out_vfieldk(ix,iy,iz, 2)%im = in_vfieldk(nx-ix+2,iy,nz-iz+2, 2)%im
                        out_vfieldk(ix,iy,iz, (/1, 3/))%im = -in_vfieldk(nx-ix+2,iy,nz-iz+2, (/1, 3/))%im

                    end if

                end do; end do

            end if
            
        end do
            
    end subroutine test_symmops_mirrory

!==============================================================================

    subroutine test_symmops_mirrorz(in_vfieldxx, out_vfieldxx)
        real(dp), intent(inout)  :: in_vfieldxx(:, :, :, :)
        real(dp), intent(out) :: out_vfieldxx(:, :, :, :)

        integer(i4) :: izz
        real(dp) :: swap(nyy, nxx)

        out_vfieldxx(:,:,:,1:3) = in_vfieldxx(:,:,:,1:3)

        ! Reflect in x
        do izz = 1, nzz / 2 + 1
            
            if (izz == 1 .or. (mod(nzz, 2) == 0 .and. izz == nzz / 2 + 1)) then
            
                ! sigma_x u(x) = - u(- x)
                out_vfieldxx(:, izz, :, 3) = -out_vfieldxx(:, izz, :, 3)
                        
            else             
            
                ! sigma_x u(x) = - u(- x)
                swap(:, :) = out_vfieldxx(:, izz, :, 3)
                out_vfieldxx(:, izz, :, 3) = - out_vfieldxx(:, nzz - (izz - 2), :, 3)
                out_vfieldxx(:, nzz - (izz - 2), :, 3) = -swap(:, :)

                swap(:, :) = out_vfieldxx(:, izz, :, 1)
                out_vfieldxx(:, izz, :, 1) = out_vfieldxx(:, nzz - (izz - 2), :, 1)
                out_vfieldxx(:, nzz - (izz - 2), :, 1) = swap(:, :)

                ! sigma_x [v,w](x) = [v,w](- x)
                swap(:, :) = out_vfieldxx(:, izz, :, 2)
                out_vfieldxx(:, izz, :, 2) = out_vfieldxx(:, nzz - (izz - 2), :, 2)
                out_vfieldxx(:, nzz - (izz - 2), :, 2) = swap(:, :)
            
            end if
            
        end do
        
    end subroutine test_symmops_mirrorz

!==============================================================================

    subroutine test_symmops_halfshiftx(in_vfieldxx, out_vfieldxx)
        
        real(dp), intent(inout) :: in_vfieldxx(:, :, :, :)
        real(dp), intent(out) :: out_vfieldxx(:, :, :, :)    

        integer(i4) :: n
        real(dp) :: swap(nyy, nzz_perproc, nxx/2)

        if (mod(nxx, 2) /= 0) then
            write(out,*) "test_symmops_halfshiftx: can't work with odd nxx."
            flush(out)
            error stop
        end if

        out_vfieldxx(:,:,:,1:3) = in_vfieldxx(:,:,:,1:3)
        
        ! Shift in y-direction by Lx/2:
        do n=1,3
            swap(:, :, :) = out_vfieldxx(:, :, nxx/2+1:nxx, n)
            out_vfieldxx(:, :, nxx/2+1:nxx, n) = out_vfieldxx(:, :, 1:nxx/2,n)
            out_vfieldxx(:,:,1:nxx/2,n) = swap(:, :, :)
        end do

    end subroutine test_symmops_halfshiftx

!==============================================================================

    subroutine test_symmops_halfshifty(in_vfieldk, out_vfieldk)
        
        complex(dpc), intent(inout) :: in_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: out_vfieldk(:, :, :, :)
        
        integer(i4) :: iy

        out_vfieldk(:,:,:,1:3) = in_vfieldk(:,:,:,1:3)
        
        do iy = 1, ny_half

            if (mod(qy(iy), 2) /= 0) then

                out_vfieldk(:, iy, :, 1:3) = -out_vfieldk(:, iy, :, 1:3)
            
            end if
            
        end do
            
    end subroutine test_symmops_halfshifty

!==============================================================================

    subroutine test_symmops_halfshiftz(in_vfieldxx, out_vfieldxx)

        real(dp), intent(inout) :: in_vfieldxx(:, :, :, :)
        real(dp), intent(out) :: out_vfieldxx(:, :, :, :)

        integer(i4) :: n
        real(dp) :: swap(nyy, nzz/2, nxx)

        if (mod(nzz, 2) /= 0) then
            write(out,*) "test_symmops_halfshiftz: can't work with odd nzz."
            flush(out)
            error stop
        end if

        out_vfieldxx(:,:,:,1:3) = in_vfieldxx(:,:,:,1:3)
        
        do n=1,3
            swap(:, :, :) = out_vfieldxx(:, nzz/2+1:nzz, :, n)
            out_vfieldxx(:, nzz/2+1:nzz, :, n) = out_vfieldxx(:, 1:nzz/2, :,n)
            out_vfieldxx(:,1:nzz/2,:,n) = swap(:, :, :)
        end do

    end subroutine test_symmops_halfshiftz

!==============================================================================

end module test_symmops
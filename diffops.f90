#include "macros.h"

module diffops
    use numbers
    use openmpi
    use io
    use parameters
    use fftw

    contains 

!==============================================================================

    complex(dpc) function nabla(ix, iy, iz, d)
        integer(i4), intent(in) :: ix, iy, iz, d

        select case (d)
            case(1)
                nabla = imag_1 * kx(ix)
            case(2)
                nabla = imag_1 * ky(iy)
            case(3)
                nabla = imag_1 * kz(iz)
        end select
    end function nabla

!==============================================================================

    complex(dpc) function laplacian(ix, iy, iz)
        integer(i4), intent(in) :: ix, iy, iz

        laplacian = kx(ix)**2 + ky(iy)**2 + kz(iz)**2
    end function laplacian

!==============================================================================

    complex(dpc) function inverse_laplacian(ix, iy, iz)
        integer(i4), intent(in) :: ix, iy, iz
        
        if (qx(ix) == 0 .and. qy(iy) == 0 .and. qz(iz) == 0) then
            inverse_laplacian = 0
        else
            inverse_laplacian = 1.0_dp / laplacian(ix,iy,iz)
        end if
    end function inverse_laplacian

!==============================================================================

    subroutine diffops_div(vfieldk, div_vfieldk)

        complex(dpc), intent(in)  :: vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: div_vfieldk(:, :, :)
        integer(i4) :: dim
        _indices
        
        div_vfieldk = 0

        do dim = 1, 3
            _loop_spec_begin
                div_vfieldk(ix, iy, iz) = div_vfieldk(ix, iy, iz) + &
                    nabla(ix, iy, iz, dim) * vfieldk(ix, iy, iz, dim)
            _loop_spec_end
        end do

    end subroutine diffops_div

!==============================================================================

    subroutine diffops_partx(sfieldk, partx_sfieldk)

        complex(dpc), intent(in)  :: sfieldk(:, :, :)
        complex(dpc), intent(out) :: partx_sfieldk(:, :, :)
        _indices
        
        _loop_spec_begin
            partx_sfieldk(ix, iy, iz) = imag_1 * kx(ix) * sfieldk(ix, iy, iz)
        _loop_spec_end
        
    end subroutine diffops_partx
    
!==============================================================================

    subroutine diffops_party(sfieldk, party_sfieldk)

        complex(dpc), intent(in)  :: sfieldk(:, :, :)
        complex(dpc), intent(out) :: party_sfieldk(:, :, :)
        _indices
        
        _loop_spec_begin
            party_sfieldk(ix, iy, iz) = imag_1 * ky(iy) * sfieldk(ix, iy, iz)
        _loop_spec_end
        
    end subroutine diffops_party
    
!==============================================================================

    subroutine diffops_partz(sfieldk, partz_sfieldk)

        complex(dpc), intent(in)  :: sfieldk(:, :, :)
        complex(dpc), intent(out) :: partz_sfieldk(:, :, :)
        _indices
        
        _loop_spec_begin
            partz_sfieldk(ix, iy, iz) = imag_1 * kz(iz) * sfieldk(ix, iy, iz)
        _loop_spec_end
        
    end subroutine diffops_partz
    
!==============================================================================  

end module diffops
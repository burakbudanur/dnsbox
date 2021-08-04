#include "macros.h"
module symred
    use numbers
    use openmpi
    use io
    use parameters
    use fieldio
    use fftw
    use symmops
    use vfield
    
    ! templates for reducing translation symmetries in x and
    ! u_xp(:,:,:,1:3) will keep the template
    ! u_txp(:,:,:,1:3) will keep the shifted version
    ! vice versa for z
    complex(dpc), allocatable, dimension(:, :, :, :) :: &
        u_xp, u_txp, u_zp, u_tzp

    ! phase angles
    real(dp) :: phi_x, phi_z
    
    ! output to write phases
    integer(i4) :: phases_ch
    character(255) :: phases_file = 'phases.gp'
    logical :: phases_written = .false.

    contains 

    subroutine symred_init

        allocate(u_xp(nx_perproc,ny_half,nz,3))
        allocate(u_txp, u_zp, u_tzp, mold=u_xp)

        ! read the templates
        fname = 'u_xp.000000'
        call fieldio_read_nocompact(u_xp(:,:,:,1:3))

        fname = 'u_zp.000000'
        call fieldio_read_nocompact(u_zp(:,:,:,1:3))

        ! prepare the shifted versions
        call symmops_shiftx(Lx/4, u_xp(:,:,:,1:3), u_txp(:,:,:,1:3))
        call symmops_shiftz(Lz/4, u_zp(:,:,:,1:3), u_tzp(:,:,:,1:3))

    end subroutine symred_init

!==============================================================================

    subroutine symred_phases_write
        logical :: there2
        
        if (my_id==0) then

            inquire(file=TRIM(phases_file), exist=there, opened=there2)
            if (.not.there) then
            open(newunit=phases_ch,file=TRIM(phases_file),form='formatted')
                write(phases_ch,"(A2,"//i4_len//","//"3"//sp_len//")") &
                    "# ", "itime", "time", "phi_x", "phi_z"
            end if
            if(there.and..not.there2) then
            open(newunit=phases_ch,file=TRIM(phases_file),position='append')
            end if
            write(phases_ch,"(A2,"//i4_f//","//"3"//sp_f//")")&
                "  ", itime, time, phi_x, phi_z

            phases_written = .true.

        end if

    end subroutine symred_phases_write

!==============================================================================

    subroutine symred_slice(in_vfieldk, out_vfieldk)
        complex(dpc), intent(in) :: in_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: out_vfieldk(:, :, :, :)
        real(dp) :: phi_x_r, phi_x_i, phi_z_r, phi_z_i

        ! slice in x

        call vfield_inprod(in_vfieldk(:,:,:,1:3), u_xp(:,:,:,1:3), phi_x_r, .true.)
        call vfield_inprod(in_vfieldk(:,:,:,1:3), u_txp(:,:,:,1:3), phi_x_i, .true.)

        phi_x = atan2(phi_x_i, phi_x_r)
        call symmops_shiftx(-(phi_x/(2*pi))*Lx, in_vfieldk(:,:,:,1:3), out_vfieldk(:,:,:,1:3))

        ! slice in z

        call vfield_inprod(out_vfieldk(:,:,:,1:3), u_zp(:,:,:,1:3), phi_z_r, .true.)
        call vfield_inprod(out_vfieldk(:,:,:,1:3), u_tzp(:,:,:,1:3), phi_z_i, .true.)

        phi_z = atan2(phi_z_i, phi_z_r)
        call symmops_shiftz(-(phi_z/(2*pi))*Lz, out_vfieldk(:,:,:,1:3), out_vfieldk(:,:,:,1:3))

    end subroutine symred_slice

end module symred

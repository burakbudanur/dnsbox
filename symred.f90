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
    integer(i4) :: phases_ch, slice_proj_ch_x, slice_proj_ch_z
    character(255) :: phases_file = 'phases.gp'
    logical :: phases_written = .false., slice_proj_written = .false.

    character(255) :: slice_proj_num_bases_str
    character(255) :: slice_proj_results_format

    type(C_PTR) :: p_projections_x, p_projections_z
    real(dp), pointer :: projections_x_rview(:), projections_z_rview(:)
    complex(dpc), pointer :: projections_x(:), projections_z(:)

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

    subroutine symred_proj_init

        write(slice_proj_num_bases_str, *) 6*ny_half - 4
        slice_proj_results_format = "(A2,"//i4_f//","//sp_f//","//TRIM(slice_proj_num_bases_str)//dp_f//")"

        ! projections_x and projections_x_rview share memory
        p_projections_x = fftw_alloc_complex(int(3*ny_half-2, i8))
        call c_f_pointer(&
            p_projections_x, projections_x, [3*ny_half-2])
        call c_f_pointer(&
            p_projections_x, projections_x_rview, [6*ny_half-4])

        p_projections_z = fftw_alloc_complex(int(3*ny_half-2, i8))
        call c_f_pointer(&
            p_projections_z, projections_z, [3*ny_half-2])
        call c_f_pointer(&
            p_projections_z, projections_z_rview, [6*ny_half-4])

        if (my_id == 0) then
            open(newunit=slice_proj_ch_x,file='slice_projections_x.gp',status='replace')
            write(slice_proj_ch_x, "(A2,"//i4_len//","//sp_len//","//dp_len//")") "# ", "itime", "time", "projections"
            
            open(newunit=slice_proj_ch_z,file='slice_projections_z.gp',status='replace')
            write(slice_proj_ch_z, "(A2,"//i4_len//","//sp_len//","//dp_len//")") "# ", "itime", "time", "projections"

            slice_proj_written = .true.
        end if
    end subroutine symred_proj_init

!==============================================================================

    subroutine symred_phases_write
        logical :: there2
        
        if (my_id==0) then

            inquire(file=TRIM(phases_file), exist=there, opened=there2)
            if (.not.there) then
            open(newunit=phases_ch,file=TRIM(phases_file),form='formatted')
                write(phases_ch,"(A2,"//i4_len//","//sp_len//","//"2"//dp_len//")") &
                    "# ", "itime", "time", "phi_x", "phi_z"
            end if
            if(there.and..not.there2) then
            open(newunit=phases_ch,file=TRIM(phases_file),position='append')
            end if
            write(phases_ch,"(A2,"//i4_f//","//sp_f//","//"2"//dp_f//")")&
                "  ", itime, time, phi_x, phi_z

            phases_written = .true.

        end if

    end subroutine symred_phases_write

!==============================================================================

    subroutine symred_slice(in_vfieldk, out_vfieldk)
        complex(dpc), intent(inout) :: in_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: out_vfieldk(:, :, :, :)
        real(dp) :: phi_x_r, phi_x_i, phi_z_r, phi_z_i

        ! slice in x

        call vfield_inprod(in_vfieldk(:,:,:,1:3), u_xp(:,:,:,1:3), phi_x_r, .true.)
        call vfield_inprod(in_vfieldk(:,:,:,1:3), u_txp(:,:,:,1:3), phi_x_i, .true.)

        phi_x = atan2(phi_x_i, phi_x_r)
        call symmops_shiftx(-(phi_x/(2.0_dp*pi))*Lx, in_vfieldk(:,:,:,1:3), out_vfieldk(:,:,:,1:3))

        ! slice in z

        call vfield_inprod(out_vfieldk(:,:,:,1:3), u_zp(:,:,:,1:3), phi_z_r, .true.)
        call vfield_inprod(out_vfieldk(:,:,:,1:3), u_tzp(:,:,:,1:3), phi_z_i, .true.)

        phi_z = atan2(phi_z_i, phi_z_r)
        call symmops_shiftz(-(phi_z/(2.0_dp*pi))*Lz, out_vfieldk(:,:,:,1:3), out_vfieldk(:,:,:,1:3))

    end subroutine symred_slice

!==============================================================================

    subroutine symred_projections(in_vfieldk)
        ! valid only if Ry
        complex(dpc), intent(in) :: in_vfieldk(:, :, :, :)
        complex(dpc) :: my_projections_x(3*ny_half-2), &
                        my_projections_z(3*ny_half-2)

        my_projections_x(:) = 0
        my_projections_z(:) = 0

        if (ix_first_p /= -1) then
            my_projections_x(1:ny_half-1) = &
                my_projections_x(1:ny_half-1) + in_vfieldk(ix_first_p,2:,1,1) / 2
            my_projections_x(ny_half:2*ny_half-2) = &
                my_projections_x(ny_half:2*ny_half-2) - imag_1 * in_vfieldk(ix_first_p,2:,1,2) / 2
            my_projections_x(2*ny_half-1) = &
                my_projections_x(2*ny_half-1) + in_vfieldk(ix_first_p,1,1,3) / 2 / sqrt(2.0_dp)
            my_projections_x(2*ny_half:) = &
                my_projections_x(2*ny_half:) + in_vfieldk(ix_first_p,2:,1,3) / 2
        end if

        if (ix_first_n /= -1) then
            my_projections_x(1:ny_half-1) = &
                my_projections_x(1:ny_half-1) + conjg(in_vfieldk(ix_first_n,2:,1,1)) / 2
            my_projections_x(ny_half:2*ny_half-2) = &
                my_projections_x(ny_half:2*ny_half-2) + imag_1 * conjg(in_vfieldk(ix_first_n,2:,1,2)) / 2
            my_projections_x(2*ny_half-1) = &
                my_projections_x(2*ny_half-1) + conjg(in_vfieldk(ix_first_n,1,1,3)) / 2 / sqrt(2.0_dp)
            my_projections_x(2*ny_half:) = &
                my_projections_x(2*ny_half:) + conjg(in_vfieldk(ix_first_n,2:,1,3)) / 2
        end if

        if (ix_zero /= -1) then
            my_projections_z(1) = &
                (in_vfieldk(ix_zero,1,2,1) + conjg(in_vfieldk(ix_zero,1,nz,1))) / 2 / sqrt(2.0_dp)
            my_projections_z(2:ny_half) = &
                (in_vfieldk(ix_zero,2:,2,1) + conjg(in_vfieldk(ix_zero,2:,nz,1))) / 2
            my_projections_z(ny_half+1:2*ny_half-1) = &
                -imag_1 * (in_vfieldk(ix_zero,2:,2,2) - conjg(in_vfieldk(ix_zero,2:,nz,2))) / 2
            my_projections_z(2*ny_half:) = &
                (in_vfieldk(ix_zero,2:,2,3) + conjg(in_vfieldk(ix_zero,2:,nz,3))) / 2
        end if

        call MPI_REDUCE(my_projections_x, projections_x, 3*ny_half-2, MPI_COMPLEX16, MPI_SUM, 0, &
        MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(my_projections_z, projections_z, 3*ny_half-2, MPI_COMPLEX16, MPI_SUM, 0, &
        MPI_COMM_WORLD, mpi_err)

        if (my_id == 0) then

            write(slice_proj_ch_x, TRIM(slice_proj_results_format)) "  ", itime, time, projections_x_rview
            write(slice_proj_ch_z, TRIM(slice_proj_results_format)) "  ", itime, time, projections_z_rview

        end if

    end subroutine symred_projections

end module symred

module projector

    use numbers
    use openmpi
    use io
    use parameters
    use fieldio
    use fftw
    use vfield

    complex(dpc), allocatable :: proj_bases(:,:,:,:,:), proj_origin(:, :, :, :)
    real(dp), allocatable :: proj_results(:)

    real(dp)    :: proj_ekin
    integer(i4) :: proj_ch, proj_uindex
    
    character(1), parameter :: proj_uname = 'u'
    character(255) :: proj_num_bases_str
    character(255) :: proj_results_format
    character(255) :: proj_dir

    logical :: proj_written = .false.

    contains

!==============================================================================

    subroutine projector_init
    
        integer :: i, un

        write(proj_num_bases_str, *) num_proj_bases
        proj_results_format = "(A2,"//i4_f//","//sp_f//","//dp_f//","//TRIM(proj_num_bases_str)//dp_f//")"

        if (my_id == 0) allocate(proj_results(num_proj_bases))
        if (subtract_origin) allocate(proj_origin(nx_perproc, ny_half, nz, 3))
        allocate(proj_bases(nx_perproc, ny_half, nz, 3,num_proj_bases))

        inquire(file='basis.in', exist=there)
        if (.not. there) then
            write(out,*) 'basis.in not found.'
            flush(out)
            stop
        else
            open(newunit=un,status='old',file='basis.in')
            read(un,'(A)') proj_dir
            close(un)

            write(out, *) 'proj_dir:'
            write(out, *) proj_dir
        end if

        ! Read the average
        if (subtract_origin) then
            proj_uindex = 0
            call projector_read(proj_origin)
        end if

        ! Read the modes
        do i=1, num_proj_bases
            proj_uindex = i
            call projector_read(proj_bases(:,:,:,:,i))
        end do
        
        ! Create the projections file, replace it if it exists
        if (my_id == 0) then
            open(newunit=proj_ch,file='projections.gp',status='replace')
            write(proj_ch, "(A2,"//i4_len//","//sp_len//","//dp_len//","//dp_len//")") "# ", "itime", "time", "proj_ekin", "projections"

            proj_written = .true.
        end if
    
    end subroutine projector_init

!==============================================================================
    
    subroutine projector_project(vfieldk)
        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        integer :: i
        real(dp) :: inprod
        complex(dpc) :: proj_vfieldk(nx_perproc, ny_half, nz, 3)

        ! Subtract the average
        if (subtract_origin) then
            proj_vfieldk(:,:,:,1:3) = vfieldk(:,:,:,1:3) - proj_origin(:,:,:,1:3)
        else
            proj_vfieldk(:,:,:,1:3) = vfieldk(:,:,:,1:3)
        end if

        do i = 1, num_proj_bases
            call vfield_inprod(proj_vfieldk, proj_bases(:,:,:,1:3,i), inprod, .false.)
            if (my_id == 0) proj_results(i) = inprod
        end do

        ! Compute the energy
        call vfield_norm2(proj_vfieldk, proj_ekin, .false.)
    
    end subroutine projector_project

!==============================================================================
    
    subroutine projector_write

        ! Write the projections
        if (my_id == 0) then
            write(proj_ch, TRIM(proj_results_format)) "  ", itime, time, proj_ekin, proj_results
        end if
    
    end subroutine projector_write

!==============================================================================

    subroutine projector_read(vfieldk)
        complex(dpc), intent(out) :: vfieldk(:, :, :, :)

        write(file_ext, "(i6.6)") proj_uindex
        fname = TRIM(proj_dir)//'/'//proj_uname//'.'//file_ext
        call fieldio_read(vfieldk)

    end subroutine projector_read

!==============================================================================

end module projector

! Module that contains working arrays

module m_work
    use m_numbers
    use m_parameters
    use m_io
    
    real(dp), allocatable :: wrk(:, :, :, :)

    contains 

!==============================================================================
    
    subroutine m_work_init
        use m_numbers
        use m_parameters
        
        
        integer :: ierr
        
        ierr = 0
        
        ! main working array, needed for FFT etc. (nz + 2, ny, nx)
        allocate(wrk(nz+2, ny, nx, 6), stat=ierr)
        
        if (ierr/=0) then
            write(out, *) "Cannot allocate work arrays, stopping"
            flush(out)
            stop
        end if
        
        write(out, "('Allocated work arrays')")

        flush(out)
                
    end subroutine m_work_init

!==============================================================================
    
    subroutine m_work_write

        !  This routine writes out wrk(:, :, :, 1:3)
        !  The file is written using the collective write (MPI-2 standard)        

        

        integer :: n, ierr

        TYPE(MPI_File) :: fh
        integer(kind=MPI_OFFSET_KIND)  :: offset

        integer :: nz1, ny1, nx1
        real(dp), allocatable :: wrk_buf(:,:,:)

        !--------------- writing process ------------------

        allocate(wrk_buf(nz + 2,ny,nx),stat=ierr)
        if (ierr/=0) then
            write(out,*) '*** m_work_write: cannot allocate wrk_buf'
            flush(out)
            stop
        end if
        wrk_buf = zero

        ! opening the file
        call MPI_INFO_CREATE(mpi_info_var, mpi_err)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(fname), &
                           MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                           mpi_info_var, fh, mpi_err)
        
        ! the master node writes the header with parameters
        if (myid==0) then
            nz1 = nz;  ny1 = ny;  nx1 = nx_all; 
            count = 1
            call MPI_FILE_WRITE(fh,  nz1, count, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  ny1, count, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  nx1, count, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh, TIME, count,    MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,   DT, count,    MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,   gamma, count, MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,   nu, count,    MPI_REAL8, mpi_status_var, mpi_err)
        end if

        ! all nodes write their stuff into the file
        
        writing_fields: do n = 1, 3

            offset = 44 + (n-1)*(nz+2)*ny*nx_all*8 + myid*(nz+2)*ny*nx*8
            count = (nz+2) * ny * nx

            wrk_buf(:, :, :) = wrk(:,:,:,n)
            call MPI_FILE_WRITE_AT(fh, offset, wrk_buf, count, MPI_REAL8, &
                                   mpi_status_var, mpi_err)

        end do writing_fields

        call MPI_FILE_CLOSE(fh, mpi_err)
        call MPI_INFO_FREE(mpi_info_var, mpi_err)
        deallocate(wrk_buf)

        ! write(out,*) '------------------------------------------------'
        ! write(out,*) 'State file written (par): '//TRIM(fname)
        ! write(out,"(' State file time = ',f15.10,i7)") time, itime/IPRINT2
        ! write(out,*) '------------------------------------------------'
        
        return
        
    end subroutine m_work_write

!==============================================================================

    subroutine m_work_write_onedim

        !  This routine writes out wrk(1:nz, :, :, 1) and rounds down to real*4
        !  The file is written using the collective write (MPI-2 standard)

        ! Intended only to write real space data / g 200805

        implicit none

        integer :: n, ierr

        TYPE(MPI_File) :: fh
        integer(kind=MPI_OFFSET_KIND)  :: offset

        integer :: nz1, ny1, nx1
        real(sp), allocatable :: wrk_buf(:,:,:)

        !--------------- writing process ------------------

        allocate(wrk_buf(nz,ny,nx),stat=ierr)
        if (ierr/=0) then
            write(out,*) '*** m_work_write_onedim: cannot allocate wrk_buf'
            flush(out)
            stop
        end if
        wrk_buf = zero

        ! opening the file
        call MPI_INFO_CREATE(mpi_info_var, mpi_err)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(fname), &
                           MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                           mpi_info_var, fh, mpi_err)
        
        ! the master node writes the header with parameters
        if (myid==0) then
            nz1 = nz;  ny1 = ny;  nx1 = nx_all; 
            count = 1
            call MPI_FILE_WRITE(fh,  nz1, count, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  ny1, count, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  nx1, count, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh, TIME, count,    MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,   DT, count,    MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,   Lz, count, MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,   nu, count,    MPI_REAL8, mpi_status_var, mpi_err)
        end if

        ! all nodes write their stuff into the file

        offset = 44 + myid*nz*ny*nx*4
        count = nz * ny * nx

        wrk_buf(:,:,:) = real(wrk(1:nz,:,:,1), 4)
        call MPI_FILE_WRITE_AT(fh, offset, wrk_buf, count, MPI_REAL4, &
                                mpi_status_var, mpi_err)
      

        call MPI_FILE_CLOSE(fh, mpi_err)
        call MPI_INFO_FREE(mpi_info_var, mpi_err)
        deallocate(wrk_buf)

        write(out,*) '------------------------------------------------'
        write(out,*) 'visualize file written: ' // TRIM(fname)
        write(out,*) '------------------------------------------------'
        
        return
        
    end subroutine m_work_write_onedim

!==============================================================================

    subroutine m_work_read

        !  This routine reads in wrk(:, :, :, 1:3)
        !  The file is written using the collective write (MPI-2 standard)
    
        

        integer    :: nz1,ny1,nx1
        real(dp)       :: TIME1, DT1, gamma1, nu1
        
        integer :: n, ierr

        integer :: un ! dummy i/o channel

        TYPE(MPI_File) :: fh
        integer(kind=MPI_OFFSET_KIND)  :: offset
        real(dp), allocatable :: wrk_buf(:,:,:)
        
        ! checking if the restart file exists

        inquire(file=TRIM(fname),exist=there)
        if(.not.there) then
            write(out,*) '*** error: Cannot find file : '//TRIM(fname)
            flush(out)
            stop
        end if

        ! write(out,*) 'Reading from the file (par): ', TRIM(fname)
        allocate(wrk_buf(nz + 2,ny,nx),stat=ierr)
        if (ierr/=0) then
            write(out,*) '*** m_work_read: cannot allocate wrk_buf'
            flush(out)
            stop
        end if
        wrk_buf = zero
        

        ! ----------------------------------------------------------------------
        ! first reading the parameters from the restart file.
        ! the root process opens it and reads the parameters, then broadcasts 
        ! the parameters.  After that it's decided if the parameters make sense,
        ! how many scalars to read etc.
        ! ----------------------------------------------------------------------

        if (myid==0) then
            open(newunit=un,file=TRIM(fname),form='unformatted',access='stream')
            read(un) nz1, ny1, nx1, TIME1, DT1, gamma1, nu1
            close(un)
        end if

        call MPI_BCAST(nz1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(ny1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(nx1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        
        call MPI_BCAST(TIME1, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(  DT1, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(  gamma1, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(  nu1, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)

        ! checking if the array sizes are the same in .in file and restart file
        if (nz/=nz1 .or. ny/=ny1 .or. nx_all/=nx1) then
            write(out,*) '*** error: Dimensions are different'
            write(out,*) '***     .in file: ',nz,ny,nx_all
            write(out,*) '*** restart file: ',nz1,ny1,nx1
            flush(out)
            stop
        end if

    !----------------------------------------------------------------------
    !  ------------ reading the stuff from the restart file ---------------
    !----------------------------------------------------------------------

        ! opening the file
        call MPI_INFO_CREATE(mpi_info_var, mpi_err)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(fname), MPI_MODE_RDONLY, mpi_info_var, &
                           fh, mpi_err)

        reading_fields: do n = 1, 3

            ! write(out,"('Reading variable # ',i3)") n

            offset = 44 + (n-1)*(nz+2)*ny*nx_all*8 + myid*(nz+2)*ny*nx*8
            count = (nz+2) * ny * nx
            call MPI_FILE_READ_AT_ALL(fh, offset, wrk_buf, count, MPI_REAL8, mpi_status_var, mpi_err)
            wrk(:,:,:,n) = wrk_buf(:,:,:)
        end do reading_fields

        call MPI_FILE_CLOSE(fh, mpi_err)
        call MPI_INFO_FREE(mpi_info_var, mpi_err)

        ! write(out,*) "State file successfully read."
        
        return        
        
    end subroutine m_work_read
    
end module m_work

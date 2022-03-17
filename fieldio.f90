#include "macros.h"

module fieldio
    use numbers
    use openmpi
    use io
    use parameters
    use fftw

    contains 

!==============================================================================

    subroutine fieldio_write_xcompact(vfieldk)
   
        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        _indices
        integer(i4) :: n, count
        TYPE(MPI_File) :: fh
        integer(MPI_OFFSET_KIND)  :: offset
        complex(dpc), allocatable, dimension(:, :, :) :: &
            buf_sfieldk, buf2_sfieldk

        if (ix_max /= -1) then
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc - 1))
        else
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc))
        end if
        if (ix_zero /= -1) allocate(buf2_sfieldk(ny_half, nz - 1, 1))

        ! opening the file
        call MPI_INFO_CREATE(mpi_info_var, mpi_err)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(fname), &
                           MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                           mpi_info_var, fh, mpi_err)

        ! the master node writes the header with parameters
        if (my_id==0) then
            call MPI_FILE_WRITE(fh,  forcing, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  nx, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  ny, 1,MPI_INTEGER4,mpi_status_var,mpi_err)
            call MPI_FILE_WRITE(fh,  nz, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  Lx, 1,    MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  Lz, 1,    MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  Re, 1, MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  tilt_angle, 1, MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  dt, 1, MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  itime, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  time, 1,    MPI_REAL8, mpi_status_var, mpi_err)
        end if

        ! Make x the slowest index and
        ! skip over the zeroed modes in x and z in all cases

        ! write only kx=0 for u
        if (ix_zero /= -1) then
            do iz = 1, nz; if(iz == iz_max) cycle; do iy = 1, ny_half;

                if (iz < iz_max) then
                    buf2_sfieldk(iy,iz,1) = vfieldk(ix_zero,iy,iz,1)
                else
                    buf2_sfieldk(iy,iz-1,1) = vfieldk(ix_zero,iy,iz,1)
                end if

            end do; end do;

            offset = 68
            count = ny_half * (nz-1)
            call MPI_FILE_WRITE_AT(fh, offset, buf2_sfieldk, count, MPI_COMPLEX16, &
                                mpi_status_var, mpi_err)
            
        end if

        ! write all of v and w
        do n = 2, 3
            _loop_spec_begin
                if (ix_max /= -1) then
                    if (ix < ix_max) then
                        if (iz < iz_max) then
                            buf_sfieldk(iy,iz,ix) = vfieldk(ix,iy,iz,n)
                        else
                            buf_sfieldk(iy,iz-1,ix) = vfieldk(ix,iy,iz,n)
                        end if
                    else
                        if (iz < iz_max) then
                            buf_sfieldk(iy,iz,ix-1) = vfieldk(ix,iy,iz,n)
                        else
                            buf_sfieldk(iy,iz-1,ix-1) = vfieldk(ix,iy,iz,n)
                        end if
                    end if
                else 
                    if (iz < iz_max) then
                        buf_sfieldk(iy,iz,ix) = vfieldk(ix,iy,iz,n)
                    else
                        buf_sfieldk(iy,iz-1,ix) = vfieldk(ix,iy,iz,n)
                    end if
                end if
            _loop_spec_end
            
            offset = 68 + ny_half*(nz-1)*16 + (n-2)*(nx-1)*ny_half*(nz-1)*16 &
                                                + my_id*nx_perproc*ny_half*(nz-1)*16
            if (my_id > my_id_ix_max) offset = offset - ny_half*(nz-1)*16
            if (ix_max /= -1) then
                count = (nx_perproc-1) * ny_half * (nz-1)
            else
                count = nx_perproc * ny_half * (nz-1)
            end if
            call MPI_FILE_WRITE_AT_ALL(fh, offset, buf_sfieldk, count, MPI_COMPLEX16, &
                                mpi_status_var, mpi_err)

        end do

        deallocate(buf_sfieldk)
        if (ix_zero /= -1) deallocate(buf2_sfieldk)

        call MPI_FILE_CLOSE(fh, mpi_err)
        call MPI_INFO_FREE(mpi_info_var, mpi_err)
        
    end subroutine fieldio_write_xcompact

!==============================================================================
    
    subroutine fieldio_write_ycompact(vfieldk)

        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        _indices
        TYPE(MPI_File) :: fh
        integer(MPI_OFFSET_KIND)  :: offset
        integer(i4) :: count
        complex(dpc), allocatable, dimension(:, :, :) :: &
            buf_sfieldk, buf2_sfieldk

        if (ix_max /= -1) then
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc - 1))
            allocate(buf2_sfieldk(1, nz - 1, nx_perproc - 1))
        else
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc))
            allocate(buf2_sfieldk(1, nz - 1, nx_perproc))
        end if

        ! opening the file
        call MPI_INFO_CREATE(mpi_info_var, mpi_err)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(fname), &
                           MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                           mpi_info_var, fh, mpi_err)

        ! the master node writes the header with parameters
        if (my_id==0) then
            call MPI_FILE_WRITE(fh,  forcing, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  nx, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  ny, 1,MPI_INTEGER4,mpi_status_var,mpi_err)
            call MPI_FILE_WRITE(fh,  nz, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  Lx, 1,    MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  Lz, 1,    MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  Re, 1, MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  tilt_angle, 1, MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  dt, 1, MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  itime, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  time, 1,    MPI_REAL8, mpi_status_var, mpi_err)
        end if

        ! Make x the slowest index and
        ! skip over the zeroed modes in x and z in all cases

        ! write u completely
        _loop_spec_begin
            if (ix_max /= -1) then
                if (ix < ix_max) then
                    if (iz < iz_max) then
                        buf_sfieldk(iy,iz,ix) = vfieldk(ix,iy,iz,1)
                    else
                        buf_sfieldk(iy,iz-1,ix) = vfieldk(ix,iy,iz,1)
                    end if
                else
                    if (iz < iz_max) then
                        buf_sfieldk(iy,iz,ix-1) = vfieldk(ix,iy,iz,1)
                    else
                        buf_sfieldk(iy,iz-1,ix-1) = vfieldk(ix,iy,iz,1)
                    end if
                end if
            else 
                if (iz < iz_max) then
                    buf_sfieldk(iy,iz,ix) = vfieldk(ix,iy,iz,1)
                else
                    buf_sfieldk(iy,iz-1,ix) = vfieldk(ix,iy,iz,1)
                end if
            end if
        _loop_spec_end
        
        offset = 68 + my_id*nx_perproc*ny_half*(nz-1)*16
        if (my_id > my_id_ix_max) offset = offset - ny_half*(nz-1)*16
        if (ix_max /= -1) then
            count = (nx_perproc-1) * ny_half * (nz-1)
        else
            count = nx_perproc * ny_half * (nz-1)
        end if
        call MPI_FILE_WRITE_AT_ALL(fh, offset, buf_sfieldk, count, MPI_COMPLEX16, &
                            mpi_status_var, mpi_err)

        ! write only ky=0 for v
        do iz = 1, nz; if(iz == iz_max) cycle; do ix = 1, nx_perproc; &
                                        if(ix_max /= -1 .and. ix == ix_max) cycle;

            if (ix_max /= -1) then
                if (ix < ix_max) then
                    if (iz < iz_max) then
                        buf2_sfieldk(1,iz,ix) = vfieldk(ix,1,iz,2)
                    else
                        buf2_sfieldk(1,iz-1,ix) = vfieldk(ix,1,iz,2)
                    end if
                else
                    if (iz < iz_max) then
                        buf2_sfieldk(1,iz,ix-1) = vfieldk(ix,1,iz,2)
                    else
                        buf2_sfieldk(1,iz-1,ix-1) = vfieldk(ix,1,iz,2)
                    end if
                end if
            else 
                if (iz < iz_max) then
                    buf2_sfieldk(1,iz,ix) = vfieldk(ix,1,iz,2)
                else
                    buf2_sfieldk(1,iz-1,ix) = vfieldk(ix,1,iz,2)
                end if
            end if
        end do; end do;
        
        offset = 68 + (nx-1)*ny_half*(nz-1)*16 + my_id*nx_perproc*(nz-1)*16
        if (my_id > my_id_ix_max) offset = offset - (nz-1)*16
        if (ix_max /= -1) then
            count = (nx_perproc-1) * (nz-1)
        else
            count = nx_perproc * (nz-1)
        end if
        call MPI_FILE_WRITE_AT_ALL(fh, offset, buf2_sfieldk, count, MPI_COMPLEX16, &
                            mpi_status_var, mpi_err)

        ! write w completely
        _loop_spec_begin
        if (ix_max /= -1) then
            if (ix < ix_max) then
                if (iz < iz_max) then
                    buf_sfieldk(iy,iz,ix) = vfieldk(ix,iy,iz,3)
                else
                    buf_sfieldk(iy,iz-1,ix) = vfieldk(ix,iy,iz,3)
                end if
            else
                if (iz < iz_max) then
                    buf_sfieldk(iy,iz,ix-1) = vfieldk(ix,iy,iz,3)
                else
                    buf_sfieldk(iy,iz-1,ix-1) = vfieldk(ix,iy,iz,3)
                end if
            end if
        else 
            if (iz < iz_max) then
                buf_sfieldk(iy,iz,ix) = vfieldk(ix,iy,iz,3)
            else
                buf_sfieldk(iy,iz-1,ix) = vfieldk(ix,iy,iz,3)
            end if
        end if
        _loop_spec_end
    
        offset = 68 + (nx-1)*(ny_half+1)*(nz-1)*16 &
                                                + my_id*nx_perproc*ny_half*(nz-1)*16
        if (my_id > my_id_ix_max) offset = offset - ny_half*(nz-1)*16
        if (ix_max /= -1) then
            count = (nx_perproc-1) * ny_half * (nz-1)
        else
            count = nx_perproc * ny_half * (nz-1)
        end if
        call MPI_FILE_WRITE_AT_ALL(fh, offset, buf_sfieldk, count, MPI_COMPLEX16, &
                            mpi_status_var, mpi_err)


        deallocate(buf_sfieldk)
        deallocate(buf2_sfieldk)

        call MPI_FILE_CLOSE(fh, mpi_err)
        call MPI_INFO_FREE(mpi_info_var, mpi_err)
        
    end subroutine fieldio_write_ycompact

!==============================================================================

    subroutine fieldio_write_zcompact(vfieldk)

        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        _indices
        integer(i4) :: n, count
        TYPE(MPI_File) :: fh
        integer(MPI_OFFSET_KIND)  :: offset
        complex(dpc), allocatable, dimension(:, :, :) :: &
            buf_sfieldk, buf2_sfieldk

        if (ix_max /= -1) then
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc - 1))
            allocate(buf2_sfieldk(ny_half, 1, nx_perproc - 1))
        else
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc))
            allocate(buf2_sfieldk(ny_half, 1, nx_perproc))
        end if

        ! opening the file
        call MPI_INFO_CREATE(mpi_info_var, mpi_err)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(fname), &
                           MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                           mpi_info_var, fh, mpi_err)

        ! the master node writes the header with parameters
        if (my_id==0) then
            call MPI_FILE_WRITE(fh,  forcing, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  nx, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  ny, 1,MPI_INTEGER4,mpi_status_var,mpi_err)
            call MPI_FILE_WRITE(fh,  nz, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  Lx, 1,    MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  Lz, 1,    MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  Re, 1, MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  tilt_angle, 1, MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  dt, 1, MPI_REAL8, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  itime, 1, MPI_INTEGER4, mpi_status_var, mpi_err)
            call MPI_FILE_WRITE(fh,  time, 1,    MPI_REAL8, mpi_status_var, mpi_err)
        end if

        ! Make x the slowest index and
        ! skip over the zeroed modes in x and z in all cases

        ! write u and v completely
        do n = 1, 2
            _loop_spec_begin
                if (ix_max /= -1) then
                    if (ix < ix_max) then
                        if (iz < iz_max) then
                            buf_sfieldk(iy,iz,ix) = vfieldk(ix,iy,iz,n)
                        else
                            buf_sfieldk(iy,iz-1,ix) = vfieldk(ix,iy,iz,n)
                        end if
                    else
                        if (iz < iz_max) then
                            buf_sfieldk(iy,iz,ix-1) = vfieldk(ix,iy,iz,n)
                        else
                            buf_sfieldk(iy,iz-1,ix-1) = vfieldk(ix,iy,iz,n)
                        end if
                    end if
                else 
                    if (iz < iz_max) then
                        buf_sfieldk(iy,iz,ix) = vfieldk(ix,iy,iz,n)
                    else
                        buf_sfieldk(iy,iz-1,ix) = vfieldk(ix,iy,iz,n)
                    end if
                end if
            _loop_spec_end
            
            offset = 68 + (n-1)*(nx-1)*ny_half*(nz-1)*16 &
                                                + my_id*nx_perproc*ny_half*(nz-1)*16
            if (my_id > my_id_ix_max) offset = offset - ny_half*(nz-1)*16
            if (ix_max /= -1) then
                count = (nx_perproc-1) * ny_half * (nz-1)
            else
                count = nx_perproc * ny_half * (nz-1)
            end if
            call MPI_FILE_WRITE_AT_ALL(fh, offset, buf_sfieldk, count, MPI_COMPLEX16, &
                                mpi_status_var, mpi_err)

        end do

        ! write only kz=0 for w
        do iy = 1, ny_half; do ix = 1, nx_perproc; if(ix_max /= -1 .and. ix == ix_max) cycle;
            if (ix_max /= -1) then
                if (ix < ix_max) then
                    buf2_sfieldk(iy,1,ix) = vfieldk(ix,iy,1,3)
                else
                    buf2_sfieldk(iy,1,ix-1) = vfieldk(ix,iy,1,3)
                end if
            else 
                buf2_sfieldk(iy,1,ix) = vfieldk(ix,iy,1,3)
            end if
        end do; end do;
        
        offset = 68 + 2*(nx-1)*ny_half*(nz-1)*16 + my_id*nx_perproc*ny_half*16
        if (my_id > my_id_ix_max) offset = offset - ny_half*16
        if (ix_max /= -1) then
            count = (nx_perproc-1) * ny_half
        else
            count = nx_perproc * ny_half
        end if
        call MPI_FILE_WRITE_AT_ALL(fh, offset, buf2_sfieldk, count, MPI_COMPLEX16, &
                            mpi_status_var, mpi_err)


        deallocate(buf_sfieldk)
        deallocate(buf2_sfieldk)

        call MPI_FILE_CLOSE(fh, mpi_err)
        call MPI_INFO_FREE(mpi_info_var, mpi_err)
        
    end subroutine fieldio_write_zcompact

!==============================================================================

    subroutine fieldio_write(vfieldk)

        complex(dpc), intent(in) :: vfieldk(:, :, :, :)

        if (nx - 1 >= ny_half .and. nx - 1 >= nz - 1) then
            call fieldio_write_xcompact(vfieldk)
        elseif (ny_half >= nx - 1 .and. ny_half >= nz - 1) then
            call fieldio_write_ycompact(vfieldk)
        else
            call fieldio_write_zcompact(vfieldk)
        end if
    end subroutine fieldio_write

!==============================================================================

    subroutine fieldio_read_xcompact(vfieldk)

        complex(dpc), intent(out) :: vfieldk(:, :, :, :)
        _indices
        integer(i4) :: forcing1, nx1, ny1, nz1, n, un, count
        TYPE(MPI_File) :: fh
        integer(MPI_OFFSET_KIND)  :: offset
        complex(dpc), allocatable, dimension(:, :, :) :: &
            buf_sfieldk, buf2_sfieldk

        inquire(file=TRIM(fname),exist=there)
        if(.not.there) then
            write(out,*) 'fieldio_read: Stopping, cannot find file : '//TRIM(fname)
            flush(out)
            error stop

        end if

        if (ix_max /= -1) then
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc - 1))
        else
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc))
        end if
        if (ix_zero /= -1) allocate(buf2_sfieldk(ny_half, nz - 1, 1))

        if (my_id==0) then
            open(newunit=un,file=TRIM(fname),form='unformatted',access='stream')
            read(un) forcing1, nx1, ny1, nz1
            close(un)
        end if

        call MPI_BCAST(forcing1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)

        call MPI_BCAST(nx1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(ny1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(nz1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)

        if (forcing /= forcing1) then
            write(out,*) 'fieldio_read: Error, forcing is different.'
            write(out,*) 'fieldio_read:     .in file: ', forcing
            write(out,*) 'fieldio_read: Restart file: ', forcing1
            flush(out)
            error stop
        end if 

        if (nx/=nx1 .or. ny/=ny1 .or. nz/=nz1) then
            write(out,*) 'fieldio_read: Eror, grid is different.'
            write(out,*) 'fieldio_read:     .in file: ',nx,ny,nz
            write(out,*) 'fieldio_read: Restart file: ',nx1,ny1,nz1
            flush(out)
            error stop
        end if

        ! opening the file
        call MPI_INFO_CREATE(mpi_info_var, mpi_err)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(fname), MPI_MODE_RDONLY, mpi_info_var, &
                           fh, mpi_err)

        ! In all cases x was made to be the slowest index while saving to disk
        ! And zeroed modes in x and z were skipped

        ! read only kx=0 for u

        if (ix_zero /= -1) then
            offset = 68
            count = ny_half * (nz-1)
            call MPI_FILE_READ_AT(fh, offset, buf2_sfieldk, count, MPI_COMPLEX16, &
                                mpi_status_var, mpi_err)

            do iz = 1, nz; if(iz == iz_max) cycle; do iy = 1, ny_half;

                if (iz < iz_max) then
                    vfieldk(ix_zero,iy,iz,1) = buf2_sfieldk(iy,iz,1)
                else
                    vfieldk(ix_zero,iy,iz,1) = buf2_sfieldk(iy,iz-1,1)
                end if

            end do; end do;
        end if

        ! read all of v and w
        do n = 2, 3
            
            offset = 68 + ny_half*(nz-1)*16 + (n-2)*(nx-1)*ny_half*(nz-1)*16 &
                                                + my_id*nx_perproc*ny_half*(nz-1)*16
            if (my_id > my_id_ix_max) offset = offset - ny_half*(nz-1)*16
            if (ix_max /= -1) then
                count = (nx_perproc-1) * ny_half * (nz-1)
            else
                count = nx_perproc * ny_half * (nz-1)
            end if
            call MPI_FILE_READ_AT_ALL(fh, offset, buf_sfieldk, count, MPI_COMPLEX16, &
                                mpi_status_var, mpi_err)

            _loop_spec_begin
                if (ix_max /= -1) then
                    if (ix < ix_max) then
                        if (iz < iz_max) then
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz,ix)
                        else
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix)
                        end if
                    else
                        if (iz < iz_max) then
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz,ix-1)
                        else
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix-1)
                        end if
                    end if
                else 
                    if (iz < iz_max) then
                        vfieldk(ix,iy,iz,n) =  buf_sfieldk(iy,iz,ix)
                    else
                        vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix)
                    end if
                end if
            _loop_spec_end

        end do

        ! compute the rest from kx ux + ky uy + kz uz = 0
        _loop_spec_begin
            if (ix_zero /= -1 .and. ix==ix_zero) cycle
            vfieldk(ix,iy,iz,1) = -(vfieldk(ix,iy,iz,2)*ky(iy) &
                                        + vfieldk(ix,iy,iz,3)*kz(iz)) / kx(ix)
        _loop_spec_end

        deallocate(buf_sfieldk)
        if (ix_zero /= -1) deallocate(buf2_sfieldk)

        if (ix_max /= -1) vfieldk(ix_max,:,:,1:3) = 0
        vfieldk(:,:,iz_max,1:3) = 0

        call MPI_FILE_CLOSE(fh, mpi_err)
        call MPI_INFO_FREE(mpi_info_var, mpi_err)


    end subroutine fieldio_read_xcompact

!==============================================================================

    subroutine fieldio_read_ycompact(vfieldk)

        complex(dpc), intent(out) :: vfieldk(:, :, :, :)
        _indices
        integer(i4) :: forcing1, nx1, ny1, nz1, un, count
        TYPE(MPI_File) :: fh
        integer(MPI_OFFSET_KIND)  :: offset
        complex(dpc), allocatable, dimension(:,:,:) :: &
            buf_sfieldk, buf2_sfieldk

        inquire(file=TRIM(fname),exist=there)
        if(.not.there) then
            write(out,*) 'fieldio_read: Stopping, cannot find file : '//TRIM(fname)
            flush(out)
            error stop
        end if

        if (ix_max /= -1) then
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc - 1))
            allocate(buf2_sfieldk(1, nz - 1, nx_perproc - 1))
        else
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc))
            allocate(buf2_sfieldk(1, nz - 1, nx_perproc))
        end if

        if (my_id==0) then
            open(newunit=un,file=TRIM(fname),form='unformatted',access='stream')
            read(un) forcing1, nx1, ny1, nz1
            close(un)
        end if

        call MPI_BCAST(forcing1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)

        call MPI_BCAST(nx1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(ny1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(nz1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)

        if (forcing /= forcing1) then
            write(out,*) 'fieldio_read: Error, forcing is different.'
            write(out,*) 'fieldio_read:     .in file: ', forcing
            write(out,*) 'fieldio_read: Restart file: ', forcing1
            flush(out)
            error stop
        end if 

        if (nx/=nx1 .or. ny/=ny1 .or. nz/=nz1) then
            write(out,*) 'fieldio_read: Eror, grid is different.'
            write(out,*) 'fieldio_read:     .in file: ',nx,ny,nz
            write(out,*) 'fieldio_read: Restart file: ',nx1,ny1,nz1
            flush(out)
            error stop
        end if

        ! opening the file
        call MPI_INFO_CREATE(mpi_info_var, mpi_err)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(fname), MPI_MODE_RDONLY, mpi_info_var, &
                           fh, mpi_err)

        ! In all cases x was made to be the slowest index while saving to disk
        ! And zeroed modes in x and z were skipped

        ! read all of u
        offset = 68 + my_id*nx_perproc*ny_half*(nz-1)*16
        if (my_id > my_id_ix_max) offset = offset - ny_half*(nz-1)*16
        if (ix_max /= -1) then
            count = (nx_perproc-1) * ny_half * (nz-1)
        else
            count = nx_perproc * ny_half * (nz-1)
        end if
        call MPI_FILE_READ_AT_ALL(fh, offset, buf_sfieldk, count, &
                                    MPI_COMPLEX16, mpi_status_var, mpi_err)

        _loop_spec_begin

            if (ix_max /= -1) then
                if (ix < ix_max) then
                    if (iz < iz_max) then
                        vfieldk(ix,iy,iz,1) = buf_sfieldk(iy,iz,ix)
                    else
                        vfieldk(ix,iy,iz,1) = buf_sfieldk(iy,iz-1,ix)
                    end if
                else
                    if (iz < iz_max) then
                        vfieldk(ix,iy,iz,1) = buf_sfieldk(iy,iz,ix-1)
                    else
                        vfieldk(ix,iy,iz,1) = buf_sfieldk(iy,iz-1,ix-1)
                    end if
                end if
            else 
                if (iz < iz_max) then
                    vfieldk(ix,iy,iz,1) = buf_sfieldk(iy,iz,ix)
                else
                    vfieldk(ix,iy,iz,1) = buf_sfieldk(iy,iz-1,ix)
                end if
            end if
        _loop_spec_end

        ! read only ky=0 for v
        offset = 68 + (nx-1)*ny_half*(nz-1)*16 + my_id*nx_perproc*(nz-1)*16
        if (my_id > my_id_ix_max) offset = offset - (nz-1)*16
        if (ix_max /= -1) then
            count = (nx_perproc-1) * (nz-1)
        else
            count = nx_perproc * (nz-1)
        end if
        call MPI_FILE_READ_AT_ALL(fh, offset, buf2_sfieldk, count, &
                                    MPI_COMPLEX16, mpi_status_var, mpi_err)

        do iz = 1, nz; if(iz == iz_max) cycle; do ix = 1, nx_perproc; &
                                            if(ix_max /= -1 .and. ix == ix_max) cycle;

            if (ix_max /= -1) then
                if (ix < ix_max) then
                    if (iz < iz_max) then
                        vfieldk(ix,1,iz,2) = buf2_sfieldk(1,iz,ix)
                    else
                        vfieldk(ix,1,iz,2) = buf2_sfieldk(1,iz-1,ix)
                    end if
                else
                    if (iz < iz_max) then
                        vfieldk(ix,1,iz,2) = buf2_sfieldk(1,iz,ix-1)
                    else
                        vfieldk(ix,1,iz,2) = buf2_sfieldk(1,iz-1,ix-1)
                    end if
                end if
            else 
                if (iz < iz_max) then
                    vfieldk(ix,1,iz,2) = buf2_sfieldk(1,iz,ix)
                else
                    vfieldk(ix,1,iz,2) = buf2_sfieldk(1,iz-1,ix)
                end if
            end if
        end do; end do;

        ! read all of w
        offset = 68 + (nx-1)*(ny_half+1)*(nz-1)*16 &
                                                + my_id*nx_perproc*ny_half*(nz-1)*16
        if (my_id > my_id_ix_max) offset = offset - ny_half*(nz-1)*16
        if (ix_max /= -1) then
            count = (nx_perproc-1) * ny_half * (nz-1)
        else
            count = nx_perproc * ny_half * (nz-1)
        end if
        call MPI_FILE_READ_AT_ALL(fh, offset, buf_sfieldk, count, &
                                    MPI_COMPLEX16, mpi_status_var, mpi_err)

        _loop_spec_begin

            if (ix_max /= -1) then
                if (ix < ix_max) then
                    if (iz < iz_max) then
                        vfieldk(ix,iy,iz,3) = buf_sfieldk(iy,iz,ix)
                    else
                        vfieldk(ix,iy,iz,3) = buf_sfieldk(iy,iz-1,ix)
                    end if
                else
                    if (iz < iz_max) then
                        vfieldk(ix,iy,iz,3) = buf_sfieldk(iy,iz,ix-1)
                    else
                        vfieldk(ix,iy,iz,3) = buf_sfieldk(iy,iz-1,ix-1)
                    end if
                end if
            else 
                if (iz < iz_max) then
                    vfieldk(ix,iy,iz,3) = buf_sfieldk(iy,iz,ix)
                else
                    vfieldk(ix,iy,iz,3) = buf_sfieldk(iy,iz-1,ix)
                end if
            end if
        _loop_spec_end

        ! compute the rest from kx ux + ky uy + kz uz = 0
        _loop_spec_begin
            if (iy==1) cycle
            vfieldk(ix,iy,iz,2) = -(vfieldk(ix,iy,iz,1)*kx(ix) &
                                        + vfieldk(ix,iy,iz,3)*kz(iz)) / ky(iy)
        _loop_spec_end

        deallocate(buf_sfieldk)
        deallocate(buf2_sfieldk)

        if (ix_max /= -1) vfieldk(ix_max,:,:,1:3) = 0
        vfieldk(:,:,iz_max,1:3) = 0

        call MPI_FILE_CLOSE(fh, mpi_err)
        call MPI_INFO_FREE(mpi_info_var, mpi_err)


    end subroutine fieldio_read_ycompact

!==============================================================================

    subroutine fieldio_read_zcompact(vfieldk)

        complex(dpc), intent(out) :: vfieldk(:, :, :, :)
        _indices
        integer(i4) :: forcing1, nx1, ny1, nz1, n, un, count
        TYPE(MPI_File) :: fh
        integer(MPI_OFFSET_KIND)  :: offset
        complex(dpc), allocatable, dimension(:, :, :) :: &
            buf_sfieldk, buf2_sfieldk

        inquire(file=TRIM(fname),exist=there)
        if(.not.there) then
            write(out,*) 'fieldio_read: Stopping, cannot find file : '//TRIM(fname)
            flush(out)
            error stop
        end if

        if (ix_max /= -1) then
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc - 1))
            allocate(buf2_sfieldk(ny_half, 1, nx_perproc - 1))
        else
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc))
            allocate(buf2_sfieldk(ny_half, 1, nx_perproc))
        end if

        if (my_id==0) then
            open(newunit=un,file=TRIM(fname),form='unformatted',access='stream')
            read(un) forcing1, nx1, ny1, nz1
            close(un)
        end if

        call MPI_BCAST(forcing1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)

        call MPI_BCAST(nx1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(ny1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(nz1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)

        if (forcing /= forcing1) then
            write(out,*) 'fieldio_read: Error, forcing is different.'
            write(out,*) 'fieldio_read:     .in file: ', forcing
            write(out,*) 'fieldio_read: Restart file: ', forcing1
            flush(out)
            error stop
        end if 

        if (nx/=nx1 .or. ny/=ny1 .or. nz/=nz1) then
            write(out,*) 'fieldio_read: Eror, grid is different.'
            write(out,*) 'fieldio_read:     .in file: ',nx,ny,nz
            write(out,*) 'fieldio_read: Restart file: ',nx1,ny1,nz1
            flush(out)
            error stop
        end if

        ! opening the file
        call MPI_INFO_CREATE(mpi_info_var, mpi_err)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(fname), MPI_MODE_RDONLY, mpi_info_var, &
                           fh, mpi_err)

        ! In all cases x was made to be the slowest index while saving to disk
        ! And zeroed modes in x and z were skipped

        ! read all of u and v
        do n = 1, 2

            offset = 68 + (n-1)*(nx-1)*ny_half*(nz-1)*16 &
                                                + my_id*nx_perproc*ny_half*(nz-1)*16
            if (my_id > my_id_ix_max) offset = offset - ny_half*(nz-1)*16
            if (ix_max /= -1) then
                count = (nx_perproc-1) * ny_half * (nz-1)
            else
                count = nx_perproc * ny_half * (nz-1)
            end if
            call MPI_FILE_READ_AT_ALL(fh, offset, buf_sfieldk, count, &
                                        MPI_COMPLEX16, mpi_status_var, mpi_err)

            _loop_spec_begin
                if (ix_max /= -1) then
                    if (ix < ix_max) then
                        if (iz < iz_max) then
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz,ix)
                        else
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix)
                        end if
                    else
                        if (iz < iz_max) then
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz,ix-1)
                        else
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix-1)
                        end if
                    end if
                else 
                    if (iz < iz_max) then
                        vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz,ix)
                    else
                        vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix)
                    end if
                end if
            _loop_spec_end

        end do

        ! read only kz=0 for w
        offset = 68 + 2*(nx-1)*ny_half*(nz-1)*16 + my_id*nx_perproc*ny_half*16
        if (my_id > my_id_ix_max) offset = offset - ny_half*16
        if (ix_max /= -1) then
            count = (nx_perproc-1) * ny_half
        else
            count = nx_perproc * ny_half
        end if
        call MPI_FILE_READ_AT_ALL(fh, offset, buf2_sfieldk, count, &
                                    MPI_COMPLEX16, mpi_status_var, mpi_err)

        do iy = 1, ny_half; do ix = 1, nx_perproc; if(ix_max /= -1 .and. ix == ix_max) cycle;

            if (ix_max /= -1) then
                if (ix < ix_max) then
                    vfieldk(ix,iy,1,3) = buf2_sfieldk(iy,1,ix)
                else
                    vfieldk(ix,iy,1,3) = buf2_sfieldk(iy,1,ix-1)
                end if
            else 
                vfieldk(ix,iy,1,3) = buf2_sfieldk(iy,1,ix)
            end if
        end do; end do;

        ! compute the rest from kx ux + ky uy + kz uz = 0
        _loop_spec_begin
            if (iz==1) cycle
            vfieldk(ix,iy,iz,3) = -(vfieldk(ix,iy,iz,1)*kx(ix) &
                                            + vfieldk(ix,iy,iz,2)*ky(iy)) / kz(iz)
        _loop_spec_end

        deallocate(buf_sfieldk)
        deallocate(buf2_sfieldk)

        if (ix_max /= -1) vfieldk(ix_max,:,:,1:3) = 0
        vfieldk(:,:,iz_max,1:3) = 0

        call MPI_FILE_CLOSE(fh, mpi_err)
        call MPI_INFO_FREE(mpi_info_var, mpi_err)


    end subroutine fieldio_read_zcompact

!==============================================================================

    subroutine fieldio_read_nocompact(vfieldk)

        complex(dpc), intent(out) :: vfieldk(:, :, :, :)
        _indices
        integer(i4) :: forcing1, nx1, ny1, nz1, n, un, count
        TYPE(MPI_File) :: fh
        integer(MPI_OFFSET_KIND)  :: offset
        complex(dpc), allocatable, dimension(:, :, :) :: &
            buf_sfieldk, buf2_sfieldk

        inquire(file=TRIM(fname),exist=there)
        if(.not.there) then
            write(out,*) 'fieldio_read: Stopping, cannot find file : '//TRIM(fname)
            flush(out)
            error stop
        end if

        if (ix_max /= -1) then
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc - 1))
            allocate(buf2_sfieldk(ny_half, 1, nx_perproc - 1))
        else
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc))
            allocate(buf2_sfieldk(ny_half, 1, nx_perproc))
        end if

        if (my_id==0) then
            open(newunit=un,file=TRIM(fname),form='unformatted',access='stream')
            read(un) forcing1, nx1, ny1, nz1
            close(un)
        end if

        call MPI_BCAST(forcing1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)

        call MPI_BCAST(nx1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(ny1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(nz1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)

        ! if (forcing /= forcing1) then
        !     write(out,*) 'fieldio_read: Error, forcing is different.'
        !     write(out,*) 'fieldio_read:     .in file: ', forcing
        !     write(out,*) 'fieldio_read: Restart file: ', forcing1
        !     flush(out)
        !     error stop
        ! end if 

        if (nx/=nx1 .or. ny/=ny1 .or. nz/=nz1) then
            write(out,*) 'fieldio_read: Eror, grid is different.'
            write(out,*) 'fieldio_read:     .in file: ',nx,ny,nz
            write(out,*) 'fieldio_read: Restart file: ',nx1,ny1,nz1
            flush(out)
            error stop
        end if

        ! opening the file
        call MPI_INFO_CREATE(mpi_info_var, mpi_err)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(fname), MPI_MODE_RDONLY, mpi_info_var, &
                           fh, mpi_err)

        ! In all cases x was made to be the slowest index while saving to disk
        ! And zeroed modes in x and z were skipped

        do n = 1, 3

            offset = 68 + (n-1)*(nx-1)*ny_half*(nz-1)*16 &
                                                + my_id*nx_perproc*ny_half*(nz-1)*16
            if (my_id > my_id_ix_max) offset = offset - ny_half*(nz-1)*16
            if (ix_max /= -1) then
                count = (nx_perproc-1) * ny_half * (nz-1)
            else
                count = nx_perproc * ny_half * (nz-1)
            end if
            call MPI_FILE_READ_AT_ALL(fh, offset, buf_sfieldk, count, &
                                        MPI_COMPLEX16, mpi_status_var, mpi_err)

            _loop_spec_begin
                if (ix_max /= -1) then
                    if (ix < ix_max) then
                        if (iz < iz_max) then
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz,ix)
                        else
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix)
                        end if
                    else
                        if (iz < iz_max) then
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz,ix-1)
                        else
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix-1)
                        end if
                    end if
                else 
                    if (iz < iz_max) then
                        vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz,ix)
                    else
                        vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix)
                    end if
                end if
            _loop_spec_end

        end do

        deallocate(buf_sfieldk)
        deallocate(buf2_sfieldk)

        if (ix_max /= -1) vfieldk(ix_max,:,:,1:3) = 0
        vfieldk(:,:,iz_max,1:3) = 0

        call MPI_FILE_CLOSE(fh, mpi_err)
        call MPI_INFO_FREE(mpi_info_var, mpi_err)


    end subroutine fieldio_read_nocompact

!==============================================================================

    subroutine fieldio_read(vfieldk)

        complex(dpc), intent(out) :: vfieldk(:, :, :, :)

        if (nx - 1 >= ny_half .and. nx - 1 >= nz - 1) then
            call fieldio_read_xcompact(vfieldk)
        elseif (ny_half >= nx - 1 .and. ny_half >= nz - 1) then
            call fieldio_read_ycompact(vfieldk)
        else
            call fieldio_read_zcompact(vfieldk)
        end if
    end subroutine fieldio_read

!==============================================================================

end module fieldio

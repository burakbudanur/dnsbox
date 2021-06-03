module m_io

    use m_openmpi
    
    character*6   :: file_ext
    character*255 :: fname
    integer       :: in, out

    contains

!================================================================================

        subroutine m_io_init

            write(fname,"('d',i4.4,'.txt')") myid
            open(newunit=out,file=fname,position="append")
            write(out,"('-------------------------------------')")
            write(out,"('Process ',i4,' of ',i4,'(',i4.4,') is alive.')") &
                  myid, numprocs, numprocs-1
            flush(out)
            return
        end subroutine m_io_init

!================================================================================

        subroutine m_io_exit
            
            write(out,"('Done.')")
            close(out)
            return
        end subroutine m_io_exit
    
end module m_io

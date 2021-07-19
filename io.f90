module io
    use numbers
    use openmpi

    character(6)   :: file_ext
    character(255) :: fname
    integer(i4)    :: out
    logical        :: there, there2

    contains

!================================================================================

    subroutine io_init

        write(fname,"('d',i4.4,'.txt')") my_id
        open(newunit=out,file=fname,position="append")
        write(out,"('Process ',i4,' of ',i4,' is alive.')") &
                my_id, num_procs-1
        write(out, '(79(''=''))')
    end subroutine io_init

!================================================================================

    subroutine io_exit
        
        write(out,*) "io_exit: Done."
        close(out)
    end subroutine io_exit

!==============================================================================

end module io

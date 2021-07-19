module openmpi
    use mpi_f08
    use numbers
    
    type(MPI_Info)   :: mpi_info_var
    type(MPI_Status) :: mpi_status_var
    integer(i4)      :: mpi_err, num_procs, my_id
    
    contains 

!==============================================================================
    
    subroutine openmpi_init
        use numbers
        
        call MPI_INIT(mpi_err)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_err)
        call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, mpi_err)

    end subroutine openmpi_init
    
!==============================================================================
    
    subroutine openmpi_exit
        
        call MPI_FINALIZE(mpi_err)
    
    end subroutine openmpi_exit

!==============================================================================
    
end module openmpi

module m_openmpi
    use mpi_f08

    ! MPI variables
    logical :: iammaster
    type(MPI_Info) :: mpi_info_var
    type(MPI_Status) :: mpi_status_var
    integer :: myid, numprocs, master, mpi_err
    integer :: id_to, id_from, tag, count
    integer :: id_root
    
    contains 

!==============================================================================
    
    subroutine m_openmpi_init
    
        ! Initiate mpi
        call MPI_INIT(mpi_err)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, mpi_err)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpi_err)
    
    end subroutine m_openmpi_init
    
!==============================================================================
    
    subroutine m_openmpi_exit
    
        call MPI_FINALIZE(mpi_err)

    end subroutine m_openmpi_exit
    
end module m_openmpi

module m_stats
    use m_numbers
    use m_parameters
    use m_work
    use x_fftw


    real(dp) :: eps_v, eta, enstrophy
    real(dp) :: norm1, norm, normrhs
    real(dp) :: sctmp, prod1, prod

    ! i/o channels
    integer :: stat1fileCh
    
    character(255) :: stat1file = 'stat.gp'

    ! Will use these to check if LAMINARIZED
    real(dp) :: e_target, prod_target, diss_target
    real(dp) :: e_diff, prod_diff, diss_diff
    integer :: KILL_SWITCH
    
    logical :: statsWritten = .false.

    ! Which index holds kF
    integer :: kkF

    contains 
            
    subroutine m_stats_init
        
        integer :: k
        real(dp) :: ky
        
        e_target = (gamma ** 2) / (4 * kF**4 * nu**2)
        prod_target =  (gamma ** 2) / (2 * kF**2 * nu)
        diss_target = prod_target
        KILL_SWITCH = 0

        ! Find which index holds kF
        kkF = -1
        do k = 1, nx
            ky = aky(k)
            if (ky == kF) kkF = k
        end do
            
        return
    
    end subroutine m_stats_init

!==============================================================================

    subroutine m_stats_compute
        ! main subroutine computes and saves stats
        use m_numbers
        use m_parameters
        use m_fields
        use m_work
        use x_fftw
        
        real(dp) :: fac
        integer :: i
        
        wrk(:, :, :, 1:3) = fields(:, :, :, 1:3)

        prod    = zero
        prod1   = zero
        prod1 = -gamma * normfac * fields(2, 1, kkF, 1)
        call m_fields_energy
                
        ! Take derivatives
        
        call x_derivative(3, 'y', 6) ! w_y
        call x_derivative(3, 'x', 5) ! w_x
        
        call x_derivative(2, 'z', 4) ! v_z
        call x_derivative(2, 'x', 3) ! v_x
        
        call x_derivative(1, 'z', 2) ! u_z
        call x_derivative(1, 'y', 1) ! u_y
        
        ! Vorticity:
        wrk(:, :, :, 3) = wrk(:, :, :, 3) - wrk(:, :, :, 1) ! omg_3 = v_x - u_y
        wrk(:, :, :, 2) = wrk(:, :, :, 2) - wrk(:, :, :, 5) ! omg_2 = u_z - w_x
        wrk(:, :, :, 1) = wrk(:, :, :, 6) - wrk(:, :, :, 4) ! omg_1 = w_y - v_z
        
        ! Transform to physical space:
        do i = 1, 3
            call xFFT3d(-1, i)
            wrk(:, :, :, i) = wrk(:, :, :, i) * normfac
        end do
        
        ! Enstrophy:
        fac = one / real(nz * ny * nx_all, 8)
        sctmp = sum(wrk(1:nz, :, :, 1:3)**2) * fac
        count = 1
        call MPI_REDUCE(sctmp, enstrophy, count, MPI_REAL8, MPI_SUM, 0, &
                        MPI_COMM_WORLD, mpi_err)
        call MPI_ALLREDUCE(sctmp, enstrophy, count, MPI_REAL8, MPI_SUM, &
                        MPI_COMM_WORLD, mpi_err)
        eps_v = enstrophy * nu

        ! reducing the Kolmogorov production from each node to master node
        call MPI_ALLREDUCE(prod1, prod, 1, MPI_REAL8, MPI_SUM, &
                        MPI_COMM_WORLD, mpi_err)
                        
        
    end subroutine m_stats_compute

!==============================================================================

    subroutine m_stats_write
        use m_numbers
        use m_parameters
        use m_fields
        use m_work
        use x_fftw        
        
        logical :: there2
        
        ! outputting statistics
        
        call m_fields_divergence
        
        if (myid==0) then
        
            ! Kolmogorov scale
            eta = (nu**3/eps_v)**0.25           
            
           ! outputting all this in the stat1 file
            inquire(file=TRIM(stat1file), exist=there, opened=there2)
            if (.not.there) then
               open(newunit=stat1fileCh,file=TRIM(stat1file),form='formatted')
               write(stat1fileCh,"(A2,"//i4_len//","//"7"//sp_len//")") "# ", "itime", "time", "energy", "prod", "diss", "eta", "normrhs", "Courant"
            end if
            if(there.and..not.there2) then
               open(newunit=stat1fileCh,file=TRIM(stat1file),position='append')
            end if
            write(stat1fileCh,"(A2,"//i4_f//","//"7"//sp_f//")") "  ", itime - 1, time, energy, prod, eps_v, eta, normrhs, courant
           
           statsWritten = .true.

    end if

    end subroutine m_stats_write

!==============================================================================

    subroutine m_stats_norm
        
        ! Compute L^2 norm of wrk(:, :, :, 1:3)
        norm  = 0.0d0
        norm1 = sum(inprodCoeffs * wrk(:,:,:,1:3)**2)

        ! summing contributions on each node
        call MPI_ALLREDUCE(norm1, norm, 1, MPI_REAL8, MPI_SUM, &
                           MPI_COMM_WORLD, mpi_err)
        
        norm = dsqrt(norm)
        
    end subroutine m_stats_norm
        
end module m_stats

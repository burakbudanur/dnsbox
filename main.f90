program main
    
    ! modules:
    use m_numbers
    use m_openmpi
    use m_io
    use m_parameters
    use m_work
    use m_fields
    use x_fftw
    use m_state
    use m_stats
    use m_step    
    use m_runs
    
    integer :: un ! dummy i/o channel

    ! Time averages
    character(255) :: averagesFile = 'averages.gp'
    integer :: averagesFileCh, periodCh, ndts
    real(dp) :: period
    real(dp) :: avgKin, avgProd, avgDiss
    real(dp) :: stdKin, stdProd, stdDiss
    real(dp), allocatable :: average(:,:,:,:)
    
    ! initialization:
    call m_openmpi_init
    call m_io_init
    call m_parameters_init
    call m_work_init
    call m_fields_init
    call x_fftw_allocate(1)
    call x_fftw_init
    call m_state_init
    
    ! Initial state
    call m_runs_init

    if (IC /= 1) then
        call m_runs_impose_symmetries
        write(out, *) "imposed symmetries"
    end if

    call m_stats_init

    ! We don't want random ICs with negative production
    call m_stats_compute
    call MPI_BCAST(prod, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
    if (prod < 0 .and. IC == 0) then
        if (SEEDNO == -1) then
            write(out, *) "Random IC started with a negative production, re-rolling until positive."
            flush(out)
            do while (prod < 0)
                call m_runs_init
                call m_runs_impose_symmetries
                call m_stats_compute
                call MPI_BCAST(prod, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
            end do
        else
            write(out, *) "Random IC started with a negative production. Use another seed or set SEEDNO to -1."
            call m_runs_exit
        end if
    end if
    
    if (IC /= 1) then
        ! Save the initial condition
        write(out, *) 'Saving the initial condition'
        flush(out)
        wrk(:, :, :, 1:3) = fields(:, :, :, 1:3)
        call m_runs_write
    end if

    ! Compute and print divergence:
    call m_fields_divergence

    write(out, *) "Starting time stepping"
    flush(out)

    if (integratePeriodic == 1) then
        open(newunit=periodCh,status='old',file='guess.in')
        read(periodCh,*) period 
        read(periodCh,*) ndts 
        close(periodCh)
        if(ndts==-1) ndts = nint(period/dt)
        ITMAX = ndts
        dt = period / real(ndts, 8)
        IPRINT1 = 1 ! For accurate time averages

        ! Initialize the energetics
        avgKin = zero
        avgProd = zero
        avgDiss = zero
        stdKin = zero
        stdProd = zero
        stdDiss = zero
        allocate(average(nz+2, ny, nx, 3))
        average = 0
    end if

    do itime = ITMIN + 1, ITMAX

        ! Time stepping
        
        ! Compute stats if it is time
        if (mod(itime - 1, IPRINT1) == 0) then
            call m_stats_compute

            if (integratePeriodic == 1 .and. itime > itmin + 1) then
                ! Update the energetics
                avgKin = avgKin + dt * energy
                stdKin = stdKin + dt * energy**2
                avgProd = avgProd + dt * prod
                stdProd = stdProd + dt * prod**2
                avgDiss = avgDiss + dt * eps_v
                stdDiss = stdDiss + dt * eps_v**2

                average(:,:,:,1:3) = average(:,:,:,1:3) + fields(:,:,:,1:3) * dt
            end if

        end if
        
        call m_step_calc_rhs

        ! Write stats
        if (mod(itime - 1, IPRINT1) == 0) then
            call m_stats_write

            ! Several kill switches and checks
            if (myid==0) then
                ! Stop if LAMINARIZED
                if (LAM_SWITCH == 1) then
                    e_diff = abs(energy - e_target)/e_target
                    prod_diff = abs(prod - prod_target)/prod_target
                    diss_diff = abs(eps_v - diss_target)/diss_target
                    
                    if (e_diff < LAM_THRESHOLD .and. &
                        prod_diff < LAM_THRESHOLD .and. diss_diff < LAM_THRESHOLD) then
                        KILL_SWITCH = 1
                        write(out, *) "Laminarized, stopping."
                        open(newunit=un,file='LAMINARIZED',position='append')
                        write(un,*) time
                        close(un)
                    end if
                    
                end if
            end if

            ! Broadcast the KILL_SWITCH status
            call MPI_BCAST(KILL_SWITCH, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
            
            ! Kill the processes if signalled
            if (KILL_SWITCH == 1) then
                call m_runs_exit
            end if
            
        end if
        
        call m_fields_dealias
        call m_step_precorr
        
        call m_fields_pressure

        time = time + dt  
      
        call m_runs_impose_symmetries

        ! Dealias
        call m_fields_dealias

        ! This is after having time-stepped, so no ITIME lag here.
        if (mod(itime, IPRINT2) == 0) then
            ! write(out, *) 'Saving state file t = ', time
            wrk(:, :, :, 1:3) = fields(:, :, :, 1:3)
            call m_runs_write
            flush(out)
            if (statsWritten) flush(stat1fileCh)
        end if
        
        ! Flush everything at the first time step
        if (myid == 0 .and. itime == ITMIN + 1) then
            flush(out)
            flush(stat1fileCh)
        end if
    end do
    
    write(out, *) "Finished time stepping. Run complete."
    flush(out)

    ! One last time compute and write things if the next step would have it.
    ! Because otherwise the very last state doesn't have anything recorded
    ! except for its state file.

    call m_step_calc_rhs
    if (mod(itime, IPRINT1) == 0) then
        call m_stats_compute
        if (integratePeriodic == 1) then
            ! Update the energetics
            avgKin = avgKin + dt * energy
            stdKin = stdKin + dt * energy**2
            avgProd = avgProd + dt * prod
            stdProd = stdProd + dt * prod**2
            avgDiss = avgDiss + dt * eps_v
            stdDiss = stdDiss + dt * eps_v**2
            average(:,:,:,1:3) = average(:,:,:,1:3) + fields(:,:,:,1:3) * dt
        end if
    end if

    if (mod(itime, IPRINT1) == 0) then
        call m_stats_write
    end if

    ! Compute the averages
    if (integratePeriodic == 1) then

        avgKin = avgKin / time
        stdKin = DSQRT(stdKin / time - avgKin ** 2)
        avgProd = avgProd / time
        stdProd = DSQRT(stdProd / time - avgProd ** 2)
        avgDiss = avgDiss / time
        stdDiss = DSQRT(stdDiss / time - avgDiss ** 2)

        if (myid == 0) then
            open(newunit=averagesFileCh,file=TRIM(averagesFile),form='formatted',status='replace')
            write(averagesFileCh,"(A2,"//"6"//dp_len//")") "# ", "avgKin", "stdKin", "avgProd", "stdProd", "avgDiss", "stdDiss"
            write(averagesFileCh,"(A2,"//"6"//dp_f//")") "  ", avgKin, stdKin, avgProd, stdProd, avgDiss, stdDiss
            close(averagesFileCh)
        end if

        average(:,:,:,1:3) = average(:,:,:,1:3) / time
        wrk(:,:,:,1:3) = average(:,:,:,1:3)
        write(file_ext, "(i6.6)") 0
        fname = 'average.'//file_ext
        call m_work_write
    end if
     
    call m_runs_exit
    
end program main


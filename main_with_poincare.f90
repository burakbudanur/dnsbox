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
    use m_pod
    
    implicit none
    
    ! Test vars
    real(dp)  :: myErrorMax, ErrorMax
    real(dp)  :: xx, yy, zz 
    
    ! Final error counters
    integer :: i, j, k, n
    integer :: un ! dummy i/o channel
    ! Kill switches
    logical :: killswitchthere

    ! Post-processing logic
    logical :: forbidDNS

    ! Time averages
    character(255) :: averagesFile = 'averages.gp'
    integer :: averagesFileCh, periodCh, ndts
    real(dp) :: period
    real(dp) :: avgKin, avgProd, avgDiss
    real(dp) :: stdKin, stdProd, stdDiss

    ! Time res increase
    logical :: timeResIncreaseThere
    integer :: timeResSteps
    real(dp)  :: timeOld
    inquire(file='timeResIncrease', exist=timeResIncreaseThere)
    if (timeResIncreaseThere) then
        ! Read visualizeMany.in
        open(newunit=un,status='old',file='timeResIncrease')
        read(un, *) timeResSteps
        close(un)
    end if

    adjoint = .false.
    
    ! initialization:
    call m_openmpi_init
    call m_io_init
    call m_parameters_init
    call m_work_init
    call m_fields_init
    call x_fftw_allocate(1)
    call x_fftw_init
    call m_state_init
    call m_step_init

    if (podSwitch == 1) then
        call m_pod_init
    end if
    
    ! Initial state
    call m_runs_init

    inquire(file='NEWTON_NOT_CONVERGED', exist=forbidDNS)
    inquire(file='EIGEN_NOT_CONVERGED', exist=forbidDNS)
    if (forbidDNS) then
        write(out, *) "Refuse to start from a non-converged Newton or Arnoldi search."
        call m_runs_exit
    end if

    ! Impose symmetries
    call m_runs_impose_symmetries
    write(out, *) "imposed symmetries"

    call m_stats_init

    ! We don't want random ICs with negative production
    call m_stats_compute
    call MPI_BCAST(prod, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
    if (prod < 0 .and. IC == 0) then
        if (SEEDNO == -1) then
            write(out, *) "Random IC started with a negative production, rolling another one."
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
    
    ! Compute and print divergence:
    call m_fields_divergence

    write(out, *) "Starting time stepping"
    flush(out)
    
    if (tAdapt == 1) then 
        
        write(out, *) "Setting the time step"
        flush(out)
            
        if (mod(itime, IPRINT1) == 0) then
            call m_step_calc_rhs
        end if
        
        ! Set the time-step
        call m_step_set_dt
        
    end if

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
            end if

        end if

        if (PoincarePD == 1) then
            ! Store the current state for possibly computing
            ! the U(t*) = Dissipation(t*) - Production(t*) = 0 intersection
            ! With directional constraint d_t U(t=t*) > 0
         
            if ( not(mod(itime - 1, IPRINT1)) == 0 ) then
                ! compute stats if haven't already
                call m_stats_compute
            end if

            UPoincare = eps_v - prod
            ! Save the field for possibly computing the Poincare
            ! section intersection
            fields(:, :, :, 25:27) = fields(:, :, :, 1:3)

        end if 

        ! Project if it is time
        if (podSwitch == 1. .and. mod(itime - 1, IPROJECT) == 0) then
            call m_pod_project 
        end if
        
        ! Calculate rhs. If using precorr, we need to know that
        !    always.
        if (tStepper == 4 .or. mod(itime - 1, IPRINT1) == 0) call m_step_calc_rhs
        
        if (mod(itime - 1, IPRINT3) == 0) call m_stats_write_spec
        ! if (podSwitch == 0 .and. mod(itime - 1, IPRINT3) == 0) call m_stats_write_spec
        
        ! Write projections if it is time
        if (podSwitch == 1. .and. mod(itime - 1, IPROJECT) == 0) then
            call m_pod_projections_write
        end if

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
                
                ! Energy and dissipation are non-negative
                if (energy < 0 .or. eps_v < 0) then
                    KILL_SWITCH = 1
                    write(out, *) "Energy or dissipation negative, stopping."
                    open(newunit=un,file='NEGATIVE_ENERGY_DISS',position='append')
                    write(un,*) time
                    close(un)
                end if
                
                ! It's interesting to have negative production during a run
                if (prod < 0) then
                    ! write(out, *) "Production is negative."
                    open(newunit=un,file='NEGATIVE_PROD_DURING',position='append')
                    write(un,*) time
                    close(un)
                    ! If the killswitch fle is there, stop the run
                    inquire(file='STOPNEGDURING', exist=killswitchthere)
                    if (killswitchthere) KILL_SWITCH = 1
                end if
            end if

            ! Broadcast the KILL_SWITCH status
            call MPI_BCAST(KILL_SWITCH, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
            
            ! Kill the processes if signalled
            if (KILL_SWITCH == 1) then
                call m_runs_exit
            end if
            
        end if
        
        if (tStepper == 1) then 
        
            call m_step_rkg
        
        else if (tStepper == 2) then
            
            call m_step_etdheun
        
        else if (tStepper == 3) then
            
            call m_step_dopri
        
        
        else if (tStepper == 4) then
            call m_fields_dealias
            call m_step_precorr
        
        end if
        
        call m_fields_pressure


        if (PoincarePD == 1) then
            ! Check if there is a Poincare section intersection         
            ! compute stats 
            call m_stats_compute
            UPoincareNext = eps_v - prod
            
            if(UPoincare < 0.0d0 .and. UPoincareNext > 0.0d0) then


                write(out, *) 'There is a Poincare section intersection'
                write(out, *) 'UPoincare = ', UPoincare
                write(out, *) 'UPoincareNext = ', UPoincareNext
                flush(out)

                ! There is an intersection
                dthold    = dt     ! Save actual time step
                fields(:, :, :, 28:30) = fields(:, :, :, 1:3) ! save current state
                
                if (dabs(UPoincare) < epsPoincare) then 
                    ! previous state satisfies Poincare section condition
                    wrk(:, :, :, 1:3) = fields(:, :, :, 25:27)
                else if (dabs(UpoincareNext) < epsPoincare) then
                    ! current state satisfies Poincare section condition
                    wrk(:, :, :, 1:3) = fields(:, :, :, 1:3)
                else
                    ! intersecting state is inbetween
                    
                    errPoincare = dmin1(Upoincare, UpoincareNext)
                    deltatPoinc = dthold
                    fields(:, :, :, 1:3) = fields(:, :, :, 25:27) ! previous tstep

                    do while(dabs(errPoincare) < epsPoincare)
                        
                        fields(:, :, :, 25:27) = fields(:, :, :, 1:3) ! repeat state back up
                        ! starting the secant loop
                        ! guess the time-step
                        dt = deltatPoinc / (1 - UPoincareNext / UPoincare) ! guessing time step
                        
                        ! Take a time step
                        if (tStepper == 1) then 
        
                            call m_step_rkg
                        
                        else if (tStepper == 2) then
                            
                            call m_step_etdheun
                        
                        else if (tStepper == 3) then
                            
                            call m_step_dopri
                        
                        
                        else if (tStepper == 4) then
                            
                            call m_step_precorr
                        
                        end if

                        ! compute production dissipation
                        call m_stats_compute
                        errPoincare = eps_v - prod
                        if (dabs(errPoincare) < epsPoincare) then
                            ! copy the state for saving
                            call m_fields_dealias
                            wrk(:, :, :, 1:3) = fields(:, :, :, 1:3)
                        else if (errPoincare < 0.0d0) then
                            UPoincare     = errPoincare
                            deltatPoinc   = deltatPoinc - dt
                        else 
                            UPoincareNext = errPoincare
                            deltatPoinc   = dt
                            fields(:, :, :, 1:3) = fields(:, :, :, 25:27) ! Restore prev. state                            
                        end if

                    end do

                end if
                
                ! write intersection
                call m_runs_write_poincare
                iPoincare = iPoincare + 1
                ! restore the current state and the time step
                fields(:, :, :, 1:3) = fields(:, :, :, 28:30)
                dt = dthold

            end if
            
        end if 


        if (.not. (tStepper == 3)) time = time + dt  
      
        call m_runs_impose_symmetries

        if (timeResIncreaseThere .and. mod(itime, IPRINT2) == 0 .and. mod(itime/IPRINT2, timeResSteps) == 0) then
            timeOld = time
            call m_runs_read
            fields(:, :, :, 1:3) = wrk(:, :, :, 1:3)
            time = timeOld
        end if

        ! Dealias
        call m_fields_dealias

        ! This is after having time-stepped, so no ITIME lag here.
        if (fullSaveDisabled == 0 .and. mod(itime, IPRINT2) == 0) then
            ! write(out, *) 'Saving state file t = ', time
            wrk(:, :, :, 1:3) = fields(:, :, :, 1:3)
            call m_runs_write
        end if
        
        ! Flush everything at the first time step
        if (myid == 0 .and. itime == ITMIN + 1) then
            flush(out)
            flush(stat1fileCh)
            flush(stat2fileCh)
            flush(esFileCh)
            if (podSwitch == 1) flush(projCh)
        end if
    end do
    
    write(out, *) "Finished time stepping. Run complete."
    flush(out)

    ! One last time compute and write things if the next step would have it.
    ! Because otherwise the very last state doesn't have anything recorded
    ! except for its state file.

    if (tStepper == 4 .or. mod(itime, IPRINT1) == 0) call m_step_calc_rhs
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
        end if
    end if
    if (mod(itime, IPRINT3) == 0) call m_stats_write_spec
    if (mod(itime, IPRINT1) == 0) then
        call m_stats_write
    end if
    if (podSwitch == 1. .and. mod(itime, IPROJECT) == 0) then
        call m_pod_project 
        call m_pod_projections_write
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
    end if

    if (IC == -2) then
      
        write(out, *) "Computing final error"
        flush(out)
  
        ! Check the final error:
        wrk(:, :, :, 1:3) = fields(:, :, :, 1:3)
        do n = 1, 3
            call xFFT3d(-1, n)
            wrk(:, :, :, n) = wrk(:, :, :, n) * normfac
        end do    
        
            
        do k = 1, nz
            xx = real(nz * myid + k - 1, 8) * dx
            do j = 1, ny
                yy = real(j - 1, 8) * dy
                do i = 1, nx
                    zz = real(i - 1, 8) * dz
                    
                    wrk(i, j, k, 4) = wrk(i, j, k, 1) &
                                    + (real(0.5d0, 8) &
                             * (dsqrt(3.0d0) * dcos(xx) * dsin(yy) * dsin(zz) & 
                              + dsin(xx) * dcos(yy) * dcos(zz))) &
                             * dexp(- 3.0d0 * time * nu) 

                    wrk(i, j, k, 5) = wrk(i, j, k, 2) &  
                                    - real(0.5d0, 8) &
                             * (dsqrt(3.0d0) * dsin(xx) * dcos(yy) * dsin(zz) & 
                              - dcos(xx) * dsin(yy) * dcos(zz)) &
                             * dexp(- 3.0d0 * time * nu) 
                                                                          
                    wrk(i, j, k, 6) = wrk(i, j, k, 3) &
                                    - (dcos(xx) * dcos(yy) * dsin(zz) &
                                     * dexp(- 3.0d0 * time * nu))  
                    
                end do
            end do
        end do
        
        wrk(:, :, :, 1) = dabs(wrk(:, :, :, 4)) &
                        + dabs(wrk(:, :, :, 5)) &
                        + dabs(wrk(:, :, :, 6))
                        
        myErrorMax = maxval(wrk(:, :, :, 1))
        call MPI_REDUCE(myErrorMax, ErrorMax, 1, MPI_REAL8, MPI_MAX, 0, &
                        MPI_COMM_WORLD, mpi_err)
        
        if (myid == 0) then
            
            write(out, *) 'Maximum error after time stepping: ', ErrorMax
            
        end if    
    
    end if
     
    call m_runs_exit

    contains
    
end program main


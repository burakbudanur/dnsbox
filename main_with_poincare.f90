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

call m_step_precorr

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
                
                call m_step_precorr

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

time = time + dt  

call m_runs_impose_symmetries

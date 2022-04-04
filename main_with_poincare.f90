if (adaptive_dt .or. (i_print_steps > 0 .and. mod(itime, i_print_steps) == 0)) then
    call timestep_courant(vel_vfieldxx_now)
end if

if (adaptive_dt) call timestep_set_dt

if (poincare) then
    ! Store the current state for possibly computing
    ! the U(t*) = Dissipation(t*) - Production(t*) = 0 intersection
    ! With directional constraint d_t U(t=t*) > 0
    
    if ((i_print_stats > 0 .and. .not. (mod(itime, i_print_stats) == 0)) &
            .or. i_print_stats <= 0) then
        call stats_compute(vel_vfieldk_now, fvel_vfieldk_now)
    end if

    U_poincare = dissip - powerin
    ! Save the field for possibly computing the Poincare
    ! section intersection
    vel_vfieldk_before = vel_vfieldk_now
    vel_vfieldxx_before = vel_vfieldxx_now
    fvel_vfieldk_before = fvel_vfieldk_now
    dt_before = dt
    time_before = time
end if

call timestep_precorr(vel_vfieldxx_now, vel_vfieldk_now, fvel_vfieldk_now)
time = time + dt
itime = itime + 1

if (poincare) then
    ! Check if there is a Poincare section intersection         
    call stats_compute(vel_vfieldk_now, fvel_vfieldk_now)
    U_poincare_next = dissip - powerin
    
    ! If there is an intersection
    if(U_poincare < 0 .and. U_poincare_next > 0) then

        ! backup timestepped state
        vel_vfieldk_after = vel_vfieldk_now
        vel_vfieldxx_after = vel_vfieldxx_now
        fvel_vfieldk_after = fvel_vfieldk_now
        dt_after = dt
        time_after = time
        itime_after = itime
        
        write(file_ext, "(i6.6)") i_poincare
        fname = 'poincare.'//file_ext
        itime = i_poincare
        if (abs(U_poincare) < eps_poincare) then 
            ! previous state satisfies Poincare section condition
            time = time_before
            dt = dt_before
            call fieldio_write(vel_vfieldk_before)
        else if (abs(U_poincare_next) < eps_poincare) then
            ! current state satisfies Poincare section condition
            call fieldio_write(vel_vfieldk_now)
        else
            ! intersecting state is inbetween

            ! restore to the state before the section
            vel_vfieldk_now = vel_vfieldk_before
            vel_vfieldxx_now = vel_vfieldxx_before
            dt = dt_before
            
            ! starting the secant loop
            dt_last = dt
            do
                ! guess the time-step
                dt = dt_last / (1 - U_poincare_next / U_poincare)
                call timestep_precorr(vel_vfieldxx_now, vel_vfieldk_now, fvel_vfieldk_now)
                time = time_before + dt
                call stats_compute(vel_vfieldk_now, fvel_vfieldk_now)
                U_poincare_next = dissip - powerin

                ! sanity check
                if (U_poincare_next < -abs(U_poincare)) then
                    write(out, *) "Poincare: Secant method failed."
                    flush(out)
                    call run_exit
                end if

                if (abs(U_poincare_next) < eps_poincare) then
                    ! we're done
                    call fieldio_write(vel_vfieldk_now)
                    exit
                else
                    ! restore to before the section
                    vel_vfieldk_now = vel_vfieldk_before
                    vel_vfieldxx_now = vel_vfieldxx_before
                    fvel_vfieldk_now = fvel_vfieldk_before
                    time = time_before
                    ! update dt guess
                    dt_last = dt
                end if

            end do

        ! restore the timestepped state
        vel_vfieldk_now = vel_vfieldk_after
        vel_vfieldxx_now = vel_vfieldxx_after
        fvel_vfieldk_now = fvel_vfieldk_after
        dt = dt_after
        time = time_after
        itime = itime_after

        end if
        
        i_poincare = i_poincare + 1

    end if
    
end if
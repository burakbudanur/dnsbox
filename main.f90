program main
    ! modules:
    use numbers
    use openmpi
    use io
    use parameters
    use fftw
    use fieldio
    use rhs
    use vfield
    use timestep
    use stats
    use run
    use solver
    use projector
    use lyap
    use symred
    
    ! initialization:
    call run_init
    if (IC == -3) then
        allocate(shapiro_vfieldk(nx_perproc, ny_half, nz, 3))
    end if
    if (compute_lyap) call lyap_init(vel_vfieldk_now)
    if (slice .and. (i_project > 0 .or. i_save_sliced_fields >0)) then
        call symred_init
        allocate(sliced_vel_vfieldk_now(nx_perproc, ny_half, nz, 3))
    end if
    if (i_project > 0) call projector_init
    if (integrate_invariant) call solver_averaging_init

    write(out, *) "Starting time stepping."

    do
        if (i_finish > 0 .and. itime > i_finish) exit

        if ((i_save_sliced_fields > 0 .and. mod(itime, i_save_sliced_fields) == 0) .or. &
            (i_project > 0 .and. mod(itime, i_project) == 0)) then

            call symred_slice(vel_vfieldk_now, sliced_vel_vfieldk_now)
            call symred_phases_write
        end if

        if (i_save_sliced_fields > 0 .and. itime > i_start .and. mod(itime, i_save_sliced_fields) == 0) then
            write(file_ext, "(i6.6)") itime/i_save_fields
            fname = 'state.'//file_ext//'_sliced'
            call fieldio_write(sliced_vel_vfieldk_now)
            call run_flush_channels
        end if

        if (i_save_fields > 0 .and. itime > i_start .and. mod(itime, i_save_fields) == 0) then
            write(file_ext, "(i6.6)") itime/i_save_fields
            fname = 'state.'//file_ext
            call fieldio_write(vel_vfieldk_now)
            call run_flush_channels
        end if

        if (adaptive_dt .or. (i_print_steps > 0 .and. mod(itime, i_print_steps) == 0)) then
            call timestep_courant(vel_vfieldxx_now)
        end if

        if (adaptive_dt) call timestep_set_dt

        if (i_print_stats > 0 .and. mod(itime, i_print_stats) == 0) then
            
            call stats_compute(vel_vfieldk_now, fvel_vfieldk_now)

            if (integrate_invariant .and. itime > i_start) call solver_averaging_update(vel_vfieldk_now)

            if (IC == -3) then
                ! expected shapiro field
                call vfield_shapiro(time, shapiro_vfieldk)
                ! its norm
                call vfield_norm(shapiro_vfieldk, shapiro_norm, .false.)
                ! difference from the field now
                shapiro_vfieldk(:,:,:,1:3) =  shapiro_vfieldk(:,:,:,1:3) - vel_vfieldk_now(:,:,:,1:3)
                ! delta's norm
                call vfield_norm(shapiro_vfieldk, shapiro_normdelta, .false.)

            end if

        end if

        ! Write stats
        if (i_print_stats > 0 .and. mod(itime, i_print_stats) == 0) then
            call stats_write

            ! Stop if laminarized
            if (my_id==0 .and. terminate_laminar) then
                e_diff = abs(ekin - ekin_lam)/ekin_lam
                input_diff = abs(powerin - powerin_lam)/powerin_lam
                diss_diff = abs(dissip - dissip_lam )/dissip_lam 
                
                if (e_diff < relerr_lam .and. &
                    input_diff < relerr_lam .and. diss_diff < relerr_lam) then
                    kill_switch = .true.
                    write(out, *) "Laminarized, stopping."
                    open(newunit=laminarized_ch,file='LAMINARIZED',position='append')
                    write(laminarized_ch,*) time
                    close(laminarized_ch)
                end if 
            end if

            if (IC == -3 .and. my_id == 0) then

                ! shapiro stats
        
                inquire(file=TRIM(shapiro_file), exist=there, opened=there2)
                if (.not.there) then
                open(newunit=shapiro_ch,file=TRIM(shapiro_file),form='formatted')
                    write(shapiro_ch,"(A2,"//i4_len//","//"3"//sp_len//")") &
                        "# ", "itime", "time", "relerr"
                end if
                if(there.and..not.there2) then
                open(newunit=shapiro_ch,file=TRIM(shapiro_file),position='append')
                end if
                write(shapiro_ch,"(A2,"//i4_f//","//"3"//sp_f//")")&
                    "  ", itime, time, shapiro_normdelta / shapiro_norm
    
                shapiro_written = .true.

            end if

        end if

        if (i_print_steps > 0 .and. mod(itime, i_print_steps) == 0) call timestep_write

        ! project onto POD basis
        if (i_project > 0 .and. mod(itime, i_project) == 0) then
            if (.not. slice) then
                call projector_project(vel_vfieldk_now)
            else
                call projector_project(sliced_vel_vfieldk_now)
            end if
            call projector_write
        end if

        ! spectrum
        if (i_print_spectrum > 0 .and. mod(itime, i_print_spectrum) == 0) call stats_spectra(vel_vfieldk_now)

        ! Compute divergence
        if (log_divergence .and. i_print_stats > 0 &
                    .and. mod(itime, i_print_stats) == 0) call stats_worst_divergence(vel_vfieldk_now)

        ! wall clock limit
        call cpu_time(cput_now)
        if (wall_clock_limit > 0 .and. cput_now - cput_start > wall_clock_limit) then
            write(out, *) "Runtime limit reached, stopping."
            kill_switch = .true.
        end if

        ! Broadcast the kill_switch status
        call MPI_BCAST(kill_switch, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_err)

        ! Flush manually if desired
        if (i_flush > 0 .and. mod(itime, i_flush) == 0) then
            call run_flush_channels
        end if

        if (kill_switch) then
            call run_exit
        end if

        call timestep_precorr(vel_vfieldxx_now, vel_vfieldk_now, fvel_vfieldk_now)

        time = time + dt
        itime = itime + 1

        if (compute_lyap) call lyap_step(vel_vfieldk_now)

    end do

    ! Compute the averages
    if (integrate_invariant) call solver_averaging_finalize

    write(out, *) "Finished time stepping."
    call run_exit
    
end program main


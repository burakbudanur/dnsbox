#include "macros.h"
module run
    use numbers
    use openmpi
    use io
    use parameters
    use fieldio
    use fftw
    use diffops
    use symmops
    use vfield
    use timestep
    use stats
    use projector
    use lyap
    use symred

    real(sp) :: cput_now, cput_start, cput_stop

    complex(dpc), allocatable, dimension(:, :, :, :) :: &
        vel_vfieldk_now, fvel_vfieldk_now, & ! u and F(u) now
        sliced_vel_vfieldk_now, & ! symmetry reduced u
        shapiro_vfieldk

    real(dp), allocatable :: vel_vfieldxx_now(:, :, :, :)

    integer(i4) :: laminarized_ch, shapiro_ch
    logical     :: kill_switch = .false., shapiro_written = .false.

    real(dp) :: e_diff, input_diff, diss_diff, shapiro_norm, shapiro_normdelta
    character(255) :: shapiro_file = 'shapiro.gp'

    contains
    
    subroutine run_init
     
        call openmpi_init
        call io_init
        call parameters_init
        call fftw_init

        allocate(vel_vfieldk_now(nx_perproc, ny_half, nz, 3))
        allocate(fvel_vfieldk_now(nx_perproc, ny_half, nz, 3))
        allocate(vel_vfieldxx_now(nyy, nzz_perproc, nxx, 3))

        ! Initial time
        itime = i_start
        time  = t_start

        if (IC == -3) then
            call vfield_shapiro(time, vel_vfieldk_now)
        elseif (IC == -2) then
            call vfield_laminar(vel_vfieldk_now)
        elseif (IC == -1) then
            ! Random initial condition
            
            call vfield_random(vel_vfieldk_now, .true.)

            write(out, *) "run: Generated a random velocity field."
            write(out, *) "run: Seed for random velocities = ", seed

            write(file_ext, "(i6.6)") 0
            fname = 'state.'//file_ext
            call fieldio_write(vel_vfieldk_now)
        else
            write(file_ext, "(i6.6)") IC
            fname = 'state.'//file_ext
            call fieldio_read(vel_vfieldk_now)
        end if

        if (IC < -1) then
            call vfield_pressure(vel_vfieldk_now)
            call vfield_galinv(vel_vfieldk_now)
            call symmopps_project_all(vel_vfieldk_now,vel_vfieldk_now)
            call vfield_solvediv(vel_vfieldk_now)    
        end if

        call fftw_vk2x(vel_vfieldk_now, vel_vfieldxx_now)
        call rhs_all(vel_vfieldxx_now, vel_vfieldk_now, fvel_vfieldk_now)

        call cpu_time(cput_start)

    end subroutine run_init

!==============================================================================

    subroutine run_exit
        real(dp)    :: run_time_wall, dt_sim
        integer(i4) :: num_time_step_sim

        call run_close_channels
        call cpu_time(cput_stop)

        run_time_wall =  cput_stop - cput_start
        dt_sim = time - t_start
        num_time_step_sim = itime - i_start
        write(out, *) "sec:", run_time_wall
        if (dt_sim > 0) write(out, *) "sec/t:", run_time_wall / dt_sim
        if (num_time_step_sim > 0) write(out, *) "sec/step:", run_time_wall / num_time_step_sim
        write(out, *) "cpusec:", run_time_wall * num_procs
        if (dt_sim > 0) write(out, *) "cpusec/t:", run_time_wall * num_procs / dt_sim
        if (num_time_step_sim > 0) write(out, *) "cpusec/step:", run_time_wall * num_procs / num_time_step_sim
        if (run_time_wall > 0) then
            write(out, *) "fft/sec:", time_fftw / run_time_wall
            write(out, *) "local_transpose/sec:", time_local_transpose / run_time_wall
            write(out, *) "global_transpose/sec:", time_global_transpose / run_time_wall
        end if

        call io_exit    
        call openmpi_exit
        stop

    end subroutine run_exit

!==============================================================================

    subroutine run_flush_channels
        ! Flush text with every state file save
        flush(out)
        if (my_id == 0) then
            if (stats_stat_written) flush(stats_stat_ch)
            if (steps_written) flush(steps_ch)
            if (stats_specs_written) then
                flush(stats_specx_ch)
                flush(stats_specy_ch)
                flush(stats_specz_ch)
            end if
            if (phases_written) flush(phases_ch)
            if (proj_written) flush(proj_ch)
            if (shapiro_written) flush(shapiro_ch)
            if (lyap_written) flush(lyap_out)
        end if
    end subroutine run_flush_channels

!==============================================================================

    subroutine run_close_channels
        if (stats_stat_written) close(stats_stat_ch)
        if (steps_written) close(steps_ch)
        if (stats_specs_written) then
            close(stats_specx_ch)
            close(stats_specy_ch)
            close(stats_specz_ch)
        end if
        if (phases_written) close(phases_ch)
        if (proj_written) close(proj_ch)
        if (shapiro_written) close(shapiro_ch)
        if (lyap_written) close(lyap_out)
    end subroutine run_close_channels

!==============================================================================

end module run
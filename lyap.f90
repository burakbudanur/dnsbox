module lyap

    use numbers
    use openmpi
    use io
    use parameters
    use fieldio
    use fftw
    use vfield
    use timestep

    complex(dpc), allocatable, dimension(:, :, :, :) :: &
        lyap_vfieldk, lyap_fvfieldk
    
    real(dp), allocatable :: lyap_vfieldxx(:, :, :, :)

    integer(i4) :: lyap_out
    real(dp)    :: lambda_max, norm_perturb, lyap_sum, norm_vel, norm_perturb_0

    character(255) :: lyap_out_format
    logical :: lyap_written = .false.

    contains

!==============================================================================
    subroutine lyap_init(vel_vfieldk)
        complex(dpc), intent(in) :: vel_vfieldk(:, :, :, :)
        logical :: lyap_exists, perturb_exists
        integer(i4) :: lyap_itime, read_stat, lyap_in
        real(dp) :: lyap_time, lyap_growth, lyap_lambda_max, time_elapsed
        character(2) :: lyap_temp
    
        lyap_out_format = "(A2,"//i4_f//","//dp_f//","//dp_f//","//dp_f//")"

        allocate(lyap_vfieldk(nx_perproc, ny_half, nz, 3))
        allocate(lyap_fvfieldk, mold=lyap_vfieldk)
        allocate(lyap_vfieldxx(nyy, nzz_perproc, nxx, 3))

        inquire(file = 'lyap.gp', exist = lyap_exists)
        if (lyap_exists .and. my_id == 0) then

            open(newunit = lyap_in, file = 'lyap.gp', status = 'unknown') 

            do 
                read(lyap_in, lyap_out_format, iostat=read_stat)  &
                lyap_temp, lyap_itime, lyap_time, lyap_growth, lyap_lambda_max
                if (read_stat == 0) then
                    lambda_max = lyap_lambda_max
                else if (read_stat < 0) then
                    exit ! end of file reached
                end if
            end do            
            
            time_elapsed  = time - trans_lyap
            lyap_sum = lambda_max * time_elapsed
            
            close(lyap_in)

            open(newunit = lyap_out, file = 'lyap.gp', status = 'old', position= 'append')
            
        else if (my_id == 0) then

            open(newunit = lyap_out, file = 'lyap.gp', status = 'replace')
            write(lyap_out, "(A2,"//i4_len//","//dp_len//","//dp_len//","//dp_len//")") &
                                 "# ",  "itime",     "time",  "growth",  "lambda_max"

            lyap_sum = 0 

        end if
        
        lyap_written = .true.


        write(file_ext, "(i6.6)") itime/i_save_fields         
        fname = 'perturb.'//file_ext
        INQUIRE(file=fname, exist=perturb_exists)

        if (perturb_exists) then 

            call fieldio_read(lyap_vfieldk)
            write(out, *) "lyap: Loading the perturbed field."

            lyap_vfieldk = lyap_vfieldk - vel_vfieldk 

            call vfield_norm(vel_vfieldk, norm_vel, .true.)
            call vfield_norm(lyap_vfieldk, norm_perturb, .true.)

            lyap_vfieldk =  norm_vel * eps_lyap * lyap_vfieldk / norm_perturb 
            norm_perturb_0 = norm_vel * eps_lyap

            lyap_vfieldk = vel_vfieldk + lyap_vfieldk
            call fieldio_write(lyap_vfieldk)
            call fftw_vk2x(lyap_vfieldk, lyap_vfieldxx)

        else

            call vfield_random(lyap_vfieldk, .true.)
            if (k_cutoff > small) then 
                
                call vfield_truncate(k_cutoff, lyap_vfieldk)

            end if


            call vfield_norm(vel_vfieldk, norm_vel, .true.)
            call vfield_norm(lyap_vfieldk, norm_perturb, .true.)

            lyap_vfieldk =  norm_vel * eps_lyap * lyap_vfieldk / norm_perturb 
            norm_perturb_0 = norm_vel * eps_lyap

            lyap_vfieldk = vel_vfieldk + lyap_vfieldk
            call fieldio_write(lyap_vfieldk)
            call fftw_vk2x(lyap_vfieldk, lyap_vfieldxx)

            write(out, *) "lyap: Generated a random perturbation field."
            write(out, *) "lyap: Seed for random velocities = ", seed
            
        end if


    end subroutine lyap_init

!==============================================================================

    subroutine lyap_step(vel_vfieldk)

        complex(dpc), intent(in) :: vel_vfieldk(:, :, :, :)
        real(dp) :: growth, time_elapsed

        ! When this subroutine is called, lyap_vfield is one time step behind 
        ! vel_vfieldk_now, so we start off by time stepping the perturbed field
        call timestep_precorr(lyap_vfieldxx, lyap_vfieldk, lyap_fvfieldk)

        if (mod(itime, i_lyap) == 0) then

            lyap_vfieldk = lyap_vfieldk - vel_vfieldk 
            call vfield_norm(vel_vfieldk, norm_vel, .true.)
            call vfield_norm(lyap_vfieldk, norm_perturb, .true.)
            
            if (time > trans_lyap) then

                growth = norm_perturb / norm_perturb_0
                lyap_sum = lyap_sum + log(growth)
                time_elapsed = time - trans_lyap
                lambda_max = lyap_sum / time_elapsed
                               
                ! Write the estimates
                if (my_id == 0) then
                    write(lyap_out, TRIM(lyap_out_format)) "  ", itime, time, growth, lambda_max
                end if            
            
            end if

            if (k_cutoff > small) then 
            
                call vfield_truncate(k_cutoff, lyap_vfieldk)
    
            end if

            lyap_vfieldk = norm_vel * eps_lyap * lyap_vfieldk / norm_perturb
            lyap_vfieldk = vel_vfieldk + lyap_vfieldk 
            norm_perturb_0 = norm_vel * eps_lyap

        end if

        write(file_ext, "(i6.6)") itime/i_save_fields

        if (mod(itime, i_save_fields) == 0) then
            fname = 'perturb.'//file_ext
            call fieldio_write(lyap_vfieldk)
        end if    

    end subroutine lyap_step

end module lyap

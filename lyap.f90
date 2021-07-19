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
    real(dp)    :: lambda_max, norm_perturbation, &
                   lyap_sum, lyap_time_elapsed 

    character(255) :: lyap_out_format
    logical :: lyap_written = .false.

    contains

!==============================================================================
    subroutine lyap_init(vel_vfieldk)
        complex(dpc), intent(in) :: vel_vfieldk(:, :, :, :)
        logical :: perturb_exists
    
        lyap_out_format = "(A2,"//i4_f//","//dp_f//","//dp_f//","//dp_f//")"

        allocate(lyap_vfieldk(nx_perproc, ny_half, nz, 3))
        allocate(lyap_fvfieldk, mold=lyap_vfieldk)
        allocate(lyap_vfieldxx(nyy, nzz_perproc, nxx, 3))

        ! Create the lyap file, replace it if it exists
        if (my_id == 0) then
            open(newunit=lyap_out,file='lyap.gp',status='replace')
            write(lyap_out, "(A2,"//i4_len//","//dp_len//","//dp_len//","//dp_len//")") &
                                 "# ",  "itime",     "time",  "growth",  "lambda_max"
            lyap_written = .true.
        end if
        
        lyap_sum = 0 
        lyap_time_elapsed = 0

        write(file_ext, "(i6.6)") itime/i_save_fields         
        fname = 'perturb.'//file_ext
        INQUIRE(file=fname, exist=perturb_exists)

        if (perturb_exists) then 

            call fieldio_read(lyap_vfieldk)
            write(out, *) "lyap: Loading the perturbed field."

        else

            call vfield_random(lyap_vfieldk, .true.)
            if (k_cutoff > small) then 
                
                call vfield_truncate(k_cutoff, lyap_vfieldk)

            end if

            call vfield_norm(lyap_vfieldk, norm_perturbation, .true.)
            lyap_vfieldk = eps_lyap * lyap_vfieldk / norm_perturbation 


            call fieldio_write(lyap_vfieldk)
            lyap_vfieldk = vel_vfieldk + lyap_vfieldk
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
            call vfield_norm(lyap_vfieldk, norm_perturbation, .true.)
            
            if (time - t_start > trans_lyap) then

                growth = norm_perturbation / eps_lyap 
                lyap_sum = lyap_sum + log(growth)
                time_elapsed = time - t_start - trans_lyap
                lambda_max = lyap_sum / time_elapsed
                               
                ! Write the estimates
                if (my_id == 0) then
                    write(lyap_out, TRIM(lyap_out_format)) "  ", itime, time, growth, lambda_max
                end if            
            
            end if

            if (k_cutoff > small) then 
            
                call vfield_truncate(k_cutoff, lyap_vfieldk)
    
            end if

            lyap_vfieldk = eps_lyap * lyap_vfieldk / norm_perturbation
            lyap_vfieldk = vel_vfieldk + lyap_vfieldk 

        end if

        write(file_ext, "(i6.6)") itime/i_save_fields

        if (mod(itime, i_save_fields) == 0) then
            fname = 'perturb.'//file_ext
            call fieldio_write(lyap_vfieldk)
        end if    

    end subroutine lyap_step

end module lyap

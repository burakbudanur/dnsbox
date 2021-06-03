module m_step
    
    ! Time stepping routines
    use m_numbers
    use m_openmpi
    use m_io
    use m_parameters
    use m_work
    use m_fields
    use x_fftw
    use m_stats
    use m_state
    
    contains

    subroutine m_step_precorr
        
        ! semi-implicit predictor-corrector method
        integer :: i, j, k, c
        real(dp)  :: ksqr1, ksqr2, invdt, error, tol
        
        ! Compute tolerance
        wrk(:, :, :, 1:3) = fields(:, :, :, 1:3)
        call m_stats_norm
        tol = atol + norm * rtol     
        
        invdt = 1.0d0/dt
        
        ! Prediction: u(n+1)^1 = ((1/dt + (1 - theta)L)u(n) + N(n)) /
        !                        (1/dt - thetaL)
        do k = 1, nx; do j = 1,ny; do i = 1, nz + 1, 2
            
            if (ialias(i, j, k) /= 0) then
                ! Set RHS to 0 for aliasing elements
                
                fields(i, j, k, 4:6) = zero
                fields(i + 1, j, k, 4:6) = zero
                
            else 
                 
                ksqr1 =  akx(j) ** 2 + aky(k) ** 2 + akz(i) ** 2
                ksqr2 =  akx(j) ** 2 + aky(k) ** 2 + akz(i+1) ** 2            
                
                fields(i    , j, k, 4:6) = ((invdt + (1 - theta) * (- nu * ksqr1)) &
                                            * fields(i, j, k, 1:3) &
                                           + state(i    , j, k, 4:6)) & !nonlintrm
                                         / (invdt - theta * (- nu * ksqr1))

                fields(i + 1, j, k, 4:6) = ((invdt + (1 - theta) * (- nu * ksqr2)) &
                                            * fields(i + 1, j, k, 1:3) &
                                           + state(i + 1, j, k, 4:6)) & !nonlintrm
                                         / (invdt - theta * (- nu * ksqr2))
            
            end if
            
        end do; end do; end do
        
        ! fields(:, :, :, 4:6) holds the prediction, u(n+1)^c
        ! fields(:, :, :, 7:9) will store N(n+1)^(c - 1) (N(n+1)^0 = N(n))
        ! fields(:, :, :, 10:12) will store corrections
        
        ! Corrector iterations
        do c = 1, 10
            ! Move nonlinear term to fields(:, :, :, 7:9) 
            fields(:, :, :, 7:9) = state(:, :, :, 4:6)

            ! Compute nonlinear term for u(n+1)^c
            state(:, :, :, 1:3) = fields(:, :, :, 4:6)
            call m_state_rhs_nonlin

            ! Now we have N(n+1)^c in state(:, :, :, 4:6)
            do k = 1, nx; do j = 1,ny; do i = 1, nz + 1, 2
            
                if (ialias(i, j, k) /= 0) then
                    ! Set RHS to 0 for aliasing elements
                    
                    fields(i, j, k, 10:12) = zero
                    fields(i + 1, j, k, 10:12) = zero
                
                else
                                            
                    ksqr1 =  akx(j) ** 2 + aky(k) ** 2 + akz(i) ** 2
                    ksqr2 =  akx(j) ** 2 + aky(k) ** 2 + akz(i+1) ** 2                   
                    
                    fields(i, j, k, 10:12) = theta &
                                           * (state(i, j, k, 4:6) &
                                            - fields(i, j, k, 7:9)) &
                                           / (invdt - theta * (- nu * ksqr1))
                    
                    fields(i + 1, j, k, 10:12) = theta &
                                               * (state(i + 1, j, k, 4:6) &
                                                - fields(i + 1, j, k, 7:9)) &
                                               / (invdt - theta * (- nu * ksqr2))

                    ! Update the prediction
                    fields(i:i+1, j, k, 4:6) = fields(i:i+1, j, k, 4:6) &
                                             + fields(i:i+1, j, k, 10:12)

                end if
                
            end do; end do; end do
            
            ! Compute the corrector norm:
            wrk(:, :, :, 1:3) = fields(:, :, :, 10:12)
            call m_stats_norm            
            error = norm

            if (error < tol) then
                fields(:, :, :, 1:3) = fields(:, :, :, 4:6) ! Accept step
                exit
            end if
            
        end do        
        
        if (error > tol) then 
            write(out, *) 'time step did not converge t = '
            write(out, *) 'terminating at t = ', time
            flush(out)
            stop            
        end if
        
    end subroutine m_step_precorr
    
!==============================================================================    
    
    subroutine m_step_calc_rhs
    
        ! We call this routine from main, seperate from steppers,
        !    to keep courant and normrhs in sync with other stats.
        
        state(:, :, :, 1:3) = fields(:, :, :, 1:3)
        
        comp_courant = .true.
        comp_normrhs = .true.
        
        call m_state_rhs_nonlin
        
        comp_courant = .false.
        comp_normrhs = .false.
        
    end subroutine m_step_calc_rhs
    
end module m_step

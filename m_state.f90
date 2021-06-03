module m_state
    
    use m_numbers
    use m_openmpi
    use m_io
    use m_parameters
    use m_work
    use x_fftw
    use m_stats
    
    real(dp), allocatable :: state(:, :, :, :), tmp(:, :, :)
    
    contains

!==============================================================================
    
    subroutine m_state_init
        
        integer :: ierr
        
        ierr = 0
        
        allocate(state(nz+2, ny, nx, 6), tmp(nz+2, ny, nx), stat=ierr)      
        
        if (ierr/=0) then
            write(out,*) "Cannot allocate state, stopping."
            flush(out)
            stop
        end if
        
        write(out, "('Allocated state arrays')")
        flush(out)

        ! Construct the inner product coefficients
        inprodFac = one / real(nz * ny * nx_all, 8) ** 2
        allocate(inprodCoeffs(nz+2, ny, nx, 3))
        inprodCoeffs = inprodFac
        
        ! akz does not hold the negative frequencies, therefore inner products
        ! need to account for it
        inprodCoeffs(1, :, :, :) = 0.5d0 * inprodFac
        inprodCoeffs(2, :, :, :) = 0.5d0 * inprodFac
        inprodCoeffs(nz + 1, :, :, :) = 0.5d0 * inprodFac
        inprodCoeffs(nz + 2, :, :, :) = 0.5d0 * inprodFac
        
    end subroutine m_state_init

!==============================================================================

    subroutine m_state_rhs_nonlin
        
        !======================================================================
        !
        ! Subroutine computes nonlinear part of the RHS for state(:, :, :, 1:3)
        ! and places the result in state(:, :, :, 4:6) 
        !
        !======================================================================        
        
        integer :: i, j, k, n
        real(dp)  :: t1(0:6), div1, div2, ksqr1, ksqr2, rtmp   
        
        ! Copy state(:, :, :, 1:3) to wrk(:, :, :, 1:3)
        wrk(:, :, :, 1:3) = state(:, :, :, 1:3)
            
        ! Go to physical space:
        do i = 1, 3
            call xFFT3d(-1, i)
            wrk(:, :, :, i) = wrk(:, :, :, i) * normfac
        end do
        
        ! This is a convenient stage to compute the Courant number
        ! Courant number
        ! getting the Courant number (on the master process only)
        
        if (comp_courant) then
            wrk(:,:,:,4) = abs(wrk(:,:,:,1)) + abs(wrk(:,:,:,2)) + abs(wrk(:,:,:,3))
            rtmp = maxval(wrk(1:nz,:,:,4))
            call MPI_REDUCE(rtmp,courant,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,mpi_err)

            courant = courant * dt / dx        
        end if
        
        ! get 6 velocity products:
        do k = 1,nx
            do j = 1,ny
                do i = 1,nz
                    
                    t1(1) = wrk(i, j, k, 1) * wrk(i, j, k, 1) ! u u
                    t1(2) = wrk(i, j, k, 1) * wrk(i, j, k, 2) ! u v
                    t1(3) = wrk(i, j, k, 1) * wrk(i, j, k, 3) ! u w
                    t1(4) = wrk(i, j, k, 2) * wrk(i, j, k, 2) ! v v
                    t1(5) = wrk(i, j, k, 2) * wrk(i, j, k, 3) ! v w
                    t1(6) = wrk(i, j, k, 3) * wrk(i, j, k, 3) ! w w
                    
                    do n = 1, 6
                        wrk(i, j, k, n) = t1(n)
                    end do
                end do
            end do
        end do
        
        ! Back to the Fourier space
        
        do n = 1, 6
            call xFFT3d(1, n)
        end do
        
        ! Now wrk(:, :, :, :, 1:6) holds products of velocity field components
        ! in Fourier space            
        
        do k = 1, nx
            do j = 1,ny
                do i = 1, nz + 1, 2
                    if (ialias(i, j, k) /= 0) then
                        ! Set RHS to 0 for aliasing elements
                        
                        state(i, j, k, 4:6) = zero
                        state(i + 1, j, k, 4:6) = zero
                        
                    else
                        
                        ! ksqr = |k|^2 
                        ksqr1 =  akx(j) ** 2 + aky(k) ** 2 + akz(i) ** 2
                        ksqr2 =  akx(j) ** 2 + aky(k) ** 2 + akz(i+1) ** 2
                                            
                        ! Nonlinear term:
                        ! 
                        
                        ! t1(1) = - Re F{u_i d_i u} 
                        
                        t1(1) =    akx(j) * wrk(i + 1, j, k, 1) & ! ikx u u 
                                    + aky(k)     * wrk(i + 1, j, k, 2) & ! iky v u
                                    + akz(i+1)     * wrk(i + 1, j, k, 3)   ! ikz w u
                        ! t1(2) = - Im F{u_i d_i u}         
                        t1(2) = - (akx(j) * wrk(i, j, k, 1) &         ! ikx u u 
                                    + aky(k) * wrk(i, j, k, 2) &         ! iky v u
                                    + akz(i) * wrk(i, j, k, 3))          ! ikz w u
                        
                        ! t1(3) = - Re F{u_i d_1 v}        
                        t1(3) =    akx(j) * wrk(i + 1, j, k, 2) & ! ikx u v 
                                    + aky(k)     * wrk(i + 1, j, k, 4) & ! iky v v 
                                    + akz(i+1)     * wrk(i + 1, j, k, 5)   ! ikz w v 
                        
                        ! t1(4) = - Im F{u_i d_1 v}        
                        t1(4) = - (akx(j)     * wrk(i, j, k, 2) &     ! ikx u v 
                                    + aky(k)     * wrk(i, j, k, 4) &     ! iky v v 
                                    + akz(i)     * wrk(i, j, k, 5))      ! ikz w v 
                        
                        ! t1(5) = - Re F{u_i d_1 w}
                        t1(5) =    akx(j) * wrk(i + 1, j, k, 3) & ! ikx u w 
                                    + aky(k)     * wrk(i + 1, j, k, 5) & ! iky v w 
                                    + akz(i+1)     * wrk(i + 1, j, k, 6)   ! ikz w w 
                        
                        ! t1(6) = - Im F{u_i d_1 w}
                        t1(6) = - (akx(j)     * wrk(i, j, k, 3) &     ! ikx u w 
                                    + aky(k)     * wrk(i, j, k, 5) &     ! iky v w 
                                    + akz(i)     * wrk(i, j, k, 6))      ! ikz w w                               
                        
                        ! Pressure terms
                        if (ksqr1==0.d0) then
                            state(i    , j, k, 4) = t1(1)
                            state(i    , j, k, 5) = t1(3)
                            state(i    , j, k, 6) = t1(5)
                        else
                            div1 = akx(j)     * t1(1) & 
                                    + aky(k)     * t1(3) &
                                    + akz(i)     * t1(5)

                            state(i    , j, k, 4) = t1(1) - akx(j) * div1 / ksqr1
                            state(i    , j, k, 5) = t1(3) - aky(k) * div1 / ksqr1
                            state(i    , j, k, 6) = t1(5) - akz(i) * div1 / ksqr1 
                        end if

                        if (ksqr2==0.d0) then
                            state(i + 1, j, k, 4) = t1(2)
                            state(i + 1, j, k, 5) = t1(4)
                            state(i + 1, j, k, 6) = t1(6)
                        else
                            div2 = akx(j) * t1(2) &
                                    + aky(k)     * t1(4) &
                                    + akz(i+1)     * t1(6)
                            
                            state(i + 1, j, k, 4) = t1(2) - akx(j) * div2/ksqr2 
                            state(i + 1, j, k, 5) = t1(4) - aky(k) * div2 / ksqr2
                            state(i + 1, j, k, 6) = t1(6) - akz(i+1) * div2 / ksqr2
                        end if

                        ! Adding the Kolmogorov force:                        
                        
                        if ((aky(k) == kF) .and. &
                            (akx(j) == zero) .and. &
                            (akz(i+1) == zero)) then
                                                        
                            state(i + 1, j, k, 4) = state(i + 1, j, k, 4) &
                                                  - 0.5d0 * gamma * real(nz * ny * nx_all, 8)
                            
                        else if ((aky(k) == - kF) .and. &
                                 (akx(j) == zero) .and. &
                                 (akz(i+1) == zero)) then
                            
                            state(i + 1, j, k, 4) = state(i + 1, j, k, 4) & 
                                                  + 0.5d0 * gamma * real(nz * ny * nx_all, 8)
                            
                        end if
                                                    
                    end if
                end do
            end do
        end do
        
        if (comp_normrhs) then 
            ! Add linear part and prepare a complete RHS for norm rhs
            
            ! Construct the right hand side in Fourier space:
            do k = 1, nx
                do j = 1,ny
                    do i = 1, nz + 1, 2
            
                        if (ialias(i, j, k) /= 0) then
                            ! Set RHS to 0 for aliasing elements
                            
                            wrk(i, j, k, 1:3) = zero
                            wrk(i + 1, j, k, 1:3) = zero
                            
                        else
                            ! Linear terms:
                            
                            ! ksqr = |k|^2 
                            ksqr1 =  akx(j) ** 2 + aky(k) ** 2 + akz(i) ** 2
                            ksqr2 =  akx(j) ** 2 + aky(k) ** 2 + akz(i+1) ** 2
                            
                            wrk(i, j, k, 1:3)     = state(i, j, k, 4:6) &
                                                  + (- nu * ksqr1)  &
                                                   * state(i, j, k, 1:3)
                            
                            wrk(i + 1, j, k, 1:3) = state(i + 1, j, k, 4:6) &
                                                  + (- nu * ksqr2)  &
                                                   * state(i + 1, j, k, 1:3)  
                        end if
                    end do
                end do
            end do
                        
            call m_stats_norm
            normrhs = norm     
        end if
        
    end subroutine m_state_rhs_nonlin

!==============================================================================

    subroutine m_state_inprod(res)
        
        ! Compute L^2 inner product <state(:, :, :, 1:3), state(:, :, :, 4:6)>
        
        
        real(dp)              :: res1
        real(dp), intent(out) :: res
        
        res = 0.d0
        res1 =  sum(state(:,:,:,1:3) * inprodCoeffs * state(:,:,:,4:6))
        
        ! summing contributions on each node
        call MPI_ALLREDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, &
                           MPI_COMM_WORLD, mpi_err)
        
    end subroutine m_state_inprod

!==============================================================================

    subroutine m_state_vorticity
        
        ! Compute vorticity for state(:, :, :, 1:3) place the result in 
        ! state(:, :, :, 4:6)
        
         
        
        ! Copy state(:, :, :, 1:3) to wrk(:, :, :, 1:3)
        wrk(:, :, :, 1:3) = state(:, :, :, 1:3)
        
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
        
        ! Now we have vorticity in wrk(:, :, :, 1:3)
        state(:, :, :, 4:6) = wrk(:, :, :, 1:3)
        
    end subroutine m_state_vorticity

!==============================================================================

end module m_state

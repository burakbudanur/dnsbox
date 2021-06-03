include "mt19937_64.f90"
module m_fields
    use m_numbers
    use m_parameters
    use m_io
    use m_work
    use x_fftw
    use mt19937_64

    real(dp), allocatable :: fields(:, :, :, :)
    real(dp) :: energy1, energy
    
    contains 
    
    subroutine m_fields_init
        use m_numbers
        use m_parameters
        
        integer :: ierr
        ierr = 0

        allocate(fields(nz+2, ny, nx, 30), stat=ierr)      
        if (ierr/=0) then
            write(out,*) "Cannot allocate fields, stopping."
            flush(out)
            stop
        end if
        
        write(out, "('Allocated field arrays')")
        flush(out)
                
    end subroutine m_fields_init

!==============================================================================

    subroutine m_fields_reflect_x
        
        integer :: j
        
        wrk(:, :, :, 1:3) = fields(:, :, :, 1:3) 
        
        ! in Fourier space, j indexes x-coordinate
        do j = 1, ny

            if (j == 1) then
                fields(:, j, :, 2:3) = wrk(:, j, :, 2:3)
                fields(:, j, :, 1) = - wrk(:, j, :, 1)

            else
                ! j-th index corresponds to (j - 1)-th positive wave number
                ! corresponding negative frequency is at N - (j - 2)
                fields(:, j, :, 2:3)          = wrk(:, ny - j + 2, :, 2:3)
                fields(:, j, :, 1)          = - wrk(:, ny - j + 2, :, 1)

            end if
            
        end do
            
    end subroutine m_fields_reflect_x

!==============================================================================

    subroutine m_fields_real_reflect_z
        
        integer :: i
        
        wrk(:, :, :, 4:6) = fields(:, :, :, 1:3) 
        
        ! Go to physical space:
        do i = 4, 6
            call xFFT3d(-1, i)
            wrk(:, :, :, i) = wrk(:, :, :, i) * normfac
        end do
        
        wrk(:, :, :, 1:3) = wrk(:, :, :, 4:6)
                
        ! Reflect in z
        do i = 1,nz
            
            if(i == 1) then
            
                ! sigma_z w(z) = - w(- z)
                wrk(i, :, :, 6) = - wrk(i, :, :, 3)

                ! sigma_z [u,v](z) = [u,v](- z)
                wrk(i, :, :, 4:5) = wrk(i, :, :, 1:2)
                        
            else             
            
                ! sigma_z w(z) = - w(- z)
                wrk(i, :, :, 6) = - wrk(nz - (i - 2), :, :, 3)

                ! sigma_z [u,v](z) = [u,v](- z)
                wrk(i, :, :, 4:5) = wrk(nz - (i - 2), :, :, 1:2)
            
            end if
            
        end do
        
        ! Back to Fourier space:
        do i = 4, 6
            call xFFT3d(1, i)
        end do        
        
        fields(:, :, :, 1:3) = wrk(:, :, :, 4:6)
        
    end subroutine m_fields_real_reflect_z

!==============================================================================

    subroutine m_fields_real_reflect_y
        
        integer :: i, j
        
        wrk(:, :, :, 4:6) = fields(:, :, :, 1:3) 
        
        ! Go to physical space:
        do i = 4, 6
            call xFFT3d(-1, i)
            wrk(:, :, :, i) = wrk(:, :, :, i) * normfac
        end do
        
        wrk(:, :, :, 1:3) = wrk(:, :, :, 4:6)
                
        ! Reflect in y
        do j = 1, ny
            
            if(j == 1) then
            
                ! sigma_y v(x) = - v(- y)
                wrk(:, j, :, 5) = - wrk(:, j, :, 2)

                ! sigma_y [u,w](x) = [u,w](- y)
                wrk(:, j, :, 4) = wrk(:, j, :, 1)
                wrk(:, j, :, 6) = wrk(:, j, :, 3)
                        
            else             
            
                ! sigma_y v(y) = - v(- y)
                wrk(:, j, :, 5) = - wrk(:, ny - (j - 2), :, 2)

                ! sigma_y [u,w](x) = [u,w](- y)
                wrk(:, j, :, 4) = wrk(:, ny - (j - 2), :, 1)
                wrk(:, j, :, 6) = wrk(:, ny - (j - 2), :, 3)
            
            end if
            
        end do
        
        ! Back to Fourier space:
        do i = 4, 6
            call xFFT3d(1, i)
        end do        
        
        fields(:, :, :, 1:3) = wrk(:, :, :, 4:6)
        
    end subroutine m_fields_real_reflect_y

!==============================================================================

    subroutine m_fields_shift_z_half
        
        ! Shift fields in z-direction by Lz/2
        
        integer :: i
        
        do i = 1, nz+1, 2

            if (mod((i-1)/2, 2) /= 0) then
            
                ! Real part:
                fields(i, :, :, 1:3) = -fields(i, :, :, 1:3)
                
                ! Imaginary part: 
                fields(i + 1, :, :, 1:3) = -fields(i + 1, :, :, 1:3)
            
            end if
            
        end do
            
    end subroutine m_fields_shift_z_half

!==============================================================================

    subroutine m_fields_shift_y_half
        
        ! Shift fields in y-direction by Ly/2
        
        integer :: j
        
        do j = 1, nx

            if (mod(int(abs(aky(j)), 8), 2) /= 0) then

                fields(:, :, j, 1:3) = -fields(:, :, j, 1:3)
            
            end if
            
        end do
            
    end subroutine m_fields_shift_y_half

!==============================================================================

    subroutine m_fields_shift_x_half

        integer :: k
        
        ! Shift fields in x-direction by Lx/2
        
        do k = 1, ny

            if (mod(int(abs(akx(k)), 8), 2) /= 0) then
            
                ! Real part:
                fields(:, k, :, 1:3) = -fields(:, k, :, 1:3)
            
            end if
            
        end do
        
    end subroutine m_fields_shift_x_half

!==============================================================================

    subroutine m_fields_kolm_Sx
        
        call m_fields_shift_y_half
        call m_fields_reflect_x
        
    end subroutine m_fields_kolm_Sx

!==============================================================================

    subroutine m_fields_real_kolm_RxRy
    
        call m_fields_real_reflect_y
        call m_fields_reflect_x
        
    end subroutine m_fields_real_kolm_RxRy

!==============================================================================

    subroutine m_fields_pressure

        !======================================================================
        !
        ! Apply pressure solver on fields(:, :, :, 1:3) 
        !
        !======================================================================    

        integer :: i, j, k
        real(dp)  :: div1, div2, lapl1, lapl2, p1, p2

        do k = 1,nx; do j = 1,ny;  do i = 1,nz+1,2

            ! getting the divergence (i*k*\hat(u))

            ! getting divergence
            div2 = akx(j) * fields(i, j, k, 1) &
                 + aky(k) * fields(i, j, k, 2) &
                 + akz(i) * fields(i, j, k, 3)
            div1 = - (akx(j) * fields(i + 1, j, k, 1) &
                    + aky(k)   * fields(i + 1, j, k, 2) &
                    + akz(i+1)   * fields(i + 1, j, k, 3))

            ! inverse laplace operator
            lapl1 =  akx(j) ** 2 + aky(k) ** 2 + akz(i  ) ** 2
            lapl2 =  akx(j) ** 2 + aky(k) ** 2 + akz(i+1) ** 2

            ! Taking derivatives of the pressure and subtracting from the 
            ! corresponding velocities
            if (lapl1/=0.d0) then
                p1 = - div1 / lapl1
                fields(i + 1, j, k, 1) = fields(i + 1, j, k, 1) - p1 * akx(j)
                fields(i + 1, j, k, 2) = fields(i + 1, j, k, 2) - p1 * aky(k)
                fields(i + 1, j, k, 3) = fields(i + 1, j, k, 3) - p1 * akz(i+1)
            end if

            if (lapl2/=0.d0) then
                p2 = - div2 / lapl2
                fields(i    , j, k, 1) = fields(i    , j, k, 1) + p2 * akx(j)
                fields(i    , j, k, 2) = fields(i    , j, k, 2) + p2 * aky(k)
                fields(i    , j, k, 3) = fields(i    , j, k, 3) + p2 * akz(i  )
            end if
    
        end do; end do; end do
                
    end subroutine m_fields_pressure

!==============================================================================

    subroutine m_fields_divergence

        !======================================================================
        ! 
        ! Compute divergence of fields(:, :, :, 1:3) and print
        ! 
        !======================================================================    
    
        real(dp) :: dmin, dmax, d1
        real(dp) :: divThreshold = 1.0d-13

        wrk(:,:,:,1:3) = fields(:,:,:,1:3)

        call x_derivative(1,'x',4)
        call x_derivative(2,'y',5)
        call x_derivative(3,'z',6)

        call xFFT3d(-1,4)
        call xFFT3d(-1,5)
        call xFFT3d(-1,6)

        wrk(1:nz, :, :, 1) = normfac*(wrk(1:nz, :, :, 4)&
                                    + wrk(1:nz, :, :, 5)&
                                    + wrk(1:nz, :, :, 6))


        d1 = minval(wrk(1:nz,:,:,1))
        call MPI_REDUCE(d1, dmin, 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, mpi_err)
        d1 = maxval(wrk(1:nz,:,:,1))
        call MPI_REDUCE(d1, dmax, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, mpi_err)

        if (myid==0 .and. (ITIME==ITMIN .or. (abs(dmin) > divThreshold .or. abs(dmax) > divThreshold))) then
            write(out,*) 'time = ', time
            write(out,*) 'divergence:',dmin,dmax
        end if
            
    end subroutine m_fields_divergence  

!==============================================================================

    subroutine m_fields_init_random
        use mt19937_64
        
        integer :: i, j, k, n
        integer :: l8

        real(dp),     allocatable     :: e_spec(:), e_spec1(:)
        
        integer :: n_shell
        
        real(dp)  :: wmag, ratio, fac, fac2
        
        ! Mersenne Twister vars
        integer(i8) :: seed
        integer :: un, istat
            
        ! Random initiation:

        allocate(e_spec(nkmax), e_spec1(nkmax), stat = ierr)
        if (ierr/=0) then
            write(out,*) "cannot allocate init_velocity arrays"
            flush(out)
            stop
        end if
        
        write(out, *) 'generating random velocities'
        flush(out)

        !----------------------------------------------------------------------
        ! Generate the velocities
        !----------------------------------------------------------------------
        
        ! Seed number generation
        !    Adapted from:
        !    https://sourceryinstitute.github.io/OpenCoarrays/proc/init_random_seed.html
        
        if (myid==0) then
                
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
                write(out,*) "OS does not provide random number generator, quitting."
                flush(out)
                stop
            end if
            
        end if

        ! If the parameters.in contains a seed number, use it
        if (SEEDNO /= -1) then
            seed = SEEDNO
        end if

        count = 1
        ! Broadcast the seed
        call MPI_BCAST(seed,count,MPI_INTEGER8,0,MPI_COMM_WORLD,mpi_err)
        
        write(out, *) "RANDOM SEED FOR VELOCITIES = ", seed
        flush(out)
        
        ! This initialises the seed.
        call init_genrand64(seed)

        ! bringing the processors to their own places in the random sequence 6 is
        ! there because we're generating six fields using seed1 because it's int*8
        
        do l8 = 1,myid*(nz+2)*ny*nx*6
            fac = genrand64_real3()
        end do
        
        ! now filling the arrays wrk1...wrk6
        do n = 1, 6
            do k = 1, nx
                do j = 1, ny
                    do i = 1, nz + 2
                        wrk(i, j, k, n) = genrand64_real3()
                    end do
                end do
            end do
        end do
        
        ! making three random arrays with Gaussian PDF out of the six arrays 
        ! we generated
        
        wrk(:, :, :, 1:3) = dsqrt(- two * dlog(wrk(:, :, :, 1:3))) &
                          * dsin(TWO_PI * wrk(:, :, :, 4:6))
                          
        ! Ensure incompressibility
        
        do n = 1,3
            call xFFT3d(1, n)
        end do
        
        
        do k = 1, nx
            do j = 1, ny
                do i = 1, nz + 1, 2
                    
                    n_shell = nint(sqrt(akx(j) ** 2 + aky(k) ** 2 + akz(i) ** 2))
            
                    if (n_shell > 0 .and. ialias(i, j, k) == 0) then
                        
                        ! Construct velocity field from a curl so that its div 0
                        wrk(i, j, k, 4) = - (aky(k) * wrk(i + 1, j, k, 3) &
                                           - akz(i+1) * wrk(i + 1, j, k, 2))
                        wrk(i + 1, j, k, 4) = aky(k) * wrk(i, j, k, 3) & 
                                            - akz(i) * wrk(i, j, k, 2)
                        
                        wrk(i, j, k, 5) = - (akz(i+1)     * wrk(i + 1, j, k, 1) &
                                           - akx(j) * wrk(i + 1, j, k, 3))
                        wrk(i + 1, j, k, 5) = akz(i) * wrk(i, j, k, 1) & 
                                            - akx(j) * wrk(i, j, k, 3)
                        
                        wrk(i, j, k, 6) = - (akx(j) * wrk(i + 1, j, k, 2) &
                                           - aky(k) * wrk(i + 1, j, k, 1))
                        wrk(i + 1, j, k, 6) = akx(j) * wrk(i, j, k, 2) & 
                                            - aky(k) * wrk(i, j, k, 1)
                                        
                    else
                    
                        wrk(i:i + 1, j, k, 4:6) = zero
                    
                    end if
                
                end do
            end do
        end do
        
        fields(:, :, :, 1:3) = wrk(:, :, :, 4:6)
        
        ! Shaping the spectrum
        
        fac = one / real(nz * ny * nx_all) ** 2
        
        e_spec1 = zero
        e_spec = zero
        
        ! assembling the total energy in each shell 
        do k = 1, nx
            do j = 1, ny
                do i = 1, nz + 2
                
                    n_shell = nint(sqrt(real(akx(j) ** 2 &
                                           + aky(k) ** 2 &
                                           + akz(i) ** 2, 4)))  
                    
                    if (n_shell > 0 .and. ialias(i, j, k) == 0) then
                        fac2 = fac * (fields(i, j, k, 1) ** 2 &
                                    + fields(i, j, k, 2) ** 2 &
                                    + fields(i, j, k, 3) ** 2)
                        
                        ! Correct double counting
                        if (i==1 .or. i==2 .or. i==nz+1 .or. i==nz+2) fac2 = fac2 * 0.5d0 
                        
                        e_spec1(n_shell) = e_spec1(n_shell) + fac2
                        
                    end if
                end do
            end do
        end do
        
        ! reducing energy to two arrays on master node
        count = nkmax
        call MPI_REDUCE(e_spec1, e_spec, count, MPI_REAL8, MPI_SUM, 0, &
                        MPI_COMM_WORLD, mpi_err)
        count = nkmax
        call MPI_BCAST(e_spec, count, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
        
        ! Define von Karman spectrum
        
        do k = 1, nkmax
        
            wmag = real(k, 8)
            ratio = wmag / 1.0d0 ! Setting initial peak wavenum arbitrarily to 2 
            
            fac  = two * PI * ratio
            e_spec1(k) = fac ** 4 / (one + fac ** 2) ** 3
                    
        end do
        
        ! normalize it 
        e_spec1 = e_spec1 / sum(e_spec1(1 : nkmax))
        
        ! multiply velocities with appropriate factors:
        
        do k = 1, nx
            do j = 1, ny
                do i = 1, nz + 2
                    
                    n_shell = nint(sqrt(real(akx(j) ** 2 &
                                           + aky(k) ** 2 &
                                           + akz(i) ** 2, 4)))
                    if (n_shell > 0 .and. ialias(i, j, k) == 0 .and. & 
                        e_spec(n_shell) > zero) then
                        
                        fields(i, j, k, 1:3) = fields(i, j, k, 1:3) &
                                             * sqrt(e_spec1(n_shell) &
                                                  / e_spec(n_shell))
                            
                    else
                    
                        fields(i, j, k, 1:3) = zero
                    
                    end if
                    
                
                end do
            end do
        end do
        
        ! dealias
        call m_fields_dealias

        write(out, *) "Generated the velocities"
        flush(out)

        ! deallocate arrays
        deallocate(e_spec, e_spec1, stat=ierr)
        return    
    
    end subroutine m_fields_init_random

!==============================================================================

    subroutine m_fields_dealias
        
        integer :: i
        ! Dealias fields(:, :, :, 1:3)

        do i=1,3
            fields(:, :, :, i) = dealias_matrix(:,:,:) * fields(:, :, :, i)
        end do
    end subroutine m_fields_dealias

!==============================================================================

    subroutine m_fields_energy
        
        
        ! Compute L^2 norm of wrk(:, :, :, 1:3)
        energy  = 0.0d0
        energy1 = sum(inprodCoeffs * fields(:,:,:,1:3)**2)

        ! summing contributions on each node
        call MPI_ALLREDUCE(energy1, energy, 1, MPI_REAL8, MPI_SUM, &
                           MPI_COMM_WORLD, mpi_err)
        
    end subroutine m_fields_energy

end module m_fields

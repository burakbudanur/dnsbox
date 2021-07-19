module parameters
    use numbers
    use openmpi
    use io

    ! DO NOT EDIT ABOVE THIS LINE
    character(7), parameter :: revision = "b829b0d"

    !# Geometry & discretization
    integer(i4) :: &
        nx = 48, ny = 48, nz = 48, & ! number of grid points
        nxx, nyy, nzz, & ! number of super-sampled grid points
        nx_half, ny_half, nz_half, & ! shorthands
        nyy_half_pad1, & ! nyy/2+1 are frequently used in r2c 
                                 ! transforms see "2.3 One-Dimensional DFTs of
                                 ! Real Data" / fftw3 manual
        nx_perproc, nzz_perproc

    integer(i4), parameter :: supsamp_fac = 3 ! supsamp_fac / 2 dealiasing

    real(dp)            :: Lx = 4.0_dp, Lz = 4.0_dp ! grid dimensions
    real(dp), parameter :: Ly = 4.0_dp ! fixed by non-dimensionalization

    !# Physics
    integer(i4) :: forcing = 1 ! none = 0, sine = 1, cosine = 2
    integer(i4), parameter :: qF = 1 ! Forcing wave number
    real(dp) :: Re = 630.0_dp, &  ! Reynolds number
                tilt_angle = 0, & ! (degrees) tilting angle
                smag_const = 0, & ! Smagorinsky constant for LES
                Delta_LES = 0
    logical :: LES = .false. , tilting = .false.

    !# Initiation
    integer(i4) :: IC = -1, &  ! Initial condition 
                               ! (-3 shapiro, -2 laminar, -1 random, 
                               !   0 state.{istart / i_save_fields)
                   i_start = 0 ! Starting time step  
    integer(i8) :: random_seed = -1 ! -1 reads from /dev/urandom
    real(dp) :: random_energy = 0.1_dp, & ! (ekin_laminar) energy of the random IC
                random_smooth = 0.9_dp, & ! smoothness factor for random IC
                t_start = 0 ! Starting time
    
    !# Results -- note i_* convention for variables in time step units
    integer(i4) :: i_print_stats = 20, &   ! Print energy, dissipation etc.
                   i_print_steps = 20, &   ! Print Courant, step error etc.
                   i_save_fields = 2000, & ! Save fields (restart files)
                   i_save_sliced_fields = -1, & ! Save symmetry reduced fields
                   i_print_spectrum = -1, &
                   i_flush = -1, &
                   i_project = -1

    !# Time stepping 
    real(dp)    :: dt = 0.025_dp, &
                   implicitness = 0.5_dp, &
                   steptol = 1.0e-9_dp, &
                   dtmax = 0.1_dp, &
                   courant_target = 0.25_dp
    integer(i4) :: ncorr = 10 
    logical :: adaptive_dt = .true., &
               integrate_invariant = .false. ! if true, use converged time 
                                             ! step of an invariant solution.
    
    !# Termination
    logical  :: terminate_laminar = .true. ! terminate on laminarization
    real(dp) :: relerr_lam = 0.1_dp ! conclude laminarization if stats within
                                    ! relerr_lam of those of laminar
    real(sp) :: wall_clock_limit = -1.0_sp
    integer(i4) :: i_finish = -1 ! time step limit
    
    !# Diagnostics
    logical :: log_divergence  = .false.
    real(dp) :: divergence_th = 1.0e-16_dp, &
                shapiro_alpha, shapiro_beta, shapiro_gamma, shapiro_lambda
    integer(i4) :: shapiro_qalpha = 1, shapiro_qbeta = 1, shapiro_qgamma = 1

    ! symmetries
    logical :: Rxy = .false., Ry = .false., Rz = .false., Sx = .false., &
               Sy = .false., &
               slice = .false.

    ! Projections
    integer(i4) :: num_proj_bases = 10
    logical :: subtract_origin = .false.

    ! Lyapunov exponent
    logical :: compute_lyap = .false.
    integer(i4) :: i_lyap 
    real(dp) :: eps_lyap, trans_lyap, k_cutoff

    integer(i4) :: itime ! simulation timestep count
    real(dp)    :: time, &  ! simulation time
                   dx, dy, dz, &! grid spacings in the expanded physical space
                   norm_fft, &! fft forward+backward needs to be normalized
                   kF ! forcing wave number

    ! laminar values
    real(dp) :: ekin_lam, powerin_lam, dissip_lam

    ! Given 3x3 symmetric matrix M, entries M_{ij} will be used
    ! TODO: perhaps create a data type for such matrices, if that's possible
    integer(i4), parameter, dimension(6) :: isym = (/1, 1, 1, 2, 2, 3/), &
                                            jsym = (/1, 2, 3, 2, 3, 3/)
    integer(i4) :: nsym(3,3)

    namelist /grid/ nx, ny, nz, Lx, Lz
    namelist /physics/ forcing, Re, tilt_angle, smag_const
    namelist /initiation/ IC, random_seed, random_energy, random_smooth, &
                          t_start, i_start
    namelist /output/ i_print_stats, i_print_steps, i_save_fields, &
                      i_print_spectrum, i_flush, i_project, i_save_sliced_fields
    namelist /time_stepping/ dt, implicitness, steptol, ncorr, adaptive_dt, &
                             dtmax, courant_target, integrate_invariant
    namelist /termination/ terminate_laminar, relerr_lam, wall_clock_limit, i_finish
    namelist /debugging/ log_divergence, divergence_th, shapiro_qalpha, shapiro_qbeta, shapiro_qgamma
    namelist /symmetries/ Rxy, Ry, Rz, Sx, Sy, slice
    namelist /lyapunov/ compute_lyap, i_lyap, eps_lyap, trans_lyap, k_cutoff
    namelist /projection/ num_proj_bases, subtract_origin

    contains 

!==============================================================================

    subroutine parameters_init

        integer(i4) :: in, i, j, n

        do n = 1,6
            i = isym(n)
            j = jsym(n)
            nsym(i,j) = n
            nsym(j,i) = n
        end do

        ! Is parameters.in here?
        inquire(file = 'parameters.in', exist=there)
        if(.not.there) then
            write(out,*) 'cannot find parameters.in'
            flush(out)
            error stop
        end if
        
        ! Read parameters from the input file
        open(newunit=in, file='parameters.in', status="old")
        read(in, nml=grid)
        read(in, nml=physics)
        read(in, nml=symmetries)
        read(in, nml=initiation)
        read(in, nml=output)
        read(in, nml=time_stepping)
        read(in, nml=termination)
        read(in, nml=debugging)
        read(in, nml=lyapunov)
        if (i_project > 0) read(in, nml=projection)
        close(in)

        write(out, *) 'dnsbox revision: ', revision
        write(out, *) 'num_procs: ', num_procs
        write(out, *) 'forcing = ', forcing
        write(out, *) 'epsilon = ', epsilon
        write(out, *) 'small = ', small

        write(out, '(79(''=''))')

        ! ---------------------------------------------------------------------

        write(out, "('nx, ny, nz', 3i4)") nx, ny, nz
        write(out, *) 'Lx = ', Lx
        write(out, *) 'Ly = ', Ly
        write(out, *) 'Lz = ', Lz

        if (mod(nx, 2) /= 0 .or. mod(ny, 2) /= 0 .or. mod(nz, 2) /= 0) then
            write(out, *) 'nx, ny and nz must be even'
            flush(out)
            error stop
        end if
        
        nx_perproc = nx / num_procs
        nx_half = nx / 2 
        ny_half = ny / 2 
        nz_half = nz / 2 

        if (nx_perproc * num_procs /= nx) then
            write(out, *) '*** wrong nx:', nx, & 
                          '*** should be divisible by num_procs:', num_procs
            flush(out)
            error stop
        end if

        ! compute the size of the supersampled grid
        nxx = supsamp_fac * nx_half
        nyy = supsamp_fac * ny_half
        nzz = supsamp_fac * nz_half

        nyy_half_pad1 = nyy / 2 + 1

        nzz_perproc  = nzz / num_procs

        ! Pad z in the physical dimension if necessary to make it divisible
        ! by num_procs
        if (nzz_perproc * num_procs /= nzz) then
            nzz_perproc = nzz_perproc + 1
            nzz = nzz_perproc * num_procs

        end if

        if (nx - 1 >= ny_half .and. nx - 1 >= nz - 1) then
            write(out, *) "using x-compact saves."
        elseif (ny_half >= nx - 1 .and. ny_half >= nz - 1) then
            write(out, *) "using y-compact saves."
        else
            write(out, *) "using z-compact saves."
        end if

        ! Not sure what happens to the transpose algorithm if in either
        ! wavenumber or physical space we end up with just one mode in the
        ! distributed dimension
        if (nx_perproc <= 1 .or. nzz_perproc <= 1) then
            write(out, *) '*** num_procs too large:', num_procs, & 
                          '*** nx_perproc:', nx_perproc, &
                          '*** nzz_perproc:', nzz_perproc
            flush(out)
            error stop
        end if

        write(out, '(79(''=''))')

        ! ---------------------------------------------------------------------

        write(out, "('nxx, nyy, nzz', 3i4)") nxx, nyy, nzz

        norm_fft = 1.0_dp / (nxx * nyy * nzz)
        
        dx = Lx / nxx
        dy = Ly / nyy
        dz = Lz / nzz
        
        write(out, '(79(''=''))')
        
        ! ---------------------------------------------------------------------

        write(out, *) 'Re = ', Re

        if (abs(tilt_angle) > small) then
            write(out, *) 'tilt_angle = ', tilt_angle
            tilting = .true.
        end if
        
        if (smag_const > small) then 
            write(out, *) 'smag_const = ', smag_const 
            LES       = .true.
            Delta_LES = Lx / nx
        end if
        
        if (LES .and. &
            ((.not. are_equal(dx, dy)) .or. (.not. are_equal(dx, dz)) &
             .or. (.not. are_equal(dy, dz)))) then
            write(out, *) 'Warning: LES simulation with a nonequal grid'
        end if 

        ! Set the forcing wavenumber
        kF = 2.0_dp * PI * qF / Ly

        ! compute laminar values
        ekin_lam    = 1.0_dp / 4.0_dp
        powerin_lam = PI**2 / (8 * Re)
        dissip_lam  = powerin_lam
        
        write(out, '(79(''=''))')   

        ! ---------------------------------------------------------------------       
        
        write(out, *) 'Rxy = ', Rxy
        write(out, *) 'Ry = ', Ry
        write(out, *) 'Rz = ', Rz
        write(out, *) 'Sx = ', Sx
        write(out, *) 'Sy = ', Sy
        write(out, *) 'slice = ', slice
        write(out, '(79(''=''))')   

        ! ---------------------------------------------------------------------  
        
        if (integrate_invariant) then
            IC = 0
            adaptive_dt = .false.
            write(out, *) 'integrate_invariant = ', integrate_invariant
        end if

        write(out, *) 'IC = ', IC
        write(out, *) 'random_seed = ', random_seed
        write(out, *) 'random_energy = ', random_energy
        write(out, *) 'random_smooth = ', random_smooth
        
        write(out, '(79(''=''))')
                       
        ! ---------------------------------------------------------------------

        if (IC /= 0) then
            i_start = 0
            t_start = 0
        end if

        write(out, *) 't_start = ', t_start
        write(out, *) 'i_start = ', i_start
        write(out, *) 'i_print_stats = ', i_print_stats
        write(out, *) 'i_print_steps = ', i_print_steps
        write(out, *) 'i_save_fields = ', i_save_fields
        write(out, *) 'i_save_sliced_fields = ', i_save_sliced_fields
        write(out, *) 'i_print_spectrum = ', i_print_spectrum
        write(out, *) 'i_flush = ', i_flush
        write(out, *) 'i_project = ', i_project

        write(out, '(79(''=''))')
        
        ! ---------------------------------------------------------------------

        write(out, *) 'dt = ', dt
        write(out, *) 'implicitness = ', implicitness
        write(out, *) 'steptol = ', steptol
        write(out, *) 'ncorr = ', ncorr
        write(out, *) 'adaptive_dt = ', adaptive_dt
        write(out, *) 'dtmax = ', dtmax
        write(out, *) 'courant_target = ', courant_target

        write(out, '(79(''=''))')
        
        ! ---------------------------------------------------------------------

        write(out, *) 'terminate_laminar = ', terminate_laminar
        write(out, *) 'relerr_lam = ', relerr_lam
        write(out, *) 'wall_clock_limit = ', wall_clock_limit
        write(out, *) 'i_finish = ', i_finish

        write(out, '(79(''=''))')
        
       ! ---------------------------------------------------------------------

        write(out, *) 'log_divergence = ', log_divergence
        write(out, *) 'divergence_th = ', divergence_th

        if (IC == -3) then
            shapiro_alpha = (2.0_dp * PI / Lx)  * shapiro_qalpha
            shapiro_beta = (2.0_dp * PI / Ly)  * shapiro_qbeta
            shapiro_gamma = (2.0_dp * PI / Lz)  * shapiro_qgamma
            shapiro_lambda = sqrt(shapiro_alpha**2 + shapiro_beta**2 + &
                                                            shapiro_gamma**2)

            write(out, *) "shapiro_qalpha = ", shapiro_qalpha
            write(out, *) "shapiro_qbeta = ", shapiro_qbeta
            write(out, *) "shapiro_qgamma = ", shapiro_qgamma
            write(out, *) "shapiro_alpha = ", shapiro_alpha
            write(out, *) "shapiro_beta = ", shapiro_beta
            write(out, *) "shapiro_gamma = ", shapiro_gamma
            write(out, *) "shapiro_lambda = ", shapiro_lambda
        end if
        write(out, '(79(''=''))')

        ! ---------------------------------------------------------------------

        if (compute_lyap) then 
            
            write(out, *) 'i_lyap = ', i_lyap
            write(out, *) 'eps_lyap = ', eps_lyap
            write(out, *) 'trans_lyap = ', trans_lyap
            write(out, *) 'k_cutoff = ', k_cutoff

        end if

        write(out, '(79(''=''))')

    end subroutine parameters_init

!==============================================================================
    
end module parameters

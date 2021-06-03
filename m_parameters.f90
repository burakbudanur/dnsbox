module m_parameters
    
    use m_numbers
    use m_openmpi
    use m_io
    
    ! Simulation parameters
    
    character(7) :: revision="82ba8a1"
    ! DO NOT EDIT ABOVE THIS LINE

    integer :: nz,ny,nx,nx_all ! Grid dimensions

    real(dp) :: time
    real(dp) :: Lz
    real(dp) :: dx, dy, dz
    integer  :: nkmax
    real(dp) :: kmax
    
    integer  :: ITIME
    integer  :: ITMIN, ITMAX, IPRINT1, IPRINT2
    real(dp) :: nu, dt
    real(dp), parameter :: gamma = 1.0_dp
    real(dp), parameter :: kF = 1.0_dp
    
    real(dp) :: Courant ! Courant number
    real(dp) :: atol, rtol ! stepping tolerances
    real(dp) :: theta ! implicitness
    
    integer :: IC, Rz, KolmSx
    
    real(dp) :: rnx, rny, rnz
    
    ! logicals etc
    logical :: there, comp_courant, comp_normrhs
    
    ! Stop when LAMINARIZED parameters
    real(dp) :: LAM_THRESHOLD
    integer :: LAM_SWITCH
    
    ! Read the seed number from parameters.in as well
    integer(i8) :: SEEDNO

    ! Time averages
    integer :: integratePeriodic

    ! Inner product coefficients
    real(dp), allocatable :: inprodCoeffs(:, :, :, :)
    real(dp) :: inprodFac

    contains 
    
    subroutine m_parameters_init
        
        ! Is parameters.in here?
        inquire(file = 'parameters.in', exist=there)
        if(.not.there) then
            write(out,*) 'cannot find parameters.in'
            flush(out)
            stop
        end if
        
        ! variable passed will be set to 0 if any of the parameters does not 
        ! make sense
        
        write(out, *) 'dnsbox revision: ', revision
        
        ! Read parameters from the input file
        open(newunit=in, file='parameters.in', form='formatted')
        read(in, *)
        read(in, *)
        read(in, *)
        
        read(in, *, ERR=9000) nz, ny, nx_all
        nx = nx_all / numprocs 

        if (nx * numprocs/=nx_all) then
            write(out, *) '*** wrong nx_all:', nx_all, & 
                          '*** should be divisible by numprocs:', numprocs
            flush(out)
            stop
        end if

        if (ny /= nx_all) then
            write(out, *) 'ny /= nx_all:', ny, nx_all
            flush(out)
            stop
        end if

        if (mod(nx_all, 2) /= 0) then
            write(out, *) 'nx_all is not even:', nx_all
            flush(out)
            stop
        end if

        if (mod(nz, 2) /= 0) then
            write(out, *) 'nz is not even:', nz
            flush(out)
            stop
        end if

        write(out, '(79(''=''))')
        write(out, "('nz, Ny, nx_All', 3i4)") nz, ny, nx_all
        write(out, "('nz, Ny, nx', 3i4)") nz, ny, nx

        read(in, *, ERR=9000, END=9000) Lz
        write(out, *) 'Lz = ', Lz
        
        dx = 2.0d0 * PI / dble(nx_all)        
        dy = 2.0d0 * PI / dble(ny)
        dz = Lz * PI / dble(nz)
        
        read(in, *)
        
        ! ---------------------------------------------------------------------
        
        read(in, *, ERR=9000, END=9000) ITMIN
        write(out, *) 'ITMIN = ', ITMIN
        
        read(in, *, ERR=9000, END=9000) ITMAX
        write(out, *) 'ITMAX = ', ITMAX
        
        read(in, *, ERR=9000, END=9000) IPRINT1
        write(out, *) 'IPRINT1 = ', IPRINT1
        
        read(in, *, ERR=9000, END=9000) IPRINT2
        write(out, *) 'IPRINT2 = ', IPRINT2
        
        read(in, *)
        write(out, '(79(''=''))')
        
        ! ---------------------------------------------------------------------

        read(in, *, ERR=9000, END=9000) nu
        
        write(out, *)  'nu = ', nu
        
        read(in, *)
        write(out, '(79(''=''))')
        

        ! ---------------------------------------------------------------------        
        
        read(in, *, ERR=9000, END=9000) dt
        write(out, *) 'dt = ', dt
        
        read(in, *)
        write(out, '(79(''=''))')
        

        ! ---------------------------------------------------------------------

        read(in, *, ERR=9000, END=9000) IC
        write(out, *) 'IC = ', IC
        
 
        read(in, *)
        write(out, '(79(''=''))')
        
                       
        ! ---------------------------------------------------------------------
                       
        read(in, *, ERR=9000, END=9000) Rz
        write(out, *) 'Rz = ', Rz

        read(in, *, ERR=9000, END=9000) KolmSx
        write(out, *) 'Sx = ', KolmSx

        read(in, *)
        write(out, '(79(''=''))')
        
        
        ! ---------------------------------------------------------------------
                        
        read(in, *, ERR=9000, END=9000) atol
        write(out, *) 'atol = ', atol
                        

        read(in, *, ERR=9000, END=9000) rtol
        write(out, *) 'rtol = ', rtol
                        

        read(in, *, ERR=9000, END=9000) theta
        write(out, *) 'theta = ', theta
                        
        read(in, *)
        write(out, '(79(''=''))')
                           
        
        ! ---------------------------------------------------------------------
        
        read(in, *, ERR=9000, END=9000) LAM_SWITCH
        write(out, *) 'LAM_SWITCH = ', LAM_SWITCH
         
        
        read(in, *, ERR=9000, END=9000) LAM_THRESHOLD
        write(out, *) 'lamthreshold = ', LAM_THRESHOLD
        
        
        read(in, *)
        write(out, '(79(''=''))')
          
   
       ! ---------------------------------------------------------------------
        
        read(in, *, ERR=9000, END=9000) SEEDNO
        write(out, *) 'SEEDNO = ', SEEDNO
        
        read(in, *)
        write(out, '(79(''=''))')
          
        
       ! ---------------------------------------------------------------------

        read(in, *, ERR=9000, END=9000) integratePeriodic
        write(out, *) 'integratePeriodic = ', integratePeriodic

        ! closing the input file
        close(in)
        write(out,'(79(''=''))') 
        

        ! If doing POD runs, only IC=1 should be permitted
        if (integratePeriodic == 1) then
            IC = 1 ! To be sure not to overwrite the IC
        end if
        
        ! defining the rest of the parameters

        rnx = 2.0d0 * real(ny / 2  - 1, 8) / 3.0d0
        write (out, *) 'rnx = ', rnx

        rny = 2.0d0 * real(nx_all / 2  - 1, 8) / 3.0d0
        write (out, *) 'rny = ', rny

        rnz = 2.0d0 * real((nz / 2) * (2.0d0 / Lz)   - 1, 8) / 3.0d0
        write (out, *) 'rnz = ', rnz

        kmax  = real(floor(dsqrt(rnx**2 + rny**2 + rnz**2)),8) ! maximum wave number
        write (out, *) 'kmax = ', kmax
        nkmax = nint(kmax)
        write (out, *) 'nkmax = ', nkmax
        flush(out)
        return

        !----------------------------------------------------------------------
        !  ERROR PROCESSING
        !----------------------------------------------------------------------

        9000 continue
            write(out,*)'An error was encountered while reading parameters.in'
            flush(out)
            stop

    end subroutine m_parameters_init
    
end module m_parameters

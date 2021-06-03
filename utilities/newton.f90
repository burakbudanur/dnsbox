!*************************************************************************
!  Example with one extra parameter T.
!  Can put parameter and constraint = 0 when not required.
!- - - - - - - - - - - - - - -
!  Newton vector:
!    x(1)   = T / scaleT
!    x(2:) = state vector x
!
!  Extra constraint:
!    (F(x)-x). dx/dt = 0 .
!       no update along direction of trajectory.
!
!  Jacobian approximation:
!    dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps
!
!  File guess.in :
!     T, ndts  initial guess for the shooting time, number of timesteps for the shoot
!     (repeat in rows)
!     Set ndts=-1 to find automatically based on dt
!
!*************************************************************************
 module orbitms
   use m_numbers
   use m_parameters
   use m_openmpi
   
   save
                    !* newton.in:
   integer          :: mgmres, nits     ! m for gmres(m), max newton its
   real(dp) :: rel_err          ! relative error
   real(dp) :: del, mndl, mxdl  ! delta for hookstep-trust-region 
   real(dp) :: gtol, epsJ       ! gmres tolerance, eps in eqn above
   integer          :: orbitType        ! 0 for no period updates, 1 for updates
   integer          :: ms               ! number of shoots
                    !* Typical newton.in:
                    !    200  100
                    !    1d-6
                    !   -1d0  1d-14  1d+7
                    !    1d-3  1d-10
                    !    1
                    !    2

   integer          :: nscalars, Nnewt_, Nnewt, nonzeros
   real(dp) :: tol, scaleT    
   real(dp) :: dt__
   real(dp), allocatable :: periods(:)
   integer :: info
   integer, allocatable :: ndtss(:)

   ! Results
   real(dp) :: period
   integer :: ndts

 end module orbitms

!*************************************************************************
include 'NewtonHook.f90'
include 'GMRESm.f90'

program newtonms
    
    ! modules:
    use m_openmpi
    use m_io
    use m_work
    use m_fields
    use x_fftw
    use m_state
    use m_stats
    use m_step
    use m_runs
    
    use orbitms
    use m_newton
    use m_solvers
    
    
    
    external                   :: getrhs, multJ, multJp, saveorbit
    real(dp), external :: dotprod, dotprod_ms
    real(dp)           :: d
    logical                    :: forbidNewton
    
    integer :: un ! dummy i/o
    integer :: ims
    logical :: fexist
    character*1 :: ims_str
    
    ! initialization:
    call m_openmpi_init
    ! call m_openmpi_test
    call m_io_init
    call m_parameters_init
    call m_work_init
    call m_fields_init
    call x_fftw_allocate(1)
    call x_fftw_init
    call m_state_init
    call m_stats_init
    
    inquire(file='NEWTON_NOT_CONVERGED', exist=forbidNewton)
    if (forbidNewton) call m_runs_exit

    inquire(file='NEWTON_DONE', exist=forbidNewton)
    if (forbidNewton) call m_runs_exit

    ! Read newton.in:
    open(newunit=un,status='old',file='newton.in')
    read(un,*) mgmres, nits
    read(un,*) rel_err
    read(un,*) del, mndl, mxdl
    read(un,*) gtol, epsJ
    read(un,*) orbitType
    read(un,*) ms
    close(un)    
    
    write(out, *) 'mgmres, nits = ', mgmres, nits
    write(out, *) 'rel_err = ', rel_err
    write(out, *) 'del, mndl, mxdl = ', del, mndl, mxdl
    write(out, *) 'gtol, epsJ = ', gtol, epsJ
    write(out, *) 'orbitType = ', orbitType
    write(out, *) 'ms = ', ms
    flush(out)

    ! Read relativeSymmetries.in
    call m_solvers_relative_symmetries_read
    
    nscalars = 1 ! look for the period

    ! "Economic" Newton /g 200805
    nonzeros = ((nz+2)*ny*nx - sum(ialias)) * 3
    ! Set the problem size
    if (myid == 0) then
        Nnewt_ = nonzeros + nscalars  
    else
        Nnewt_ = nonzeros
    end if
    
    Nnewt = ms * Nnewt_
    
    ! allocate vectors:
    allocate(new_x(Nnewt))
    allocate(new_fx(Nnewt))
    new_x = zero
    
    allocate(periods(ms))
    allocate(ndtss(ms))
    ! Initial guess of period and number of steps:
    open(newunit=un,status='old',file='guess.in')
    do ims = 1, ms
        read(un,*) periods(ims)
        read(un,*) ndtss(ims)
        if (ndtss(ims) == -1) ndtss(ims) = nint(periods(ims)/dt)
    end do
    close(un)
    
    if (myid == 0) then
        do ims = 0, ms - 1
            new_x(ims*Nnewt_ + 1) = periods(ims + 1)
        end do
    end if
    
    ! Load the state    
    ! Initial time
    itime = ITMIN
    if (ms > 1) then
        write(file_ext, "(i6.6)") itime/IPRINT2
        fname = 'state.'//file_ext//'-0'
        inquire(file=fname, exist=fexist)
        if (fexist) then
            ! read from disk
            do ims = 0, ms - 1
                write(ims_str, "(i1.1)") ims
                fname = 'state.'//file_ext//'-'//ims_str
                call m_work_read
                fields(:, :, :, 1:3) = wrk(:, :, :, 1:3)
                call m_fields_dealias
                call fields2svec(1, ims, new_x)
            end do
        else
            ! construct from time zero
            call m_runs_read
            dt = periods(1) / real(ndtss(1), 8) ! reading state.000000 overwrites dt, restore
            fields(:, :, :, 1:3) = wrk(:, :, :, 1:3)
            call m_fields_dealias
            call fields2svec(1, 0, new_x)
            scaleT = 1d0
            do ims = 1, ms - 1
                call steporbit(ndtss(ims), ims-1, new_x)
                call fields2svec(1, ims, new_x)
            end do
        end if
    else
        call m_runs_read
        fields(:, :, :, 1:3) = wrk(:, :, :, 1:3)
        call m_fields_dealias
        call fields2svec(1, 0, new_x)
    end if
    time = 0d0

    ! Set the scales
    d = dotprod(-1,0,new_x,new_x)
    if(d==0d0) d=1d0
    scaleT = periods(1) / dsqrt(d)
    if (myid == 0) then
        do ims  = 0, ms - 1
            new_x(ims*Nnewt_+1) = new_x(ims*Nnewt_+1) / scaleT
        end do
    end if

    d = dotprod_ms(-1,new_x,new_x)
    if(d==0d0) d=1d0
    tol  = rel_err * dsqrt(d)
    del  = del     * dsqrt(d)
    mndl = mndl    * dsqrt(d)
    mxdl = mxdl    * dsqrt(d)    
    
    info = 1
    call newtonhook(getrhs, multJ, multJp, saveorbit, dotprod_ms, &
                    mgmres, Nnewt, gtol, tol, del, mndl, mxdl, nits, info)

    if (info == 0) then
        call newton_signal_converged
    else
        call newton_signal_not_converged
    end if
     
    call m_runs_exit

    contains    
    
end program newtonms

!==============================================================================
subroutine fields2svec(ni, ims, x)

    ! Move fields(:,:,:,ni:ni+2) to x
    
    use orbitms
    use m_fields
    
    integer,    intent(in)  :: ni
    integer,    intent(in)  :: ims
    real(dp) ,    intent(out)  :: x(Nnewt)
    
    integer :: i, j, k, n, ix
    
    ix = 0
    do k = 1,nx; do j = 1,ny; do i = 1,nz+2; do n = 0, 2

        ! "Economic" Newton /g 200805
        if (ialias(i,j,k) == 0) then
        
            ix = ix + 1
            if (myid == 0) then
                ! master processor also holds scalars being searched for
                x(ims*Nnewt_ + nscalars + ix) = fields(i, j, k, ni + n)
            else
                x(ims*Nnewt_ + ix) = fields(i, j, k, ni + n)
            end if

        end if
        
    end do; end do; end do; end do
    
end subroutine fields2svec

!==============================================================================
subroutine svec2fields(ni, ims, x)
    
    ! Move x to fields(:,:,:,ni:ni+2) 
    
    use orbitms
    use m_fields
    
    integer,    intent(in)  :: ni
    integer,    intent(in)  :: ims
    real(dp) ,    intent(in)  :: x(Nnewt)
    
    integer :: i, j, k, n, ix
    
    ix = 0
    do k = 1,nx; do j = 1,ny; do i = 1,nz+2; do n = 0, 2
        
        ! "Economic" Newton /g 200805
        if (ialias(i,j,k) == 0) then

            ix = ix + 1
            if (myid == 0) then
                ! master processor also holds scalars being searched for
                fields(i, j, k, ni + n) = x(ims*Nnewt_ + nscalars + ix)
            else
                fields(i, j, k, ni + n) = x(ims*Nnewt_ + ix)
            end if

        else
            fields(i, j, k, ni + n) = zero
        end if
        
    end do; end do; end do; end do
    
end subroutine svec2fields

!==============================================================================

subroutine getrhs(n_, x, y)

    !  function to be minimised   
    use orbitms
    use m_newton
    use m_solvers
    
    integer,          intent(in)  :: n_
    real(dp), intent(in)  :: x(Nnewt)
    real(dp), intent(out) :: y(Nnewt)
    real(dp), external :: dotprod
    real(dp) :: y_(Nnewt)
    integer :: ims, ims_

    time=0d0
    do ims = 0, ms -1
        ims_ = modulo(ims+1,ms)
        call steporbit(ndtss(ims+1), ims, x)
        if (ims == ms -1) call m_solvers_relative_symmetries_apply
        call fields2svec(1, ims_, y_)
    end do
    y = y_ - x            ! diff
    
    if (myid == 0) then
        do ims=0, ms-1
            y(ims*Nnewt_+1:ims*Nnewt_+nscalars) = 0d0 ! constraints, rhs=0
        end do
    end if
    
end subroutine getrhs

!==============================================================================

subroutine multJ(n_,x, y)
    
    !  Jacobian of function + lhs of constraints on update
    
    use m_newton
    use orbitms
    use m_io
    
    integer,          intent(in)     :: n_
    real(dp), intent(in)     :: x(Nnewt)
    real(dp), intent(out)    :: y(Nnewt)   
    real(dp), external :: dotprod, dotprod_ms
    real(dp) :: eps, s(Nnewt), d
    integer :: ims
                    ! (F(x0+eps.x)-F(x0))/eps
    
    write(out, *) 'Jacobian call '
    
    
    eps = dsqrt(dotprod_ms(1,x,x))
    
    write(out, *) 'eps = ', eps
    if(eps==0d0)  then
        write(out,*) 'multJ: eps=0 (1)'
        flush(out)
        stop
    end if
    eps = epsJ * dsqrt(dotprod_ms(1,new_x,new_x)) / eps
    write(out, *) 'eps_ = ', eps
    if(eps==0d0)  then 
        write(out,*) 'multJ: eps=0 (2)'
        flush(out)
        stop
    end if
    y = new_x + eps*x
    ! This is for debugging
    d = dotprod_ms(-1, y, y)
    write(out, *) 'Jacobian: normy2 = ', d
    
    call getrhs(n_,y, s)
    ! This is different from openpipeflow
    call getrhs(n_,new_x, new_fx)
    y = (s - new_fx) / eps

    do ims = 0, ms -1
        if (myid == 0) then 
            dt__ = new_x(ims*Nnewt_ + 1) * scaleT / real(ndtss(ims + 1), 8)
        end if
        
        call MPI_BCAST(dt__, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)

        write(out, *) 'MultJ: dt__ = ', dt__

                        ! contstraint, 
                        ! no update in trajectory direction   

        call steporbit(1,ims,new_x)
        call fields2svec(1, ims, s)
        s(ims*Nnewt_+1:(ims+1)*Nnewt_) = (s(ims*Nnewt_+1:(ims+1)*Nnewt_) - new_x(ims*Nnewt_+1:(ims+1)*Nnewt_)) / dt__
        d = dotprod(-1,ims,s,x)
        
        if (myid == 0) then
            y(ims*Nnewt_ + 1) = d
            if (orbitType == 0) y(ims*Nnewt_ + 1) = x(ims*Nnewt_ + 1)
        end if
    end do
    
end subroutine multJ

!==============================================================================

subroutine multJp(n, x)
    use m_numbers
    !  preconditioner for multJ.  Empty - no preconditioner required
    
    
    integer,          intent(in)    :: n
    real(dp), intent(inout) :: x(n)

end subroutine multJp

!==============================================================================

subroutine saveorbit

    use m_newton
    use orbitms
    use m_work
    use m_fields
    use m_runs
    
    real(dp) :: norm_x
    real(dp), external :: dotprod_ms
    integer :: KILL_SWITCH_PERIOD = 0
    real(dp) :: period1
    integer :: un
    integer :: ims
    character*1 :: ims_str

    ndts = sum(ndtss)
    norm_x = dsqrt(dotprod_ms(-1,new_x,new_x))
    
    if (myid == 0) then
        open(newunit=un,status='unknown',access='append',file='newton.dat')
        if(new_nits==0) then
            write(un,"(A2,"//"4"//i4_f//")") "# ", ndts, mgmres, Nnewt, ms
            write(un, "(A2,"//"2"//i4_len//","//"5"//sp_len//")") "# ", "nits", "gits", &
                                        "rel_err", "tol_ratio", "del", "tol", "norm_x"
        end if
        write(un, "(A2,"//"2"//i4_f//","//"5"//sp_f//")") "  ", new_nits, new_gits, &
                            new_tol / norm_x, new_tol / tol, new_del, new_tol, norm_x
        close(un)
    end if
    
    if (myid == 0) then
        period1 = 0d0
        do ims = 0, ms -1
            period1 = period1 + new_x(ims*Nnewt_+1) * scaleT
            ! Kill switch for negative period guesses
            if (new_x(ims*Nnewt_+1) * scaleT < 0) then
                KILL_SWITCH_PERIOD = 1
            end if
        end do

        open(newunit=un,status='unknown',access='append',file='guesses.dat')
        if (ms > 1) then
            write(un,"(A2,"//dp_f//","//i4_f//")") "# ", period1, ndts
            do ims = 0, ms -1
                write(un, "("//i4_f//","//dp_f//","//i4_f//")") new_nits, new_x(ims*Nnewt_+1) * scaleT, ndtss(ims+1)
            end do
        else
            if(new_nits==0) write(un,"(A2,"//i4_f//","//i4_f//","//i4_f//")") "# ", ndts
            write(un, "("//i4_f//","//dp_f//","//dp_f//")") new_nits, period1
        end if
        close(un)

    end if

    ! Broadcast the KILL_SWITCH status
    call MPI_BCAST(KILL_SWITCH_PERIOD, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)

    if (KILL_SWITCH_PERIOD == 1) then
        write(out, *) "Period guess is negative, stopping."
        call newton_signal_not_converged
        call m_runs_exit
    end if
    
    ! Save the state file:
    write(file_ext, "(i6.6)") new_nits
    do ims = 0, ms - 1
        write(ims_str, "(i1.1)") ims
        if (ms > 1) then
            fname = 'newton.'//file_ext//'-'//ims_str
        else
            fname = 'newton.'//file_ext
        end if
        call svec2fields(4, ims, new_x)
        wrk(:, :, :, 1:3) = fields(:, :, :, 4:6)
        call m_work_write ! Write wrk(:, :, :, 1:3)
    end do

end subroutine saveorbit

!==============================================================================

real(dp) function dotprod(n_,ims,a,b)
    ! dot product
    use orbitms
    use m_fields
    use m_state
    
    integer,          intent(in) :: n_        ! Redundant
    integer,          intent(in) :: ims
    real(dp), intent(in) :: a(Nnewt), b(Nnewt)
    real(dp) :: d
    
    call svec2fields(4, ims, a)
    call svec2fields(7, ims, b)
    state(:, :, :, 1:3) = fields(:, :, :, 4:6)
    state(:, :, :, 4:6) = fields(:, :, :, 7:9)
    
    call m_state_inprod(d)  ! L2 inner product between state 1:3 & 4:6
    
    dotprod = d
end function dotprod     

!==============================================================================

real(dp) function dotprod_ms(n_,a,b)
use orbitms

integer,          intent(in) :: n_
real(dp), intent(in) :: a(Nnewt), b(Nnewt)
real(dp), external :: dotprod
real(dp) :: d 
integer :: ims

d = 0d0
do ims = 0, ms-1
   d = d + dotprod(n_,ims,a,b)
end do
dotprod_ms = d

end function dotprod_ms

!==============================================================================

subroutine steporbit(ndts_, ims, x)
    
    ! Time-stepper interface
    use m_numbers
    use m_openmpi
    use m_parameters
    use orbitms
    use m_work
    use m_stats
    use m_step
    use m_fields
    use m_runs
    
    
    integer,          intent(in)  :: ndts_
    integer,          intent(in)  :: ims
    real(dp), intent(in)  :: x(Nnewt)
    
    if (myid == 0 .and. orbitType /= 0 .and. ndts_ /= 1) then 
        dt = x(ims*Nnewt_ + 1) * scaleT / real(ndts_, 8)
    end if
            
    call MPI_BCAST(dt, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)

    write(out, *) 'dt = ', dt
            
    call svec2fields(1, ims, x)

    call m_runs_impose_symmetries
                    
    ! Time-stepping block copied from main.f90:
    do itime = 1, ndts_
        
        call m_step_calc_rhs
        
        call m_step_precorr
        
        call m_fields_pressure
        
        time = time + dt
      
        call m_runs_impose_symmetries

        call m_fields_dealias
        
        ! Flush everything at the first time step
        if (myid == 0 .and. itime == ITMIN + 1) then
            flush(out)
        end if
    end do
    
end subroutine steporbit

subroutine newton_signal_converged
    use m_openmpi
    
    integer :: un
    if (myid == 0) then
        open(newunit=un,file='NEWTON_CONVERGED',position='append')
        write(un,*)
        close(un)
    end if
 end subroutine newton_signal_converged

 subroutine newton_signal_not_converged
    use m_openmpi
    
    integer :: un
    if (myid == 0) then
        open(newunit=un,file='NEWTON_NOT_CONVERGED',position='append')
        write(un,*)
        close(un)
    end if
 end subroutine newton_signal_not_converged

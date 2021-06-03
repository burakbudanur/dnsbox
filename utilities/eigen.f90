!*************************************************************************
! 
!   Arnoldi method for computing eigenvalues/eigenvectors associated with
!   an invariant solution
!
!  File eigen.in :
!     T       period
!     ndts    number of timesteps taken in one period
!     ncgd    number of eigenvalues to be computed
!     
!*************************************************************************
 module orbit
   
   use m_numbers
   use m_parameters
   use m_openmpi
   
   save
                    !* newton.in:
   integer          :: mgmres, nits, ncgd     !m for gmres(m), max newton its
   real(dp) :: rel_err        !relative error
   real(dp) :: del, mndl, mxdl  !delta for hookstep-trust-region 
   real(dp) :: gtol, epsJ    !gmres tolerance, eps in eqn above
                    !* Typical newton.in:
                    !    100  100
                    !    1d-7
                    !   -1d0  1d-7  1d+7
                    !    1d-3  1d-6
                    !    0

   integer          :: nscalars, Nnewt, orbitType, nonzeros
   real(dp) :: tol    
   real(dp) :: dt__, period
   integer :: info, ndts

    
 end module orbit

module m_newton
    use m_numbers
    
    save
    real(dp), allocatable :: new_x(:), new_fx(:)
    real(dp)              :: new_tol,  new_del
    integer                       :: new_nits, new_gits
end module m_newton
 
!*************************************************************************
include 'arnoldi.f'

program eigen
    
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
    
    use orbit
    use m_newton
    use m_solvers
    
    external                   :: getrhs, multJ, multJp, steporbit
    external                   :: multA
    real(dp), external :: dotprod
    real(dp)           :: d
    real(dp), allocatable :: qq(:,:), sv(:), b(:,:), wr(:), wi(:)
    real(dp), allocatable :: h(:,:), y0(:)
    real(dp) :: epsA, rerr
    real(dp) :: norm_x
    integer :: i, j, k, j_, p, ifail
    
    integer :: un ! dummy i/o

    logical :: forbidNewton
    
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
    if (forbidNewton) then
        write(out, *) "Refuse to start from a non-converged Newton search."
        call m_runs_exit
    end if
    
    ! Read relativeSy   mmetries.in
    call m_solvers_relative_symmetries_read

    ! Read newton.in:
    open(newunit=un,status='old',file='newton.in')
    read(un,*) mgmres, nits
    read(un,*) rel_err
    read(un,*) del, mndl, mxdl
    read(un,*) gtol, epsJ
    read(un,*) orbitType
    close(un)

    ! Eigen depends on Newton output
    ! Let's have
    rel_err = rel_err * 10
    epsJ = epsJ * 10
    ! So that we don't try to be as / more precise as / than Newton
    
    write(out, *) 'mgmres, nits = ', mgmres, nits
    write(out, *) 'rel_err = ', rel_err
    write(out, *) 'del, mndl, mxdl = ', del, mndl, mxdl
    write(out, *) 'gtol, epsJ = ', gtol, epsJ
    write(out, *) 'orbitType = ', orbitType
    flush(out)
     
    ! Initial guess of period and number of steps:
    open(newunit=un,status='old',file='guess.in')
    read(un,*) period
    read(un,*) ndts
    if (ndts == -1) ndts = nint(period/dt)
    close(un)

    ! eigen.in:
    open(newunit=un,status='old',file='eigen.in')
    read(un,*) ncgd
    close(un)
    
    nscalars = 1 ! look for the period
    
    ! "Economic" Newton / g 200805
    nonzeros = ((nz+2)*ny*nx - sum(ialias)) * 3
    ! Set the problem size
    if (myid == 0) then
        Nnewt = nonzeros + nscalars  
    else
        Nnewt = nonzeros
    end if
           
    ! Allocate arrays:
    allocate(qq(Nnewt, mgmres), sv(Nnewt), b(Nnewt, ncgd * 2), &
             wr(ncgd * 3), wi(ncgd * 3), h(mgmres, mgmres), y0(Nnewt))
    
    ! allocate vectors:
    allocate(new_x(Nnewt))
    allocate(new_fx(Nnewt))
    new_x = zero

    if (myid == 0) new_x(1) = period
    
    ! Load the state    
    ! Initial time
    itime = ITMIN
    call m_runs_read                 ! Read from file state."(i6.6)", ITIME
    fields(:, :, :, 1:3) = wrk(:, :, :, 1:3)
    call m_fields_dealias
    
    dt = period / real(ndts, 8) ! reading state.000000 overwrites dt, restore
    
    ! Move fields(:,:,:,1:1+2) to new_x
    call fields2svec(1, new_x) 
    ! Scale parameters by norm of state
    d = dotprod(-1,new_x,new_x)
    if (d == 0.0d0) d = 1.0
    tol  = rel_err * dsqrt(d)
    del  = del     * dsqrt(d)
    mndl = mndl    * dsqrt(d)
    mxdl = mxdl    * dsqrt(d)  
    
    norm_x = dsqrt(d)
    
    ! Arnoldi parameters:
    epsA = epsJ
    rerr = min(rel_err, gtol)
    
    call steporbit(ndts, new_x)
    call m_solvers_relative_symmetries_apply
    call fields2svec(1, y0)
    
    ! Initial input:
    call m_fields_init_random ! Mersenne twister

    ! call m_runs_impose_symmetries

    call fields2svec(1, sv)
    
    k = 0
    
    ! main arnoldi loop
    do while(.true.)
        if (myid == 0) sv(1) = 0d0
        
        call arnold(Nnewt, k, mgmres, ncgd, dotprod, rerr, &
                    sv, h, qq, b, wr, wi, ifail)
        if (k > ncgd + 2) then
            write(out, *) 'arnoldi: k = ', k
            do i = 1, ncgd
                write(out, *) 'arnoldi: ', &
                              real(dlog(dsqrt(wr(i)**2 + wi(i) **2)) / period)
            end do
            
        end if
        
        if (ifail == -1) then
            write(out, *) 'arnoldi converged!'
            flush(out)
            call eigen_signal_converged
            exit
        else if(ifail == 0) then
            call multA(sv, y0, epsA, sv)
        else if(ifail == 1) then
            write(out, *) 'arnoldi: reached max its'
            flush(out)
            call eigen_signal_not_converged
            call m_runs_exit
        else if(ifail >= 2) then
            write(out, *) 'arnoldi error'
            deallocate(qq)
            call eigen_signal_not_converged
            call m_runs_exit
        end if        
    end do
                                ! magnitude and angle of Floq mult
    wr(ncgd*2+1:) = dsqrt(wr(:ncgd)**2+wi(:ncgd)**2)
    wi(ncgd*2+1:) = datan2(wi(:ncgd),wr(:ncgd)) !in (-pi,pi]
                                ! convert Floquet to growth rates
    wr(ncgd+1:ncgd*2) = wr(:ncgd)
    wi(ncgd+1:ncgd*2) = wi(:ncgd)
    call logeig(Nnewt,ncgd,b,wr,wi,period,ifail)
    deallocate(qq)
                    ! save eigvals and Floq multipliers
    if(myid == 0) then
    open(newunit=un,status='unknown',file='arnoldi.dat') ! '# its  = ', k
    write(un,"(A2,A4,A2,"//i4_f//")") "# ", "its ", "= ", k
    write(un,"(A2,A5,A2,"//sp_f//")") "# ", "epsA ", "= ", epsA
    write(un,"(A2,A2,A2,"//dp_f//")") "# ", "T ", "= ", period
    write(un, "(A2,"//i2_len//","//"6"//sp_len//")") "# ", "i", "Re(Exponent)", &
    "Im(Exponent)", "Re(Multiplier)", "Im(Multiplier)", "Magnitude", "Angle"
    do i = 1, ncgd
                ! Floq exp, Floq mult, magn & angle
        write(un, "(A2,"//i2_f//","//"6"//sp_f//")") "  ", i, wr(i), wi(i), &
                                wr(ncgd+i), wi(ncgd+i), wr(ncgd*2+i), wi(ncgd*2+i)
    end do
    close(un) 
    end if    

    ! make norm of eigvecs same as that
    ! of converged state and save    
    
    if (norm_x == 0.0d0) norm_x = 1.0d0   ! Workaround for trivial solution
    
    j = 0
    do i = 1, ncgd
        p = 1
        if(wi(i)/=0d0) p = 2
        
        ! Compute norm of the eigenvector
        d = zero
        do j_ = 1, p
            d = d + dotprod(-1, b(:, j+j_), b(:, j+j_))
        end do
        d = dsqrt(d)
        
        ! Rescale eigenvector and save
        do j_ = 1, p
            if(d > 0d0) b(:, j+j_) = b(:, j+j_) * dsqrt(d)
            ! ITIME = (800000 + j + j_) * IPRINT2
            write(file_ext, "(i6.6)") (j + j_)
            fname = 'eigenvector.'//file_ext
            call svec2fields(4, b(:, j+j_))
            wrk(:, :, :, 1:3) = fields(:, :, :, 4:6)
            call m_work_write ! Write wrk(:, :, :, 1:3)
        end do
        j = j + p
    end do
    
    
    ! finalization:
    call m_runs_exit

    contains    
    
end program eigen


!==============================================================================

subroutine fields2svec(ni, x)
    
    ! Move fields(:,:,:,ni:ni+2) to x
    
    use orbit
    use m_fields
    
    integer,    intent(in )  :: ni
    real(dp) ,    intent(out)  :: x(Nnewt)
    
    integer :: i, j, k, n, ix
    
    ix = 0
    do k = 1,nx; do j = 1,ny; do i = 1,nz+2; do n = 0, 2

        ! Experimenting with "economic" Newton / g 200805
        if (ialias(i,j,k) == 0) then
            ix = ix + 1
            if (myid == 0) then
                
                ! master processor also holds scalars being searched for
                x(nscalars + ix) = fields(i, j, k, ni + n)
        
            else
        
                x(ix) = fields(i, j, k, ni + n)
        
            end if
        end if
        
    end do; end do; end do; end do
    
end subroutine fields2svec

!==============================================================================

subroutine svec2fields(ni, x)
    
    ! Move x to fields(:,:,:,ni:ni+2) 
    
    use orbit
    use m_fields
    
    integer,    intent(in)  :: ni
    real(dp) ,    intent(in)  :: x(Nnewt)
    
    integer :: i, j, k, n, ix
    
    ix = 0
    do k = 1,nx; do j = 1,ny; do i = 1,nz+2; do n = 0, 2
        
        ! Experimenting with "economic" Newton / g 200805
        if (ialias(i,j,k) == 0) then

            ix = ix + 1
            if (myid == 0) then
                
                ! master processor also holds scalars being searched for
                fields(i, j, k, ni + n) = x(nscalars + ix)
        
            else
        
                fields(i, j, k, ni + n) = x(ix)
        
            end if

        else
            fields(i, j, k, ni + n) = zero
        end if
        
    end do; end do; end do; end do
    
end subroutine svec2fields

!==============================================================================

subroutine getrhs(n_, x, y)

    !  function to be minimised   
    use orbit
    use m_newton
    use m_solvers
    
    integer,          intent(in)  :: n_
    real(dp), intent(in)  :: x(Nnewt)
    real(dp), intent(out) :: y(Nnewt)
    real(dp), external :: dotprod
    real(dp) :: y_(Nnewt)
    
    time = 0d0
    call steporbit(ndts,x)
    call m_solvers_relative_symmetries_apply
    call fields2svec(1, y_)
    y = (y_ - x)            ! diff
    
    if (myid == 0) then
        y(1) = 0d0                 ! constraints, rhs=0
    end if
    
end subroutine getrhs

!==============================================================================

subroutine multJ(n_,x, y)
    
    !  Jacobian of function + lhs of constraints on update
    
    use m_newton
    use orbit !,    only : Nnewt, epsJ, dt__
    use m_io
    
    integer,          intent(in)     :: n_
    real(dp), intent(inout)  :: x(Nnewt)
    real(dp), intent(out)    :: y(Nnewt)   
    real(dp), external :: dotprod
    real(dp) :: eps, s(Nnewt), d
                    ! (F(x0+eps.x)-F(x0))/eps
    
    write(out, *) 'Jacobian call '
    
    
    eps = dsqrt(dotprod(1,x,x))
    
    write(out, *) 'eps = ', eps
    if(eps==0d0)  then
        write(out,*) 'multJ: eps=0 (1)'
        flush(out)
        stop
    end if
    eps = epsJ * dsqrt(dotprod(1,new_x,new_x)) / eps
    write(out, *) 'eps_ = ', eps
    if(eps==0d0)  then
        write(out,*) 'multJ: eps=0 (2)'
        flush(out)
        stop
    end if
    y = new_x + eps*x
    d = dotprod(-1, y, y)
    write(out, *) 'Jacobian: normy2 = ', d
    
    
    call getrhs(n_,y, s)
    y = (s - new_fx) / eps
                    ! contstraint, 
                    ! no update in trajectory direction
    call steporbit(1,new_x)
    call fields2svec(1, s)
    s = (s - new_x) / dt__
    d = dotprod(-1,s,x)
    
    if (myid == 0) then
        y(1) = d
        if (orbitType == 0) y(1) = x(1)
    end if
    
end subroutine multJ

!==============================================================================

subroutine multJp(n, x)

    !  preconditioner for multJ.  Empty - no preconditioner required
    use m_numbers
    
    integer,          intent(in)    :: n
    real(dp), intent(inout) :: x(n)

end subroutine multJp

!==============================================================================

real(dp) function dotprod(n_,a,b)
    ! dot product
    use orbit
    use m_fields
    use m_state
    
    integer,          intent(in) :: n_        ! Redundant
    real(dp), intent(inout) :: a(Nnewt), b(Nnewt)
    real(dp) :: d
    
    call svec2fields(4, a)
    call svec2fields(7, b)
    state(:, :, :, 1:3) = fields(:, :, :, 4:6)
    state(:, :, :, 4:6) = fields(:, :, :, 7:9)
    
    call m_state_inprod(d)  ! L2 inner product between state 1:3 & 4:6
    
    dotprod = d
end function dotprod     

!==============================================================================

subroutine steporbit(ndts_, x)
    
    ! Time-stepper interface
    use m_numbers
    use m_openmpi
    use m_parameters
    use orbit
    use m_work
    use m_stats
    use m_step
    use m_runs
    
    
    integer,          intent(in)  :: ndts_
    real(dp), intent(in)  :: x(Nnewt)

    ! integer :: un
    
    if (myid == 0 .and. orbitType /= 0 .and. ndts_ /= 1) then 
        dt = x(1) / real(ndts_, 8)
    end if
            
    call MPI_BCAST(dt, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
            
    write(out, *) 'dt = ', dt
          
    call svec2fields(1, x)
    
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

!-------------------------------------------------------------------------
!  linearised time step for arnoldi
!-------------------------------------------------------------------------
subroutine multA(x,y0,epsA, y)
    use m_numbers
    use m_parameters
    use m_newton, only : new_x
    use orbit,  only : ndts, Nnewt
    use m_solvers
    
    real(dp), intent(in)  :: x(Nnewt), y0(Nnewt), epsA
    real(dp), intent(out) :: y(Nnewt)
    real(dp), external :: dotprod
    real(dp) :: eps
   
    y = x   

    eps = dsqrt(dotprod(-1, y, y))
    if(eps==0d0)  then
        write(out,*) 'multA: eps=0 (1)'
        flush(out)
        stop
    end if
    eps = epsA * dsqrt(dotprod(-1, new_x, new_x)) / eps
    if(eps==0d0)  then
        write(out,*) 'multA: eps=0 (2)'
        flush(out)
        stop
    end if
    y = new_x + eps*y
    y(1) = new_x(1)
    call steporbit(ndts,y)
    call m_solvers_relative_symmetries_apply
    call fields2svec(1, y)
    y = (y - y0) / eps
    
end subroutine multA

subroutine eigen_signal_converged
    use m_openmpi
    
    integer :: un
    if (myid == 0) then
        open(newunit=un,file='EIGEN_CONVERGED',position='append')
        write(un,*)
        close(un)
    end if
end subroutine eigen_signal_converged

subroutine eigen_signal_not_converged
    use m_openmpi
    
    integer :: un
    if (myid == 0) then
        open(newunit=un,file='EIGEN_NOT_CONVERGED',position='append')
        write(un,*)
        close(un)
    end if
end subroutine eigen_signal_not_converged
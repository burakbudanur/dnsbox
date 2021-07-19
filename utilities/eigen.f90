#include "../macros.h"
module meigen
    use numbers
    use openmpi
    use io
    use run
    use solver

    contains

    include 'arnoldi.f90' 

!==============================================================================

    subroutine multA(x,y0, y)
        !-------------------------------------------------------------------------
        !  linearised time step for arnoldi
        !-------------------------------------------------------------------------
        use numbers
        use parameters
        
        use solver
        
        real(dp), intent(inout) :: x(nnewt)
        real(dp), intent(in)  :: y0(nnewt)
        real(dp), intent(out) :: y(nnewt)
        real(dp) :: eps
    
        y = x   

        eps = sqrt(solver_dotprod(0,y, y))
        if(abs(eps)<small)  then
            write(out,*) 'multA: eps=0 (1)'
            flush(out)
            stop
        end if
        eps = epsA * sqrt(solver_dotprod(0,new_x, new_x)) / eps
        if(abs(eps)<small)  then
            write(out,*) 'multA: eps=0 (2)'
            flush(out)
            stop
        end if
        y = new_x + eps*y

        call solver_steporbit(ndts,0,y)
        call solver_relative_symmetries_apply(vel_vfieldxx_now, vel_vfieldk_now)
        call solver_vectorize(vel_vfieldk_now, 0, y)
        y = (y - y0) / eps
        
    end subroutine multA

!==============================================================================

    subroutine eigen_signal_converged
        use openmpi
        
        integer(i4) :: un
        if (my_id == 0) then
            open(newunit=un,file='EIGEN_CONVERGED',position='append')
            write(un,*)
            close(un)
        end if
    end subroutine eigen_signal_converged

!==============================================================================

    subroutine eigen_signal_not_converged
        use openmpi
        
        integer(i4) :: un
        if (my_id == 0) then
            open(newunit=un,file='EIGEN_NOT_CONVERGED',position='append')
            write(un,*)
            close(un)
        end if
    end subroutine eigen_signal_not_converged

end module meigen

!==============================================================================

program eigen
    !*************************************************************************
    ! 
    !   Arnoldi method for computing eigenvalues/eigenvectors associated with
    !   an invariant solution
    !     
    !*************************************************************************
    
    ! modules:
    use numbers
    use openmpi
    use io
    use parameters
    use fftw
    use symmops
    use fieldio
    use vfield
    use timestep
    use run
    use solver
    use meigen
    
    real(dp)           :: d
    real(dp), allocatable :: qq(:,:), sv(:), b(:,:), wr(:), wi(:)
    real(dp), allocatable :: h(:,:), y0(:)
    real(dp) :: norm_random
    integer(i4)  :: i, j, k, j_, p, ifail, un
    complex(dpc), allocatable :: random_vfieldk(:, :, :, :)

    logical :: forbidNewton

    call run_init
    IC = 0
    adaptive_dt = .false.
    
    inquire(file='NEWTON_NOT_CONVERGED', exist=forbidNewton)
    if (forbidNewton) then
        write(out, *) "Refuse to start from a non-converged Newton search."
        call run_exit
    end if
    
    call solver_data_read
    ms = 1
    scaleT = 1.0_dp

    ndts = ndtss(1)
    period = periods(1)
    if(ndts==-1) ndts = nint(period/dt)
    if (find_period) dt = period / ndts
    find_period = .false.

    call solver_set_problem_size
               
    ! Allocate arrays:
    allocate(qq(nnewt, aits), sv(nnewt), b(nnewt, ncgd * 2), &
             wr(ncgd * 3), wi(ncgd * 3), h(aits, aits), y0(nnewt))

    allocate(random_vfieldk(nx_perproc, ny_half, nz, 3))
    
    ! allocate vectors:
    allocate(new_x(nnewt))
    allocate(new_fx(nnewt))
    new_x = 0
    
    ! Load the state    
    ! Initial time
    write(file_ext, "(i6.6)") 0
    fname = 'state.'//file_ext
    call fieldio_read(vel_vfieldk_now)
    call solver_vectorize(vel_vfieldk_now, 0, new_x)
    
    call solver_steporbit(ndts, 0, new_x)
    call solver_relative_symmetries_apply(vel_vfieldxx_now, vel_vfieldk_now)
    call solver_vectorize(vel_vfieldk_now, 0, y0)
    
    ! Initial input
    call vfield_random(random_vfieldk, .true.)
    call vfield_norm(random_vfieldk, norm_random, .true.)
    random_vfieldk = random_vfieldk / norm_random
    call solver_vectorize(random_vfieldk, 0, sv)
    
    k = 0
    
    ! main arnoldi loop
    do while(.true.)
        if (my_id == 0) sv(1) = 0
        
        call arnold(nnewt, k, aits, ncgd, solver_dotprod, abs_err, &
                    sv, h, qq, b, wr, wi, ifail)
        if (k > ncgd + 2) then
            write(out, *) 'arnoldi: k = ', k
            do i = 1, ncgd
                write(out, *) 'arnoldi: ', &
                              real(log(sqrt(wr(i)**2 + wi(i) **2)) / period)
            end do
            
        end if
        
        if (ifail == -1) then
            write(out, *) 'arnoldi converged!'
            flush(out)
            call eigen_signal_converged
            exit
        else if(ifail == 0) then
            call multA(sv, y0, sv)
        else if(ifail == 1) then
            write(out, *) 'arnoldi: reached max its'
            flush(out)
            call eigen_signal_not_converged
            call run_exit
        else if(ifail >= 2) then
            write(out, *) 'arnoldi error'
            deallocate(qq)
            call eigen_signal_not_converged
            call run_exit
        end if        
    end do
                                ! magnitude and angle of Floq mult
    wr(ncgd*2+1:) = sqrt(wr(:ncgd)**2+wi(:ncgd)**2)
    wi(ncgd*2+1:) = atan2(wi(:ncgd),wr(:ncgd)) !in (-pi,pi]
                                ! convert Floquet to growth rates
    wr(ncgd+1:ncgd*2) = wr(:ncgd)
    wi(ncgd+1:ncgd*2) = wi(:ncgd)
    call logeig(ncgd,wr,wi,period,ifail)
    deallocate(qq)
                    ! save eigvals and Floq multipliers
    if(my_id == 0) then
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

    ! save eigvecs
    
    j = 0
    do i = 1, ncgd
        p = 1
        if(abs(wi(i))>small) p = 2
        
        ! Compute norm of the eigenvector
        d = 0
        do j_ = 1, p
            d = d + solver_dotprod(0,b(:, j+j_), b(:, j+j_))
        end do
        
        ! Rescale eigenvector to unit norm and save
        do j_ = 1, p
            if(d > 0) b(:, j+j_) = b(:, j+j_) / sqrt(d)
            write(file_ext, "(i6.6)") (j + j_)
            fname = 'eigenvector.'//file_ext
            call solver_tensorize(vel_vfieldk_now, 0, b(:, j+j_))
            call fieldio_write(vel_vfieldk_now)
        end do
        j = j + p
    end do
    
    ! finalization:
    call run_exit 
    
end program eigen

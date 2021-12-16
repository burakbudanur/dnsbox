#include "macros.h"
module solver
    use io
    use parameters
    use stats
    use run

    real(dp), allocatable    :: periods(:), new_x(:), new_fx(:)
    integer(i4), allocatable :: ndtss(:)

    logical :: find_period = .false., &
        find_shift_x = .false., find_shift_z = .false.,&
        rel_Rxy = .false., rel_Ry = .false., rel_Rz = .false., &
        rel_Sx = .false., rel_Sy = .false., rel_TxHalf = .false., &
        rel_TzHalf = .false.

    integer(i4)  :: ndts, &
                    new_nits, new_gits, &
                    ms=1, & ! number of shoots
                    mgmres=50, nits=50, &! m for gmres(m), max newton its
                    ncgd=10, aits=50, & ! arnoldi params
                    i_find_period, nnewt_pershot, nnewt, &
                    averages_ch
    
    real(dp)     :: period, scaleT, shift_x, scale_dex, shift_z, scale_dez, &
                    rel_err=1.0e-11_dp, &
                    del=-1.0_dp, & ! delta for hookstep-trust-region
                    mndl=1.0e-13_dp, &
                    mxdl=1.0e+7_dp, &
                    gtol=1.0e-3_dp, epsJ=1.0e-13_dp, tol, & ! gmres tolerance
                    abs_err=1.0e-10_dp, epsA=1.0e-12_dp, & ! arnoldi tolerances
                    new_tol, new_del, &
                    avg_ekin = 0, avg_powerin = 0, avg_dissip = 0, &
                    std_ekin = 0, std_powerin = 0, std_dissip = 0

    character(255) :: averages_file = 'averages.gp'
    complex(dpc), allocatable :: average(:, :, :, :)

    namelist /invariants/ find_period, find_shift_x, find_shift_z, &
        periods, shift_x, shift_z, ndtss, &
        rel_Rxy, rel_Ry, rel_Rz, rel_Sx, rel_Sy, rel_TxHalf, rel_TzHalf
    namelist /newton_params/ ms, mgmres, nits, rel_err, del, mndl, mxdl, gtol, epsJ
    namelist /arnoldi_params/ ncgd, aits, abs_err, epsA

    contains

    subroutine solver_data_read
        integer(i4) :: in

        write(out, '(79(''=''))')
        write(out, *) "Reading solver parameters."
        write(out, '(79(''=''))')

        open(newunit=in, file='parameters.in', status="old")
        read(in, nml=newton_params)
        read(in, nml=arnoldi_params)
        allocate(periods(ms), ndtss(ms))
        read(in, nml=invariants)
        close(in)

        write(out, *) 'integrate_invariant = ', integrate_invariant

        write(out, *) "Invariant data:"
        write(out, *) 'find_period = ', find_period
        write(out, *) 'find_shift_x = ', find_shift_x
        write(out, *) 'find_shift_z = ', find_shift_z
        write(out, *) "periods", periods
        write(out, *) "shift_x", shift_x
        write(out, *) "shift_z", shift_z
        write(out, *) "ndtss", ndtss
        write(out, *) "rel_Rxy", rel_Rxy
        write(out, *) "rel_Ry", rel_Ry
        write(out, *) "rel_Rz", rel_Rz
        write(out, *) "rel_Sx", rel_Sx
        write(out, *) "rel_Sy", rel_Sy
        write(out, *) "rel_TxHalf", rel_TxHalf
        write(out, *) "rel_TzHalf", rel_TzHalf
        write(out, '(79(''=''))')

        write(out, *) "Newton parameters:"
        write(out, *) 'ms = ', ms
        write(out, *) 'mgmres, nits = ', mgmres, nits
        write(out, *) 'rel_err = ', rel_err
        write(out, *) 'del, mndl, mxdl = ', del, mndl, mxdl
        write(out, *) 'gtol, epsJ = ', gtol, epsJ
        write(out,'(79(''=''))')

        write(out, *) "Arnoldi parameters:"
        write(out, *) 'ncgd = ', ncgd
        write(out, *) 'aits = ', aits
        write(out, *) 'abs_err = ', abs_err
        write(out, *) 'epsA = ', epsA
        write(out,'(79(''=''))')

        flush(out)

    end subroutine solver_data_read

!==============================================================================

    subroutine solver_averaging_init
        call solver_data_read

        if (ms /= 1) then
            write(out, *) "main: exact integration requires ms = 1."
            call run_exit
        end if
        ndts = ndtss(1)
        period = periods(1)
        if(ndts==-1) ndts = nint(period/dt)
        if (find_period) dt = period / ndts
        i_finish = ndts
        i_print_stats = 1
        allocate(average(nx_perproc, ny_half, nz, 3))
        average = 0
    end subroutine solver_averaging_init

!==============================================================================

    subroutine solver_averaging_update(vfieldk)
        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        ! Update the energetics
        avg_ekin = avg_ekin + dt * ekin
        std_ekin = std_ekin + dt * ekin**2
        avg_powerin = avg_powerin + dt * powerin
        std_powerin = std_powerin + dt * powerin**2
        avg_dissip = avg_dissip + dt * dissip
        std_dissip = std_dissip + dt * dissip**2

        average(:,:,:,1:3) = average(:,:,:,1:3) + vfieldk(:,:,:,1:3) * dt
    end subroutine solver_averaging_update

!==============================================================================

    subroutine solver_averaging_finalize
        avg_ekin = avg_ekin / time
        std_ekin = sqrt(std_ekin / time - avg_ekin ** 2)
        avg_powerin = avg_powerin / time
        std_powerin = sqrt(std_powerin / time - avg_powerin ** 2)
        avg_dissip = avg_dissip / time
        std_dissip = sqrt(std_dissip / time - avg_dissip ** 2)

        ! integral of the periodic orbit
        average(:,:,:,1:3) = average(:,:,:,1:3) / time
        write(file_ext, "(i6.6)") 0
        fname = 'average.'//file_ext
        call fieldio_write(average(:,:,:,1:3))

        if (my_id == 0) then
            open(newunit=averages_ch,file=TRIM(averages_file),&
                                            form='formatted',status='replace')
            write(averages_ch,"(A2,"//"6"//dp_len//")") &
                "# ", "avg_ekin", "std_ekin", "avg_powerin", "std_powerin", "avg_dissip", "std_dissip"
            write(averages_ch,"(A2,"//"6"//dp_f//")") &
                            "  ", avg_ekin, std_ekin, avg_powerin, std_powerin, avg_dissip, std_dissip
            close(averages_ch)
        end if

    end subroutine solver_averaging_finalize

!==============================================================================

    subroutine solver_set_problem_size

        i_find_period = 0
        i_find_shift_x = 0
        i_find_shift_z = 0
        if (my_id == 0) then
            if (find_period) i_find_period = 1
            if (find_shift_x) i_find_shift_x = 1
            if (find_shift_z) i_find_shift_z = 1
        end if

        ! all of two dimensions

        if (ix_max /= -1) then
            nnewt_pershot = 2*(nx_perproc-1)*ny_half*(nz-1) * 2 + i_find_period
        else
            nnewt_pershot = 2*nx_perproc*ny_half*(nz-1) * 2 + i_find_period
        end if

        ! only ki=0 of the i'th dimension

        if (nx - 1 >= ny_half .and. nx - 1 >= nz - 1) then
            if (ix_zero /= -1) nnewt_pershot = nnewt_pershot + 2*(nz-1)*ny_half

        elseif (ny_half >= nx - 1 .and. ny_half >= nz - 1) then
            if (ix_max /= -1) then
                nnewt_pershot = nnewt_pershot + 2*(nx_perproc-1)*(nz-1)
            else
                nnewt_pershot = nnewt_pershot + 2*nx_perproc*(nz-1)
            end if

        else
            if (ix_max /= -1) then
                nnewt_pershot = nnewt_pershot + 2*(nx_perproc-1)*ny_half
            else
                nnewt_pershot = nnewt_pershot + 2*nx_perproc*ny_half
            end if

        end if

        nnewt = ms * nnewt_pershot + i_find_shift_s + i_find_shift_z

    end subroutine solver_set_problem_size

!==============================================================================

    subroutine solver_relative_symmetries_apply(vfieldxx, vfieldk)
        use symmops

        real(dp), intent(out) :: vfieldxx(:, :, :, :)
        complex(dpc), intent(inout) :: vfieldk(:, :, :, :)
        real(dp) :: solver_vfieldxx(nyy,nzz_perproc,nxx,3)

        if (rel_TxHalf) call symmops_halfshiftx(vfieldk,vfieldk)
    
        if (rel_TzHalf) call symmops_halfshiftz(vfieldk,vfieldk)

        if (rel_Rz) call symmops_mirrorz(vfieldk,vfieldk)

        if (rel_Rxy .or. rel_Ry .or. rel_Sx .or. rel_Sy) then

            call fftw_vk2x(vfieldk, solver_vfieldxx)

            if (rel_Rxy) call symmops_Rxy(solver_vfieldxx,solver_vfieldxx)

            if (rel_Ry) call symmops_mirrory(solver_vfieldxx,solver_vfieldxx)

            if (rel_Sx) call symmops_Sx(solver_vfieldxx,solver_vfieldxx)
            
            if (rel_Sy) call symmops_Sy(solver_vfieldxx,solver_vfieldxx)

            call fftw_vx2k(solver_vfieldxx, vfieldk)

        end if

        if (rel_TxHalf .or. rel_TzHalf .or. rel_Rz .or. rel_Rxy .or. rel_Ry &
                                                .or. rel_Sx .or. rel_Sy) then

            call fftw_vk2x(vfieldk, vfieldxx)

        end if

    end subroutine solver_relative_symmetries_apply

!==============================================================================

    subroutine solver_vectorize(vfieldk, ims, x)
    
        complex(dpc),   intent(in)  :: vfieldk(:, :, :, :)
        integer(i4),    intent(in)  :: ims
        real(dp),       intent(out) :: x(:)
        
        integer(i4) :: n, ivec_0, i_vec
        _indices

        ivec_0 = ims * nnewt_pershot + i_find_period &
        + ims * (i_find_shift_x + i_find_shift_z) ! shifts in 0th shot

        ivec = 1

        if (nx - 1 >= ny_half .and. nx - 1 >= nz - 1) then
            ! u at kx=0

            if (ix_zero /= -1) then
                do iz = 1, nz; if(iz == iz_max) cycle; do iy = 1, ny_half;
                    
                    x(ivec_0 + ivec) = vfieldk(ix_zero, iy, iz, 1)%re
                    x(ivec_0 + ivec + 1) = vfieldk(ix_zero, iy, iz, 1)%im

                    ivec = ivec + 2
                
                end do; end do;
            end if

            ! v, w
            do n = 2, 3
                _loop_spec_begin
                    
                    x(i_vec_0 + ivec) = vfieldk(ix, iy, iz, n)%re
                    x(i_vec_0 + ivec + 1) = vfieldk(ix, iy, iz, n)%im

                    ivec = ivec + 2
                    
                _loop_spec_end
            end do

        elseif (ny_half >= nx - 1 .and. ny_half >= nz - 1) then
            ! u
            _loop_spec_begin
                
                x(i_vec_0 + ivec) = vfieldk(ix, iy, iz, 1)%re
                x(i_vec_0 + ivec + 1) = vfieldk(ix, iy, iz, 1)%im

                ivec = ivec + 2
                
            _loop_spec_end

            ! v at ky=0
            do iz = 1, nz; if(iz == iz_max) cycle; do ix = 1, nx_perproc; if(ix_max /= -1 .and. ix == ix_max) cycle;
                    
                x(i_vec_0 + ivec) = vfieldk(ix, 1, iz, 2)%re
                x(i_vec_0 + ivec + 1) = vfieldk(ix, 1, iz, 2)%im

                ivec = ivec + 2
            
            end do; end do;

            ! w
            _loop_spec_begin
                
                x(i_vec_0 + ivec) = vfieldk(ix, iy, iz, 3)%re
                x(i_vec_0 + ivec + 1) = vfieldk(ix, iy, iz, 3)%im

                ivec = ivec + 2
                
            _loop_spec_end

        else
            ! u, v
            do n = 1, 2
                _loop_spec_begin
                    
                    x(i_vec_0 + ivec) = vfieldk(ix, iy, iz, n)%re
                    x(i_vec_0 + ivec + 1) = vfieldk(ix, iy, iz, n)%im

                    ivec = ivec + 2
                    
                _loop_spec_end
            end do
        
            ! w at kz=0
            do iy = 1, ny_half; do ix = 1, nx_perproc; if(ix_max /= -1 .and. ix == ix_max) cycle;
                    
                x(i_vec_0 + ivec) = vfieldk(ix, iy, 1, 3)%re
                x(i_vec_0 + ivec + 1) = vfieldk(ix, iy, 1, 3)%im

                ivec = ivec + 2
            
            end do; end do;

        end if
        
    end subroutine solver_vectorize
    
!==============================================================================

    subroutine solver_tensorize(vfieldk, ims, x)
    
        complex(dpc),   intent(out) :: vfieldk(:, :, :, :)
        integer(i4),    intent(in)  :: ims
        real(dp) ,      intent(in)  :: x(:)
        
        integer(i4) :: n, ivec_0, i_vec
        _indices

        ivec_0 = ims * nnewt_pershot + i_find_period &
        + ims * (i_find_shift_x + i_find_shift_z) ! shifts in 0th shot

        ivec = 1

        if (nx - 1 >= ny_half .and. nx - 1 >= nz - 1) then
            ! u at kx=0
            if (ix_zero /= -1) then
                do iz = 1, nz; if(iz == iz_max) cycle; do iy = 1, ny_half;

                    vfieldk(ix_zero, iy, iz, 1)%re = x(i_vec_0 + ivec)
                    vfieldk(ix_zero, iy, iz, 1)%im = x(i_vec_0 + ivec + 1)

                    ivec = ivec + 2
                end do; end do;
            end if

            ! v, w
            do n = 2, 3
                _loop_spec_begin
        
                    vfieldk(ix, iy, iz, n)%re = x(i_vec_0 + ivec)
                    vfieldk(ix, iy, iz, n)%im = x(i_vec_0 + ivec + 1)

                    ivec = ivec + 2
                
                _loop_spec_end
            end do

            ! compute the rest from kx ux + ky uy + kz uz = 0
            _loop_spec_begin
                if (ix_zero/=-1 .and. ix == ix_zero) cycle
                vfieldk(ix,iy,iz,1) = -(vfieldk(ix,iy,iz,2)*ky(iy) &
                                            + vfieldk(ix,iy,iz,3)*kz(iz)) / kx(ix)
            _loop_spec_end

        elseif (ny_half >= nx - 1 .and. ny_half >= nz - 1) then
            ! u
            _loop_spec_begin

                vfieldk(ix, iy, iz, 1)%re = x(i_vec_0 + ivec)
                vfieldk(ix, iy, iz, 1)%im = x(i_vec_0 + ivec + 1)

                ivec = ivec + 2
        
            _loop_spec_end

            ! v at ky=0
            do iz = 1, nz; if(iz == iz_max) cycle; do ix = 1, nx_perproc; if(ix_max /= -1 .and. ix == ix_max) cycle;

                vfieldk(ix, 1, iz, 2)%re = x(i_vec_0 + ivec)
                vfieldk(ix, 1, iz, 2)%im = x(i_vec_0 + ivec + 1)

                ivec = ivec + 2
            end do; end do;

            ! w
            _loop_spec_begin

                vfieldk(ix, iy, iz, 3)%re = x(i_vec_0 + ivec)
                vfieldk(ix, iy, iz, 3)%im = x(i_vec_0 + ivec + 1)

                ivec = ivec + 2
        
            _loop_spec_end

            ! compute the rest from kx ux + ky uy + kz uz = 0
            _loop_spec_begin
                if (iy==1) cycle
                vfieldk(ix,iy,iz,2) = -(vfieldk(ix,iy,iz,1)*kx(ix) &
                                            + vfieldk(ix,iy,iz,3)*kz(iz)) / ky(iy)
            _loop_spec_end

        else
            ! u, v
            do n = 1, 2
                _loop_spec_begin
    
                    vfieldk(ix, iy, iz, n)%re = x(i_vec_0 + ivec)
                    vfieldk(ix, iy, iz, n)%im = x(i_vec_0 + ivec + 1)

                    ivec = ivec + 2
                
                _loop_spec_end
            end do
        
            ! w at kz=0
            do iy = 1, ny_half; do ix = 1, nx_perproc; if(ix_max /= -1 .and. ix == ix_max) cycle;

                vfieldk(ix, iy, 1, 3)%re = x(i_vec_0 + ivec)
                vfieldk(ix, iy, 1, 3)%im = x(i_vec_0 + ivec + 1)

                ivec = ivec + 2
            end do; end do;

            ! compute the rest from kx ux + ky uy + kz uz = 0
            _loop_spec_begin
                if (iz==1) cycle
                vfieldk(ix,iy,iz,3) = -(vfieldk(ix,iy,iz,1)*kx(ix) &
                                            + vfieldk(ix,iy,iz,2)*ky(iy)) / kz(iz)
            _loop_spec_end

        end if

        ! zero maximal modes
        if (ix_max /= -1) vfieldk(ix_max,:,:,1:3) = 0
        vfieldk(:,:,iz_max,1:3) = 0
        
    end subroutine solver_tensorize
    
!==============================================================================

    subroutine solver_steporbit(ndts_, ims, x)
    
        ! Time-stepper interface
        use numbers
        use openmpi
        use io
        use parameters
        use fftw
        use fieldio
        use rhs
        use vfield
        use timestep
        use stats
        use run
         
        integer(i4),          intent(in)  :: ndts_
        integer(i4),          intent(in)  :: ims
        real(dp), intent(in)  :: x(nnewt)
        
        if (my_id == 0 .and. ndts_ /= 1 .and. find_period) then 
            dt = x(ims*nnewt_pershot + 1) * scaleT / ndts_
        end if
                
        call MPI_BCAST(dt, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
    
        write(out, *) 'dt = ', dt
                
        call solver_tensorize(vel_vfieldk_now, ims, x)
        ! get the physical space version too
        call fftw_vk2x(vel_vfieldk_now, vel_vfieldxx_now)
                        
        time = 0
        do itime = 0, ndts_ - 1
    
            call timestep_precorr(vel_vfieldxx_now, vel_vfieldk_now, fvel_vfieldk_now)
            time = time + dt
    
        end do
        
    end subroutine solver_steporbit

!==============================================================================

    real(dp) function solver_dotprod(ims,a,b)
        
        integer(i4),          intent(in) :: ims
        real(dp), intent(in) :: a(nnewt), b(nnewt)
        real(dp) :: d
        complex(dpc) :: vfieldk_a(nx_perproc, ny_half, nz, 3), vfieldk_b(nx_perproc, ny_half, nz, 3)
        
        call solver_tensorize(vfieldk_a, ims, a)
        call solver_tensorize(vfieldk_b, ims, b)
        
        call vfield_inprod(vfieldk_a, vfieldk_b, d, .true.)
        
        solver_dotprod = d
    end function solver_dotprod     

!==============================================================================

    real(dp) function solver_dotprod_ms(a,b)

        real(dp), intent(in) :: a(nnewt), b(nnewt)
        real(dp) :: d 
        integer(i4) :: ims

        d = 0
        do ims = 0, ms-1
        d = d + solver_dotprod(ims,a,b)
        end do
        solver_dotprod_ms = d

    end function solver_dotprod_ms

!==============================================================================

end module solver

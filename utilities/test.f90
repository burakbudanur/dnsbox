#include "../macros.h"
module test_helpers
    use numbers
    use openmpi
    use io
    use parameters
    use fftw
    use vfield
    use symmops
    use test_symmops

    complex(dpc), allocatable, dimension(:, :, :, :) :: &
        routine_in, routine_out, difference, backup, reference

    real(dp), allocatable, dimension(:, :, :, :) :: &
        routine_in_phys, routine_out_phys, difference_phys, backup_phys, &
        reference_phys, in_vfieldxx, out_vfieldxx

    contains

    subroutine test_init
        allocate(routine_in(nx_perproc, ny_half, nz, 3))
        allocate(routine_out, difference, backup, reference, mold=routine_in)

        allocate(routine_in_phys(nyy,nzz_perproc,nxx,3))
        allocate(routine_out_phys, difference_phys, backup_phys, &
                    reference_phys, in_vfieldxx, out_vfieldxx, &
                    mold=routine_in_phys)

    end subroutine test_init

    subroutine test_sym(name, routine)
        character(*), intent(in) :: name
        external                 :: routine

        real(dp) :: delta, refnorm, relerr

        routine_in(:, :, :, :) = reference(:, :, :, :)

        ! check 2-periodicity
        call routine
        routine_in(:, :, :, :) = routine_out(:, :, :, :)
        call routine

        difference(:, :, :, :) = routine_out(:, :, :, :) - reference(:, :, :, :)
        call vfield_norm(difference, delta, .false.)
        call vfield_norm(reference, refnorm, .false.)

        if (my_id == 0) then
            relerr = delta / refnorm
            write(*,"('"//trim(name)//" hit twice, relative error: '"//sp_f//")") relerr
        end if

        routine_in(:, :, :, :) = reference(:, :, :, :)

        ! check idempotency
        call routine
        routine_out(:, :, :, :) = 0.5d0 * routine_in(:, :, :, :) &
                                                + 0.5d0 * routine_out(:, :, :, :)
        ! backup the projection
        backup(:, :, :, :) = routine_out(:, :, :, :)
        ! project again
        routine_in(:, :, :, :) = routine_out(:, :, :, :)
        call routine
        routine_out(:, :, :, :) = 0.5d0 * routine_in(:, :, :, :) &
                                                + 0.5d0 * routine_out(:, :, :, :)
        
        difference(:, :, :, :) = routine_out(:, :, :, :) - backup(:, :, :, :)
        call vfield_norm(difference, delta, .false.)
        call vfield_norm(backup, refnorm, .false.)

        if (my_id == 0) then
            relerr = delta / refnorm
            write(*,"('"//trim(name)//" projected twice, relative error: '"//sp_f//")") &
                                                                                relerr
        end if

    end subroutine test_sym

    subroutine test_sym_phys(name, routine)
        character(*), intent(in) :: name
        external                 :: routine
        real(dp) :: delta, refnorm, relerr

        routine_in_phys(:, :, :, :) = reference_phys(:, :, :, :)

        ! check 2-periodicity
        call routine
        routine_in_phys(:, :, :, :) = routine_out_phys(:, :, :, :)
        call routine

        difference_phys(:, :, :, :) = routine_out_phys(:, :, :, :) - reference_phys(:, :, :, :)
        call test_norm_phys(difference_phys, delta)
        call test_norm_phys(reference_phys, refnorm)

        if (my_id == 0) then
            relerr = delta / refnorm
            write(*,"('"//trim(name)//" hit twice, relative error: '"//sp_f//")") relerr
        end if

        routine_in_phys(:, :, :, :) = reference_phys(:, :, :, :)

        ! check idempotency
        call routine
        routine_out_phys(:, :, :, :) = 0.5d0 * routine_in_phys(:, :, :, :) &
                                                + 0.5d0 * routine_out_phys(:, :, :, :)
        ! backup the projection
        backup_phys(:, :, :, :) = routine_out_phys(:, :, :, :)
        ! project again
        routine_in_phys(:, :, :, :) = routine_out_phys(:, :, :, :)
        call routine
        routine_out_phys(:, :, :, :) = 0.5d0 * routine_in_phys(:, :, :, :) &
                                                + 0.5d0 * routine_out_phys(:, :, :, :)
        
        difference_phys(:, :, :, :) = routine_out_phys(:, :, :, :) - backup_phys(:, :, :, :)
        call test_norm_phys(difference_phys, delta)
        call test_norm_phys(backup_phys, refnorm)

        if (my_id == 0) then
            relerr = delta / refnorm
            write(*,"('"//trim(name)//" projected twice, relative error: '"//sp_f//")") &
                                                                                relerr
        end if

    end subroutine test_sym_phys

    subroutine test_sym_compare(name, routine, routine_test)
        character(*), intent(in) :: name
        external                 :: routine

        real(dp) :: delta, refnorm, relerr
    
        routine_in(:, :, :, :) = reference(:, :, :, :)

        ! run original
        call routine
        backup(:, :, :, :) = routine_out(:, :, :, :)

        ! run alternative
        routine_in(:, :, :, :) = reference(:, :, :, :)
        call routine_test

        difference(:, :, :, :) = routine_out(:, :, :, :) - backup(:, :, :, :)
        call vfield_norm(difference, delta, .false.)
        call vfield_norm(backup, refnorm, .false.)

        if (my_id == 0) then
            relerr = delta / refnorm
            write(*,"('"//trim(name)//" compared with test, relative error: '"//sp_f//")") &
                                                                                relerr
        end if

    end subroutine test_sym_compare

    subroutine test_projection(name, routine)
        character(*), intent(in) :: name
        external                 :: routine
        real(dp) :: delta, refnorm, relerr

        routine_in(:, :, :, :) = reference(:, :, :, :)

        call routine
        ! backup the projection
        backup(:, :, :, :) = routine_out(:, :, :, :)
        ! apply constraints again
        routine_in(:, :, :, :) = routine_out(:, :, :, :)
        call routine

        difference(:, :, :, :) = routine_out(:, :, :, :) - backup(:, :, :, :)
        call vfield_norm(difference, delta, .false.)
        call vfield_norm(backup, refnorm, .false.)

        if (my_id == 0) then
            relerr = delta / refnorm
            write(*,"('"//trim(name)//" projected twice, relative error: '"//sp_f//")") &
                                                                                relerr
        end if

    end subroutine test_projection

    subroutine test_norm_phys(vfieldxx,res)
        real(dp), intent(in) :: vfieldxx(:, :, :, :)
        real(dp), intent(out) :: res
        real(dp)      :: res1
        real(dp)      :: dummy
        _indicess
        integer(i4) :: n

        res1 = 0
        do n=1,3
            _loop_phys_begin
                dummy = vfieldxx(iyy,izz,ixx,n)**2
                res1 = res1 + dummy
            _loop_phys_end
        end do

        ! normalize
        res1 = (res1 * norm_fft) / 2.0_dp

        call MPI_REDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, 0, &
        MPI_COMM_WORLD, mpi_err)

        res = sqrt(res)
    end subroutine test_norm_phys

    ! wrappers for routines to be tested with arguments to call them with the 
    ! above testers

    ! all of below could be simplified using the "abstract" of fortran

    ! all take as input routine_in, and put the result to routine_out
    ! physical-space versions do the same for _phys arrays

    ! symmmetries
    ! the order these are presented is:
    !   multi-core version (used runtime), alternative single-core version,
    !   physical-space version, spectral space version

    subroutine mirrorx
        call fftw_vk2x(routine_in,in_vfieldxx)
        call symmops_mirrorx(in_vfieldxx, out_vfieldxx)
        call fftw_vx2k(out_vfieldxx,routine_out)
    end subroutine mirrorx

    subroutine mirrorx_test
        call test_symmops_mirrorx(routine_in, routine_out)
    end subroutine mirrorx_test

    subroutine mirrorx_phys
        call symmops_mirrorx(routine_in_phys, routine_out_phys)
    end subroutine mirrorx_phys

    subroutine mirrorx_spec
        call test_symmops_mirrorx(routine_in, routine_out)
    end subroutine mirrorx_spec

    !

    subroutine mirrory
        call fftw_vk2x(routine_in,in_vfieldxx)
        call symmops_mirrory(in_vfieldxx, out_vfieldxx)
        call fftw_vx2k(out_vfieldxx,routine_out)
    end subroutine mirrory

    subroutine mirrory_test
        call test_symmops_mirrory(routine_in, routine_out)
    end subroutine mirrory_test

    subroutine mirrory_phys
        call symmops_mirrory(routine_in_phys, routine_out_phys)
    end subroutine mirrory_phys

    subroutine mirrory_spec
        call test_symmops_mirrory(routine_in, routine_out)
    end subroutine mirrory_spec

    !

    subroutine mirrorz
        call symmops_mirrorz(routine_in, routine_out)
    end subroutine mirrorz

    subroutine mirrorz_test
        call fftw_vk2x(routine_in,in_vfieldxx)
        call test_symmops_mirrorz(in_vfieldxx, out_vfieldxx)
        call fftw_vx2k(out_vfieldxx,routine_out)
    end subroutine mirrorz_test

    subroutine mirrorz_phys
        call test_symmops_mirrorz(routine_in_phys, routine_out_phys)
    end subroutine mirrorz_phys

    subroutine mirrorz_spec
        call symmops_mirrorz(routine_in, routine_out)
    end subroutine mirrorz_spec

    !

    subroutine shiftx
        call symmops_shiftx(Lx/2, routine_in, routine_out)
    end subroutine shiftx

    !

    subroutine shiftz
        call symmops_shiftz(Lz/2, routine_in, routine_out)
    end subroutine shiftz

    !

    subroutine halfshiftx
        call symmops_halfshiftx(routine_in, routine_out)
    end subroutine halfshiftx

    subroutine halfshiftx_test
        call fftw_vk2x(routine_in,in_vfieldxx)
        call test_symmops_halfshiftx(in_vfieldxx, out_vfieldxx)
        call fftw_vx2k(out_vfieldxx,routine_out)
    end subroutine halfshiftx_test

    subroutine halfshiftx_phys
        call test_symmops_halfshiftx(routine_in_phys, routine_out_phys)
    end subroutine halfshiftx_phys

    subroutine halfshiftx_spec
        call symmops_halfshiftx(routine_in, routine_out)
    end subroutine halfshiftx_spec

    !

    subroutine halfshifty
        call fftw_vk2x(routine_in,in_vfieldxx)
        call symmops_halfshifty(in_vfieldxx, out_vfieldxx)
        call fftw_vx2k(out_vfieldxx,routine_out)
    end subroutine halfshifty

    subroutine halfshifty_test
        call test_symmops_halfshifty(routine_in, routine_out)
    end subroutine halfshifty_test

    subroutine halfshifty_phys
        call symmops_halfshifty(routine_in_phys, routine_out_phys)
    end subroutine halfshifty_phys

    subroutine halfshifty_spec
        call test_symmops_halfshifty(routine_in, routine_out)
    end subroutine halfshifty_spec

    !

    subroutine halfshiftz
        call symmops_halfshiftz(routine_in, routine_out)
    end subroutine halfshiftz

    subroutine halfshiftz_test
        call fftw_vk2x(routine_in,in_vfieldxx)
        call test_symmops_halfshiftz(in_vfieldxx, out_vfieldxx)
        call fftw_vx2k(out_vfieldxx,routine_out)
    end subroutine halfshiftz_test

    subroutine halfshiftz_phys
        call test_symmops_halfshiftz(routine_in_phys, routine_out_phys)
    end subroutine halfshiftz_phys

    subroutine halfshiftz_spec
        call symmops_halfshiftz(routine_in, routine_out)
    end subroutine halfshiftz_spec

    ! other routines to test

    subroutine pressure
        routine_out(:, :, :, :) = routine_in(:, :, :, :)
        call vfield_pressure(routine_out)
    end subroutine pressure

    subroutine galinv
        routine_out(:, :, :, :) = routine_in(:, :, :, :)
        call vfield_galinv(routine_out)
    end subroutine galinv

    subroutine solvediv
        routine_out(:, :, :, :) = routine_in(:, :, :, :)
        call vfield_solvediv(routine_out)
    end subroutine solvediv

end module test_helpers

program test

    use numbers
    use openmpi
    use io
    use parameters
    use fieldio
    use fftw
    use diffops
    use symmops
    use vfield
    use test_helpers
    use solver

    real(dp) :: delta, refnorm, relerr, my_div, div
    complex(dpc), allocatable :: div_sfieldk(:, :, :)
    real(dp), allocatable :: div_sfieldxx(:, :, :)
    real(dp), allocatable :: vector(:)
    
    call openmpi_init
    call io_init
    call parameters_init
    call fftw_init
    call test_init
    
    ! get random initial condition with no symmetry
    ! might be good to fix the seed
    call vfield_random(reference, .false.)
    call fftw_vk2x(reference,reference_phys)
    
    ! check state file writing ################################################
    routine_in(:, :, :, :) = reference(:, :, :, :)

    ! save the random initial condition
    write(file_ext, "(i6.6)") 0
    fname = 'dummy.'//file_ext
    call fieldio_write(routine_in)
    ! read again
    call fieldio_read(routine_in)
    difference(:, :, :, :) = routine_in(:, :, :, :) - reference(:, :, :, :)
    call vfield_norm(difference, delta, .false.)
    call vfield_norm(reference, refnorm, .false.)
    if (my_id == 0) then
        relerr = delta / refnorm
        write(*,"('"//"i/o"//", relative error: '"//sp_f//")") relerr
    end if

    ! #########################################################################

    ! check fft ###############################################################
    routine_in(:, :, :, :) = reference(:, :, :, :)

    call fftw_vk2x(routine_in,in_vfieldxx)
    call fftw_vx2k(in_vfieldxx,routine_out)
    difference(:, :, :, :) = routine_out(:, :, :, :) - routine_in(:, :, :, :)
    call vfield_norm(difference, delta, .false.)
    call vfield_norm(routine_in, refnorm, .false.)
    if (my_id == 0) then
        relerr = delta / refnorm
        write(*,"('"//"fft"//", relative error: '"//sp_f//")") relerr
    end if

    ! ######################################################################### 

    ! check divergence

    allocate(div_sfieldk(nx_perproc,ny_half,nz))
    allocate(div_sfieldxx(nyy,nzz_perproc,nxx))

    call diffops_div(reference, div_sfieldk)
    call fftw_sk2x(div_sfieldk, div_sfieldxx)

    my_div = maxval(abs(div_sfieldxx(:, :, :)))
    call MPI_REDUCE(my_div, div, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, mpi_err)

    if (my_id == 0) then
        write(*,"('"//"worst divergence: '"//sp_f//")") relerr
    end if

    deallocate(div_sfieldk, div_sfieldxx)

    ! #########################################################################

    ! check vectorize / tensorize #############################################
    routine_in(:, :, :, :) = reference(:, :, :, :)

    ms = 1
    find_period = .false.
    call solver_set_problem_size
    allocate(vector(nnewt))

    call solver_vectorize(routine_in, 0, vector)
    call solver_tensorize(routine_out, 0, vector)
    difference(:, :, :, :) = routine_out(:, :, :, :) - routine_in(:, :, :, :)
    call vfield_norm(difference, delta, .false.)
    call vfield_norm(routine_in, refnorm, .false.)
    if (my_id == 0) then
        relerr = delta / refnorm
        write(*,"('"//"vectorize/tensorize"//", relative error: '"//sp_f//")") relerr
    end if

    deallocate(vector)

    ! #########################################################################

    ! check each symmetry action
    ! multi-core, runtime versions first
    call test_sym_phys("mirrorx", mirrorx_phys)
    call test_sym_phys("mirrory", mirrory_phys)
    call test_sym("mirrorz", mirrorz_spec)
    call test_sym("halfshiftx", halfshiftx_spec)
    if (mod(nyy, 2) == 0) call test_sym_phys("halfshifty", halfshifty_phys)
    call test_sym("halfshiftz", halfshiftz_spec)
    call test_sym("shiftx(Lx/2)", shiftx)
    call test_sym("shiftz(Lz/2)", shiftz)

    ! check each symmetry action against the alternative implementation
    if (num_procs == 1) then

        write(*,*) "checking single-core versions of the symmetries"

        call test_sym("mirrorx", mirrorx_spec)
        call test_sym("mirrory", mirrory_spec)
        call test_sym_phys("mirrorz", mirrorz_phys)
        if (mod(nxx, 2) == 0) call test_sym_phys("halfshiftx", halfshiftx_phys)
        call test_sym("halfshifty", halfshifty_spec)
        if (mod(nzz, 2) == 0) call test_sym_phys("halfshiftz", halfshiftz_phys)

        write(*,*) "comparing multi-core implementations to single-core"

        call test_sym_compare("mirrorx", mirrorx, mirrorx_test)
        call test_sym_compare("mirrory", mirrory, mirrory_test)
        call test_sym_compare("mirrorz", mirrorz, mirrorz_test)
        if (mod(nxx, 2) == 0) call test_sym_compare("halfshiftx", halfshiftx, halfshiftx_test)
        if (mod(nyy, 2) == 0) call test_sym_compare("halfshifty", halfshifty, halfshifty_test)
        if (mod(nzz, 2) == 0) call test_sym_compare("halfshiftz", halfshiftz, halfshiftz_test)
    
    else
        write(*,*) "run the tests single core to check symmetries against"
        write(*,*) "their single core implementations."
    end if

    ! other routines to test

    call test_projection("pressure", pressure)
    call test_projection("galinv", galinv)
    call test_projection("solvediv", solvediv)

    call io_exit    
    call openmpi_exit

end program test


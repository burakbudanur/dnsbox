#include "macros.h"
program main
    ! modules:
    use numbers
    use openmpi
    use io
    use parameters
    use fftw
    use diffops 
    use vfield
    use rhs
    use timestep
    use run
    
    integer(i4) :: d, iter
    real(dp) :: local_error, relative_error, norm_velfield, norm_test_fieldk, max_error
    complex(dpc), allocatable :: velfieldk(:, :, :, :), div_velfield(:, :, :), test_fieldk(:, :, :, :)
    real(dp), allocatable :: test_fieldxx(:, :, :, :)
    
    _indices

    ! initialization:
    call openmpi_init
    call io_init
    call parameters_init
    call fftw_init

    write(out, '(79(''=''))')
    write(out, *) "# Tests begin."
    
    allocate(div_velfield(nx_perproc, ny_half, nz), stat=ierr)
    allocate(velfieldk(nx_perproc, ny_half, nz, 3), stat=ierr)
    allocate(test_fieldk(nx_perproc, ny_half, nz, 3), stat=ierr)
    allocate(test_fieldxx(nyy, nzz_perproc, nxx, 3), stat=ierr)

    if (ierr/=0) then
        write(out,*) "Cannot allocate test variables, stopping."
        flush(out)
        error stop
    end if

    relative_error = 0
    test_fieldxx = 0
    test_fieldk = 0

    ! initialize from a random field with no symmetries [gokhan 21-06-24]
    call vfield_random(velfieldk, .false.)

    call vfield_norm(velfieldk, norm_velfield, .true.)
    
    ! write(out, *) "Random state generator relative_error = ", relative_error / norm_velfield
    
    write(out, '(79(''=''))')
    write(out, *) "Testing forward-backward Fourier transforms ..."

    call fftw_vk2x(velfieldk, test_fieldxx)
    call fftw_vx2k(test_fieldxx, test_fieldk)
    call vfield_norm(test_fieldk, norm_test_fieldk, .true.)
    call vfield_norm(velfieldk - test_fieldk, relative_error, .true.)
    relative_error = relative_error / norm_velfield

    write(out, *) "norm_velfield = ", norm_velfield
    write(out, *) "norm_test_fieldk = ", norm_test_fieldk 
    write(out, *) "Relative error after k2x and x2k transforms = ",  relative_error

    if(relative_error < epsilon) then
        write(out, *) "forward-backward Fourier transform -- Pass"
    else
        write(out, *) "forward-backward Fourier transform -- Fail"
        flush(out)
        error stop
    end if

    write(out, '(79(''=''))')
    write(out, *) "Computing the velocity field's divergence ..."

    call diffops_div(velfieldk, div_velfield)
    max_error = maxval(abs(div_velfield))
    write(out, *) "max_error = ", max_error

    if(max_error < epsilon) then
        write(out, *) "divergence test -- Pass"
    else
        write(out, *) "divergence test -- Fail"
        flush(out)
        error stop
    end if

    
    write(out, '(79(''=''))')
    write(out, *) "# Tests over. -- Passed all."
    write(out, '(79(''=''))')
    call run_exit
    
end program main
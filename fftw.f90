#include "macros.h"

module fftw_include
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
end module fftw_include

!==============================================================================

module fftw
    use numbers
    use openmpi
    use io
    use parameters
    use fftw_include

    ! field arrays for fft
    complex(dpc), pointer, dimension(:, :, :) :: &
        fftw_sfieldk, fftw_sfieldkx, fftw_sfieldkx_tmp1, fftw_sfieldkx_tmp1_t, &
        fftw_sfieldkx_tmp2
    real(dp), pointer :: fftw_sfieldxx(:, :, :)

    type(C_PTR) :: p_fftw_sfieldk, p_fftw_sfieldxx, p_fftw_sfieldkx_tmp1, &
        p_fftw_sfieldkx_tmp2

    real(dp), allocatable, dimension(:) :: kx, ky, kz ! wave numbers
    integer(i4), allocatable, dimension(:) :: qx, qy, qz ! integer wave numbers

    integer(i4) :: &
        ix_zero = -1, iy_force = -1, & ! indices of ky=kF and kx=0
        ix_first_p = -1, ix_first_n = -1, &
        ix_max = -1, iz_max = -1, & ! indices of max(kx) and max(kz)
        my_id_ix_max, my_id_ix_max1=-1 ! id of process with max(kx)

    ! fftw plans
    type(C_PTR) :: p_plan_x_forward,  p_plan_y_forward,  p_plan_z_forward, &
                   p_plan_x_backward, p_plan_y_backward, p_plan_z_backward

    ! benchmarking
    real(sp)    :: time_fftw = 0, time_local_transpose = 0, &
        time_global_transpose = 0

    contains

!==============================================================================
    
    subroutine fftw_allocate

        allocate(kx(nx_perproc), ky(ny_half), kz(nz), &
                 qx(nx_perproc), qy(ny_half), qz(nz))

        ! allocating fftw_sfieldk arrays with fftw_alloc_complex and pointing 
        ! to them with c_f_pointer to make fftw happy with regards to memory 
        ! alignment
        ! this also allows several arrays to share the same memory

        p_fftw_sfieldk = fftw_alloc_complex(int(nx_perproc * ny_half * nz, i8))
        ! fftw_sfieldk has its own memory
        call c_f_pointer(p_fftw_sfieldk, fftw_sfieldk, [nx_perproc, ny_half, nz]) 

        p_fftw_sfieldxx = fftw_alloc_complex(& 
            int(nyy_half_pad1  * nzz_perproc * nxx, i8))
        ! fftw_sfieldxx and fftw_sfieldkx share memory
        call c_f_pointer(&
            p_fftw_sfieldxx, fftw_sfieldxx, [2 * nyy_half_pad1 , nzz_perproc, nxx])
        call c_f_pointer(&
            p_fftw_sfieldxx, fftw_sfieldkx, [nyy_half_pad1 , nzz_perproc, nxx])
            
        p_fftw_sfieldkx_tmp1 = fftw_alloc_complex(int(ny_half * nx_perproc * nzz, i8))
        ! fftw_sfieldkx_tmp1 and fftw_sfieldkx_tmp1_t share memory
        call c_f_pointer(p_fftw_sfieldkx_tmp1, fftw_sfieldkx_tmp1,   &
                         [ny_half, nx_perproc,  nzz])
        call c_f_pointer(p_fftw_sfieldkx_tmp1, fftw_sfieldkx_tmp1_t, &
                         [ny_half, nzz_perproc, nx])

        p_fftw_sfieldkx_tmp2 = fftw_alloc_complex(int(ny_half * nzz_perproc * nxx, i8))
        ! fftw_sfieldkx_tmp2 has its own memory
        call c_f_pointer(p_fftw_sfieldkx_tmp2, fftw_sfieldkx_tmp2, [ny_half, nzz_perproc, nxx])

        fftw_sfieldk = 0
        fftw_sfieldxx = 0
        fftw_sfieldkx_tmp1 = 0
        fftw_sfieldkx_tmp2 = 0
        
    end subroutine fftw_allocate

!==============================================================================

    subroutine fftw_init
     
        _indices
        
        call fftw_allocate()
        
        ! fftw plans
        
        ! forward, real to complex

        p_plan_y_forward = fftw_plan_many_dft_r2c(&
            1, [nyy], nzz_perproc * nxx, & 
            fftw_sfieldxx, [2 * nyy_half_pad1], 1, 2 * nyy_half_pad1, &
            fftw_sfieldkx, [nyy_half_pad1], 1, nyy_half_pad1 , &
            FFTW_ESTIMATE)

        p_plan_x_forward = fftw_plan_many_dft(&
            1, [nxx], ny_half * nzz_perproc, &
            fftw_sfieldkx_tmp2(:, :, :), [nxx], ny_half * nzz_perproc, 1, &
            fftw_sfieldkx_tmp2(:, :, :), [nxx], ny_half * nzz_perproc, 1, &
            FFTW_FORWARD, FFTW_ESTIMATE)

        p_plan_z_forward = fftw_plan_many_dft(&
            1, [nzz], ny_half * nx_perproc, &
            fftw_sfieldkx_tmp1(:, :, :), [nzz], ny_half * nx_perproc, 1, &
            fftw_sfieldkx_tmp1(:, :, :), [nzz], ny_half * nx_perproc, 1, &
            FFTW_FORWARD, FFTW_ESTIMATE)

        ! backward, complex to real

        p_plan_z_backward = fftw_plan_many_dft(&
            1, [nzz], ny_half * nx_perproc, &
            fftw_sfieldkx_tmp1(:, :, :), [nzz], ny_half * nx_perproc, 1, &
            fftw_sfieldkx_tmp1(:, :, :), [nzz], ny_half * nx_perproc, 1, &
            FFTW_BACKWARD, FFTW_ESTIMATE)

        p_plan_x_backward = fftw_plan_many_dft(&
            1, [nxx], ny_half * nzz_perproc, &
            fftw_sfieldkx_tmp2(:, :, :), [nxx], ny_half * nzz_perproc, 1, &
            fftw_sfieldkx_tmp2(:, :, :), [nxx], ny_half * nzz_perproc, 1, &
            FFTW_BACKWARD, FFTW_ESTIMATE)

        p_plan_y_backward = fftw_plan_many_dft_c2r(&
            1, [nyy], nzz_perproc * nxx, &
            fftw_sfieldkx, [nyy_half_pad1 ], 1, nyy_half_pad1 , &
            fftw_sfieldxx, [2 * nyy_half_pad1 ], 1, 2 * nyy_half_pad1 , &
            FFTW_ESTIMATE)

        ! wavenumbers

        ! kx are shared among processes
        do ix = 1, nx_perproc
            qx(ix) = my_id * nx_perproc + ix - 1
            if (qx(ix) > nx_half) then
                qx(ix) = qx(ix) - nx
            end if
            kx(ix) = (2.0_dp * PI / Lx)  * qx(ix)
        end do

        ! ky keeps only the non-negative modes (rest comes from hermitian symmetry)
        do iy = 1, ny_half
            qy(iy) = iy - 1
            ky(iy) = (2.0_dp * PI / Ly)  * qy(iy)
        end do
        
        ! kz
        do iz = 1, nz
            qz(iz) = iz - 1
            if (qz(iz) > nz_half) then 
                qz(iz) = qz(iz) - nz
            end if
            kz(iz) = (2.0_dp * PI / Lz)  * qz(iz)
        end do

        ! Find which index holds kF
        if (forcing /= 0) then
            do iy = 1, ny_half
                if (are_equal(ky(iy), kF)) then
                    iy_force = iy
                end if
            end do
        end if

        ! Find which index holds the zero mode in the distributed dimension
        do ix = 1, nx_perproc
            if (qx(ix) == 0) then
                ix_zero = ix
            else if (qx(ix) == 1) then
                ix_first_p = ix
            else if (qx(ix) == -1) then
                ix_first_n = ix
            end if
        end do

        ! Find which indices holds the maximal modes in x and z
        do ix = 1, nx_perproc
            if (qx(ix) == nx_half) then
                ix_max = ix
            end if
        end do

        ! broadcast which processor holds this index
        if (ix_max /= -1) my_id_ix_max1 = my_id
        call MPI_ALLREDUCE(my_id_ix_max1, my_id_ix_max, 1, MPI_INTEGER4, MPI_PROD, &
        MPI_COMM_WORLD, mpi_err)
        my_id_ix_max = abs(my_id_ix_max)

        do iz = 1, nz
            if (qz(iz) == nz_half) then
                iz_max = iz
            end if
        end do
            
    end subroutine fftw_init


!==============================================================================

    subroutine fftw_transpose(L,M,N,in,out)
        ! Adapted from nsCouette, mod_myMpi.f90, Florian Merz (IBM,RZG)
        ! https://gitlab.mpcdf.mpg.de/mjr/nscouette
    
        integer(i4),  intent(in)  :: L,M,N
        complex(dpc), intent(in)  :: in(L,M/num_procs,N)
        complex(dpc), intent(out) :: out(L,N/num_procs,M)
        complex(dpc) :: temp(L,M/num_procs,N/num_procs,num_procs)
        integer(i4)  :: blocksize, n_bl, pe, i, d1
        real(sp)     :: timer_start, timer_stop
    
        blocksize=L*(M/num_procs)*(N/num_procs)
        n_bl=M/num_procs

        ! WARNING: Might not work for when N=num_procs or M=num_procs
        ! parameters is set to abort at such values

        ! TODO: handle the above case.

        ! processor transpose (=alltoall), reordering and  local transpose
        call cpu_time(timer_start)
        call mpi_alltoall(in, blocksize, MPI_COMPLEX16,&
            temp, blocksize, MPI_COMPLEX16,&
            MPI_COMM_WORLD)
        call cpu_time(timer_stop)
        time_global_transpose = time_global_transpose + timer_stop - timer_start

        call cpu_time(timer_start)
        do pe=1,num_procs
            do i=1,N/num_procs
                do d1=1,L
                    out(d1,i,(pe-1)*n_bl+1:pe*n_bl)=temp(d1,:,i,pe)
                end do
            end do
        end do
        call cpu_time(timer_stop)
        time_local_transpose = time_local_transpose + timer_stop - timer_start
    
    end subroutine fftw_transpose

!==============================================================================

    subroutine fftw_transform(flag)

        integer(i4), intent(in) :: flag
        real(sp)                :: timer_start, timer_stop
        _indices
        
        if (flag == 1) then
            ! physical to spectral, end result is rescaled with norm_fft
            ! fftw_sfieldxx -> fftw_sfieldk

            ! do y: fftw_sfieldxx -> fftw_sfieldkx, in-place
            call cpu_time(timer_start)
            call fftw_execute_dft_r2c(p_plan_y_forward, &
                fftw_sfieldxx(:, :, :), fftw_sfieldkx(:, :, :))
            call cpu_time(timer_stop)
            time_fftw = time_fftw + timer_stop - timer_start

            ! shrink y: fftw_sfieldkx -> fftw_sfieldkx_tmp2
            call cpu_time(timer_start)
            do ix = 1, nxx
                do iz = 1, nzz_perproc
                    do iy = 1, ny_half
                        ! non-negative
                        fftw_sfieldkx_tmp2(iy, iz, ix) = fftw_sfieldkx(iy, iz, ix)
                    end do
                end do
            end do
            call cpu_time(timer_stop)
            time_local_transpose = time_local_transpose + timer_stop - timer_start

            ! do x: fftw_sfieldkx_tmp2 -> fftw_sfieldkx_tmp2, in-place
            call cpu_time(timer_start)
            call fftw_execute_dft(p_plan_x_forward, fftw_sfieldkx_tmp2, fftw_sfieldkx_tmp2)
            call cpu_time(timer_stop)
            time_fftw = time_fftw + timer_stop - timer_start

            ! shrink x: fftw_sfieldkx_tmp2 -> fftw_sfieldkx_tmp1_t
            call cpu_time(timer_start)
            ! maximum x mode is zeroed
            fftw_sfieldkx_tmp1_t(:, :, nx_half+1) = 0
            do ix = 1, nx_half 
                do iz = 1, nzz_perproc
                    do iy = 1, ny_half
                        ! non-negative
                        fftw_sfieldkx_tmp1_t(iy, iz, ix) = fftw_sfieldkx_tmp2(iy, iz, ix)
                        ! negative
                        if (ix <= nx_half - 1) then
                            fftw_sfieldkx_tmp1_t(iy, iz, nx_half + 1 + ix) &
                                = fftw_sfieldkx_tmp2(iy, iz, 1 + nxx - nx_half + ix)
                        end if
                    end do
                end do
            end do
            call cpu_time(timer_stop)
            time_local_transpose = time_local_transpose + timer_stop - timer_start

            ! transpose: fftw_sfieldkx_tmp1_t -> fftw_sfieldkx_tmp1, in-place with temporary arrays
            call fftw_transpose(ny_half, nzz, nx, fftw_sfieldkx_tmp1_t, fftw_sfieldkx_tmp1)

            ! do z: fftw_sfieldkx_tmp1 -> fftw_sfieldkx_tmp1, in-place
            call cpu_time(timer_start)
            call fftw_execute_dft(p_plan_z_forward, fftw_sfieldkx_tmp1, fftw_sfieldkx_tmp1)
            call cpu_time(timer_stop)
            time_fftw = time_fftw + timer_stop - timer_start

            ! shrink z: fftw_sfieldkx_tmp1 -> fftw_sfieldk
            call cpu_time(timer_start)
            ! maximum z mode is zeroed
            fftw_sfieldk(:, :, nz_half+1) = 0
            do iz = 1, nz_half
                do iy = 1, ny_half
                    do ix = 1, nx_perproc
                        if (ix_max /= -1 .and. ix == ix_max) cycle
                        ! non-negative
                        fftw_sfieldk(ix, iy, iz) = fftw_sfieldkx_tmp1(iy, ix, iz)
                        ! negative
                        if (iz <= nz_half - 1) then
                            fftw_sfieldk(ix, iy, nz_half + 1 + iz) &
                                = fftw_sfieldkx_tmp1(iy, ix, 1 + nzz - nz_half + iz)
                        end if
                    end do
                end do
            end do
            call cpu_time(timer_stop)
            time_local_transpose = time_local_transpose + timer_stop - timer_start

            ! normalize
            fftw_sfieldk(:, :, :) = fftw_sfieldk(:, :, :) * norm_fft

        else if (flag == -1) then

            ! spectral to physical
            ! fftw_sfieldk -> fftw_sfieldxx

            ! expand z: fftw_sfieldk -> fftw_sfieldkx_tmp1
            call cpu_time(timer_start)
            fftw_sfieldkx_tmp1(:, :, :) = 0
            do iz = 1, nz_half
                do iy = 1, ny_half
                    do ix = 1, nx_perproc
                        if (ix_max /= -1 .and. ix == ix_max) cycle
                        ! non-negative
                        fftw_sfieldkx_tmp1(iy, ix, iz) = fftw_sfieldk(ix, iy, iz)
                        ! negative
                        if (iz <= nz_half - 1) then
                            fftw_sfieldkx_tmp1(iy, ix, 1 + nzz - nz_half + iz) &
                                                = fftw_sfieldk(ix, iy, nz_half + 1 + iz)
                        end if
                    end do
                end do
            end do
            call cpu_time(timer_stop)
            time_local_transpose = time_local_transpose + timer_stop - timer_start
            
            ! do z: fftw_sfieldkx_tmp1 -> fftw_sfieldkx_tmp1, in-place
            call cpu_time(timer_start)
            call fftw_execute_dft(p_plan_z_backward, fftw_sfieldkx_tmp1, fftw_sfieldkx_tmp1)
            call cpu_time(timer_stop)
            time_fftw = time_fftw + timer_stop - timer_start

            ! transpose: fftw_sfieldkx_tmp1 -> fftw_sfieldkx_tmp1_t, in-place with temporary arrays
            call fftw_transpose(ny_half, nx, nzz, fftw_sfieldkx_tmp1, fftw_sfieldkx_tmp1_t)

            ! expand x: fftw_sfieldkx_tmp1_t -> fftw_sfieldkx_tmp2
            call cpu_time(timer_start)
            fftw_sfieldkx_tmp2(:, :, :) = 0
            do ix = 1, nx_half
                do iz = 1, nzz_perproc
                    do iy = 1, ny_half
                        ! non-negative
                        fftw_sfieldkx_tmp2(iy, iz, ix) = fftw_sfieldkx_tmp1_t(iy, iz, ix)
                        ! negative
                        if (ix <= nx_half - 1) then
                            ! x has one zeroed mode (maximum positive)
                            fftw_sfieldkx_tmp2(iy, iz, 1 + nxx - nx_half + ix) &
                                            = fftw_sfieldkx_tmp1_t(iy, iz, nx_half + 1 + ix)
                        end if
                    end do
                end do
            end do
            call cpu_time(timer_stop)
            time_local_transpose = time_local_transpose + timer_stop - timer_start

            ! do x: fftw_sfieldkx_tmp2 -> fftw_sfieldkx_tmp2, in-place
            call cpu_time(timer_start)
            call fftw_execute_dft(p_plan_x_backward, fftw_sfieldkx_tmp2, fftw_sfieldkx_tmp2)
            call cpu_time(timer_stop)
            time_fftw = time_fftw + timer_stop - timer_start

            ! expand y: fftw_sfieldkx_tmp2 -> fftw_sfieldkx
            call cpu_time(timer_start)
            fftw_sfieldkx(:, :, :) = 0
            do ix = 1, nxx
                do iz = 1, nzz_perproc
                    do iy = 1, ny_half
                        ! non-negative
                        fftw_sfieldkx(iy, iz, ix) = fftw_sfieldkx_tmp2(iy, iz, ix)
                    end do
                end do
            end do
            call cpu_time(timer_stop)
            time_local_transpose = time_local_transpose + timer_stop - timer_start

            ! do y: fftw_sfieldkx -> fftw_sfieldxx, in-place
            call cpu_time(timer_start)
            call fftw_execute_dft_c2r(p_plan_y_backward, &
                    fftw_sfieldkx(:, :, :), fftw_sfieldxx(:, :, :))
            call cpu_time(timer_stop)
            time_fftw = time_fftw + timer_stop - timer_start

        end if
        
    end subroutine fftw_transform

!==============================================================================

    subroutine fftw_sx2k(sfieldxx, sfieldk)
        
        real(dp), intent(in) :: sfieldxx(:, :, :)
        complex(dpc), intent(out) :: sfieldk(:, :, :)
        _indices
        _indicess

        _loop_phys_begin        
            fftw_sfieldxx(iyy, izz, ixx) = sfieldxx(iyy, izz, ixx)
        _loop_phys_end

        
        call fftw_transform(1)

        _loop_spec_begin
            sfieldk(ix, iy, iz) = fftw_sfieldk(ix, iy, iz)
        _loop_spec_end

    end subroutine fftw_sx2k

!==============================================================================

    subroutine fftw_sk2x(sfieldk, sfieldxx)
        
        complex(dpc), intent(in) :: sfieldk(:, :, :)
        real(dp), intent(out) :: sfieldxx(:, :, :)
        _indices
        _indicess
        
        _loop_spec_begin
            fftw_sfieldk(ix, iy, iz) = sfieldk(ix, iy, iz)
        _loop_spec_end

        call fftw_transform(-1)

        _loop_phys_begin     
            sfieldxx(iyy, izz, ixx) = fftw_sfieldxx(iyy, izz, ixx)
        _loop_phys_end

    end subroutine fftw_sk2x

!==============================================================================

    subroutine fftw_vx2k(sfieldxx, sfieldk)
        
        real(dp), intent(in) :: sfieldxx(:, :, :, :)
        complex(dpc), intent(out) :: sfieldk(:, :, :, :)
        integer(i4) :: n

        do n=1,3
            call fftw_sx2k(sfieldxx(:,:,:,n), sfieldk(:,:,:,n))
        end do

    end subroutine fftw_vx2k

!==============================================================================

    subroutine fftw_vk2x(sfieldk, sfieldxx)
        
        complex(dpc), intent(in) :: sfieldk(:, :, :, :)
        real(dp), intent(out) :: sfieldxx(:, :, :, :)
        integer(i4) :: n

        do n=1,3
            call fftw_sk2x(sfieldk(:,:,:,n), sfieldxx(:,:,:,n))
        end do

    end subroutine fftw_vk2x

!==============================================================================
    
end module fftw

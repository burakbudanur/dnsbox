program visualize
    ! modules:
    use m_numbers
    use m_openmpi
    use m_io
    use m_parameters
    use m_work
    use x_fftw
    use m_state
    use m_fields
    use m_stats
    use m_runs
    
    integer :: i
    
    ! initialization:
    call m_openmpi_init
    ! call m_openmpi_test
    call m_io_init
    call m_parameters_init
    call m_work_init
    call x_fftw_allocate(1)
    call x_fftw_init
    call m_state_init
    call m_fields_init
    call m_stats_init

    write(file_ext, "(i6.6)") 0
    fname = 'state.'//file_ext
    call m_work_read

    state(:, :, :, 1:3) = wrk(:, :, :, 1:3)
    fields(:, :, :, 1:3) = wrk(:, :, :, 1:3)
        
    call m_stats_compute
    call m_stats_write
    
    wrk(:, :, :, 1:3) = fields(:, :, :, 1:3)
    
    ! Convert to physical space:
    do i = 1, 3
        call xFFT3d(-1, i)
        wrk(:, :, :, i) = wrk(:, :, :, i) * normfac
    end do
    
    fname = 'velocity.'//file_ext
    call m_work_write_onedim

    ! State file is in state(:, :, :, 1:3)
    call m_state_vorticity
        
    wrk(:, :, :, 1:3) = state(:, :, :, 4:6)
    ! Convert to physical space:
    do i = 1, 3
        call xFFT3d(-1, i)
        wrk(:, :, :, i) = wrk(:, :, :, i) * normfac
    end do
    
    fname = 'vorticity.'//file_ext
    call m_work_write_onedim
    
    call m_runs_exit

end program visualize

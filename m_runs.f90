module m_runs
    
    ! modules:
    use m_numbers
    use m_openmpi
    use m_io
    use m_parameters
    use m_work
    use m_fields
    use x_fftw
    use m_state

    contains
    
    subroutine m_runs_init

        ! Initial time
        itime = ITMIN

        if (IC == 1) then
            call m_runs_read

            fields(:, :, :, 1:3) = wrk(:, :, :, 1:3)
            time = itime*dt
        elseif (IC == 0) then
            

            ! Random initial condition
            time = itime*dt
            call m_fields_init_random ! Mersenne twister
            
        end if

        call m_fields_dealias
        

    end subroutine m_runs_init
    
!==============================================================================

    subroutine m_runs_impose_symmetries

        if (Rz == 1) then
            
            state(:, :, :, 1:3) = fields(:, :, :, 1:3)
            call m_fields_real_reflect_z
            if (itime==ITMIN) then
                write(out, *) 'imposing Rz-symmetry on the field'
                flush(out)
            end if
            fields(:, :, :, 1:3) = 0.5d0 * fields(:, :, :, 1:3) &
                                + 0.5d0 * state(:, :, :, 1:3) 
        end if

        if (KolmSx == 1) then
            
            state(:, :, :, 1:3) = fields(:, :, :, 1:3)
            call m_fields_kolm_Sx
            if (itime==ITMIN) then
                write(out, *) 'imposing Kolm-Sx symmetry on the field'
                flush(out)
            end if
            fields(:, :, :, 1:3) = 0.5d0 * fields(:, :, :, 1:3) &
                                + 0.5d0 * state(:, :, :, 1:3) 
        end if

    end subroutine m_runs_impose_symmetries

!==============================================================================

    subroutine m_runs_exit

        if (statsWritten) then
            close(stat1FileCh)
        end if

        call x_fftw_allocate(-1)
        call m_io_exit    
        call m_openmpi_exit
        stop

    end subroutine m_runs_exit
    
!==============================================================================
    
    subroutine m_runs_set_fname
        
        write(file_ext, "(i6.6)") itime/IPRINT2
        fname = 'state.'//file_ext
    end subroutine m_runs_set_fname

!==============================================================================

    subroutine m_runs_write
        
        call m_runs_set_fname
        call m_work_write
    end subroutine m_runs_write

!==============================================================================

    subroutine m_runs_read
        
        call m_runs_set_fname
        call m_work_read
    end subroutine m_runs_read

end module m_runs
! Calculate the recurrence data from full state files
! Final state is decided with the laminarization time,
! initial state is set to state.000000

program recurrence
    
    ! modules:
    use m_numbers
    use m_openmpi
    use m_io
    use m_parameters
    use m_work
    use m_fields
    use x_fftw
    use m_state
    use m_stats
    use m_step
    use m_runs
    
    real(dp) :: recInitialTime, recFinalTime, stateiTime, statefTime
    integer :: statei, deltaState, statef, ishift, NRecur
    integer :: stateHead, stateTail, i
    real(dp)  :: d1, d2
    real(dp), allocatable :: states(:,:,:,:,:)
    real(dp), allocatable :: headState(:,:,:,:)
    real(sp), allocatable :: dist(:)

    character(*), parameter :: recFile   = 'recurrence.dat'
    
    integer :: un, recCh ! i/o channels
    adjoint = .false.
    
    ! initialization:
    call m_openmpi_init
    call m_io_init
    call m_parameters_init
    call m_work_init
    call m_fields_init
    call x_fftw_allocate(1)
    call x_fftw_init
    call m_state_init
    call m_step_init
    
    ! ! Read recurrence.in
    ! open(newunit=un,status='old',file='recurrence.in')
    ! read(un, *) stateiTime
    ! read(un, *) statefTime
    ! read(un, *) recInitialTime
    ! read(un, *) recFinalTime
    ! close(un)

    ! statei = FLOOR((stateiTime/dt)/IPRINT2, i8)

    ! More hardcoded version

    ! Read LAMINARIZED
    ! Note: This reads the very first line, so the Python script
    ! should also do so
    open(newunit=un,status='old',file='LAMINARIZED')
    read(un, *) statefTime
    close(un)

    ! Read recurrence.in
    open(newunit=un,status='old',file='recurrence.in')
    read(un, *) recInitialTime
    read(un, *) recFinalTime
    close(un)

    statei = 0

    deltaState = FLOOR((recInitialTime/dt)/IPRINT2, i8)
    NRecur = FLOOR((recFinalTime/dt)/IPRINT2, i8)
    statef = FLOOR((statefTime/dt)/IPRINT2, i8) - NRecur

    if (statei + deltaState >= statef) then
        write(*, *) 'Recurrence: Boundary state numbers are off, quitting.'
        call m_runs_exit
    end if
    
    if (myid == 0) allocate(dist(NRecur - deltaState + 1))
    allocate(states(nz+2,ny,nx,3,NRecur - deltaState + 1))
    allocate(headState(nz+2,ny,nx,3))

    ! Open the stream
    if (myid == 0) then
        open(newunit=recCh,file=recFile,status='replace',access="stream")
    end if

    ! run over all states
    do stateHead = statei, statef

        ! Load the head state
        itime = stateHead * IPRINT2
        call m_runs_read! Load a(t)
        headState(:, :, :, :) = wrk(:, :, :, 1:3)

        ! shift the tail
        ! if this is the first time, load the tail
        if (stateHead == statei) then
            do stateTail = statei + deltaState, statei + NRecur
                itime = stateTail * IPRINT2
                call m_runs_read! Load a(t)
                states(:, :, :, :,stateTail - statei - deltaState + 1) = wrk(:, :, :, 1:3)
            end do
        else
            ! shift the tail
            do i = 1, NRecur - deltaState
                states(:, :, :, :,i) = states(:, :, :, :,i+1)
            end do
            ! read the end of the tail
            itime = (stateHead + NRecur) * IPRINT2
            call m_runs_read ! Load a(t)
            states(:, :, :, :, NRecur - deltaState + 1) = wrk(:, :, :, 1:3)
        end if

        ! calculate recurrences

        state(:, :, :, 1:3) = headState(:, :, :, :)
        state(:, :, :, 4:6) = state(:, :, :, 1:3)
    
        call m_state_inprod(d1) ! |a(t)|^2

        do i = 1, NRecur - deltaState + 1
            state(:, :, :, 1:3) = states(:, :, :, 1:3,i) - headState(:, :, :, :)
            state(:, :, :, 4:6) = state(:, :, :, 1:3)

            call m_state_inprod(d2)
            if (myid == 0) dist(i) = real(dsqrt(d2/d1), sp)

        end do

        if (myid == 0) write(recCh) dist

    end do
    
    if (myid == 0) close(recCh)

    call m_runs_exit
        
end program recurrence


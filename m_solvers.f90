module m_solvers

    integer :: TxHalf_, TzHalf_, RxRy_

    contains

    subroutine m_solvers_relative_symmetries_read
        use m_io

        logical :: there

        ! Is relativeSymmetries.in here?
        inquire(file = 'relativeSymmetries.in', exist=there)
        if(.not.there) then
            write(out,*) 'cannot find relativeSymmetries.in'
            flush(out)
            stop
        end if
        
        ! Read parameters from the input file
        open(newunit=in, file='relativeSymmetries.in', form='formatted')
        read(in, *)
        read(in, *)
        read(in, *)

        write(out, '(79(''=''))')
        write(out, *) "Reading relative symmetries"
        write(out, '(79(''=''))')
                       
        read(in, *, ERR=9000, END=9000) TxHalf_
        write(out, *) 'TxHalf_ = ', TxHalf_

        read(in, *, ERR=9000, END=9000) TzHalf_
        write(out, *) 'TzHalf_ = ', TzHalf_

        read(in, *, ERR=9000, END=9000) RxRy_
        write(out, *) 'RxRy_ = ', RxRy_

        close(in)
        write(out,'(79(''=''))') 

        return
        !----------------------------------------------------------------------
        !  ERROR PROCESSING
        !----------------------------------------------------------------------

        9000 continue
            write(out,*)'An error was encountered while reading relativeSymmetries.in'
            flush(out)
            stop

    end subroutine m_solvers_relative_symmetries_read

!==============================================================================

    subroutine m_solvers_relative_symmetries_apply
        use m_fields
        
        if (TxHalf_ == 1) call m_fields_shift_x_half
    
        if (TzHalf_ == 1) call m_fields_shift_z_half

        if (RxRy_ == 1) call m_fields_real_kolm_RxRy

    end subroutine m_solvers_relative_symmetries_apply

end module m_solvers

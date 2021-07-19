#include "macros.h"
module stats
    use numbers
    use openmpi
    use io
    use parameters
    use fieldio
    use fftw
    use vfield
    use rhs
    use timestep

    real(dp) :: ekin, powerin, enstrophy, dissip, norm_rhs

    integer(i4) :: stats_stat_ch, stats_specx_ch, stats_specy_ch, &
        stats_specz_ch
    logical :: stats_stat_written = .false., stats_specs_written = .false.
    
    character(255) :: stats_stat_file = 'stat.gp', &
                      stats_specx_file = 'specs_x.gp', &
                      stats_specy_file = 'specs_y.gp', &
                      stats_specz_file = 'specs_z.gp'

    contains 

!==============================================================================

    subroutine stats_compute_powerin(vfieldk)
        complex(dpc), intent(in)  :: vfieldk(:, :, :, :)
        real(dp) :: my_powerin

        my_powerin = 0
        if (ix_zero /= -1) then
            if (tilting) then
                if (forcing == 1) then ! sine
                    my_powerin = -cos(tilt_angle * PI / 180.0_dp) * (PI**2 / (4.0_dp * Re)) &
                                    * vfieldk(ix_zero,iy_force,1,1)%im &
                                 -sin(tilt_angle * PI / 180.0_dp) * (PI**2 / (4.0_dp * Re)) &
                                    * vfieldk(ix_zero,iy_force,1,3)%im
                elseif (forcing == 2) then ! cosine
                    ! This may require thinking in the presence of drag
                    my_powerin = cos(tilt_angle * PI / 180.0_dp) * (PI**2 / (4.0_dp * Re)) &
                                                        * vfieldk(ix_zero,iy_force,1,1)%re  &
                                 + sin(tilt_angle * PI / 180.0_dp) * (PI**2 / (4.0_dp * Re)) &
                                                        * vfieldk(ix_zero,iy_force,1,3)%re
                end if
            else
                if (forcing == 1) then ! sine
                    my_powerin = -(PI**2 / (4.0_dp * Re)) * vfieldk(ix_zero,iy_force,1,1)%im
                elseif (forcing == 2) then ! cosine
                    ! This may require thinking in the presence of drag
                    my_powerin = (PI**2 / (4.0_dp * Re)) * vfieldk(ix_zero,iy_force,1,1)%re 
                end if
            end if
        end if

        call MPI_REDUCE(my_powerin, powerin, 1, MPI_REAL8, MPI_SUM, 0, &
        MPI_COMM_WORLD, mpi_err)

    end subroutine stats_compute_powerin

!==============================================================================

    subroutine stats_compute(vfieldk, fvfieldk)
        complex(dpc), intent(in), dimension(:, :, :, :)  :: &
             vfieldk, fvfieldk

        ! Kinetic energy
        call vfield_norm2(vfieldk, ekin, .false.)

        ! Power input
        call stats_compute_powerin(vfieldk)
                
        ! Dissipation
        call vfield_enstrophy(vfieldk, enstrophy, .false.)
        dissip = 2.0_dp * enstrophy / Re

        ! norm of rhs
        call vfield_norm(fvfieldk,norm_rhs,.false.)

    end subroutine stats_compute

!==============================================================================

    subroutine stats_spectra(vfieldk)
        complex(dpc), intent(in)  :: vfieldk(:, :, :, :)
        real(dp) :: specx(nx_half), my_specx(nx_half), &
                    specy(ny_half), my_specy(ny_half), &
                    specz(nz_half), my_specz(nz_half)
        complex(dpc) :: spec_
        character(255) :: formatStr

        _indices

        my_specx(:) = 0
        my_specy(:) = 0
        my_specz(:) = 0

        _loop_spec_begin
            spec_ = sum(conjg(vfieldk(ix,iy,iz,1:3))*vfieldk(ix,iy,iz,1:3))
            if (iy==1) spec_ = spec_ / 2
            my_specx(abs(qx(ix))+1) = my_specx(abs(qx(ix))+1) + spec_%re
            my_specy(abs(qy(iy))+1) = my_specy(abs(qy(iy))+1) + spec_%re
            my_specz(abs(qz(iz))+1) = my_specz(abs(qz(iz))+1) + spec_%re
        _loop_spec_end

        call MPI_REDUCE(my_specx, specx, nx_half, MPI_REAL8, MPI_SUM, 0, &
        MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(my_specy, specy, ny_half, MPI_REAL8, MPI_SUM, 0, &
        MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(my_specz, specz, nz_half, MPI_REAL8, MPI_SUM, 0, &
        MPI_COMM_WORLD, mpi_err)

        if (my_id==0) then
            
            ! outputting all this in the stat file
 
             inquire(file=TRIM(stats_specx_file), exist=there, opened=there2)
             if (.not.there) then
                open(newunit=stats_specx_ch,file=TRIM(stats_specx_file),form='formatted')
                 write(stats_specx_ch,"(A2,"//i4_len//","//"2"//sp_len//")") &
                     "# ", "itime", "time", "specs_x"
             end if
             if(there.and..not.there2) then
                open(newunit=stats_specx_ch,file=TRIM(stats_specx_file),position='append')
             end if
             write(formatStr,*) nx_half + 1
             write(stats_specx_ch,"(A2,"//i4_f//","//TRIM(formatStr)//sp_f//")")&
                 "  ", itime, time, specx

            inquire(file=TRIM(stats_specy_file), exist=there, opened=there2)
            if (.not.there) then
            open(newunit=stats_specy_ch,file=TRIM(stats_specy_file),form='formatted')
                write(stats_specy_ch,"(A2,"//i4_len//","//"2"//sp_len//")") &
                    "# ", "itime", "time", "specs_y"
            end if
            if(there.and..not.there2) then
            open(newunit=stats_specy_ch,file=TRIM(stats_specy_file),position='append')
            end if
            write(formatStr,*) 1 + ny_half
            write(stats_specy_ch,"(A2,"//i4_f//","//TRIM(formatStr)//sp_f//")")&
                "  ", itime, time, specy

            inquire(file=TRIM(stats_specz_file), exist=there, opened=there2)
            if (.not.there) then
            open(newunit=stats_specz_ch,file=TRIM(stats_specz_file),form='formatted')
                write(stats_specz_ch,"(A2,"//i4_len//","//"2"//sp_len//")") &
                    "# ", "itime", "time", "specs_z"
            end if
            if(there.and..not.there2) then
            open(newunit=stats_specz_ch,file=TRIM(stats_specz_file),position='append')
            end if
            write(formatStr,*) 1 + nz_half
            write(stats_specz_ch,"(A2,"//i4_f//","//TRIM(formatStr)//sp_f//")")&
                "  ", itime, time, specz
 
            stats_specs_written = .true.
 
         end if
    end subroutine stats_spectra

!==============================================================================

    subroutine stats_write
        
        ! outputting statistics
        
        if (my_id==0) then
            
           ! outputting all this in the stat file

            inquire(file=TRIM(stats_stat_file), exist=there, opened=there2)
            if (.not.there) then
            open(newunit=stats_stat_ch,file=TRIM(stats_stat_file),form='formatted')
                write(stats_stat_ch,"(A2,"//i4_len//","//"5"//sp_len//")") &
                    "# ", "itime", "time", "ekin", "powerin", "dissip", "norm_rhs"
            end if
            if(there.and..not.there2) then
            open(newunit=stats_stat_ch,file=TRIM(stats_stat_file),position='append')
            end if
            write(stats_stat_ch,"(A2,"//i4_f//","//"5"//sp_f//")")&
                "  ", itime, time, ekin, powerin, dissip, norm_rhs

           stats_stat_written = .true.

        end if

    end subroutine stats_write

!==============================================================================

    subroutine stats_worst_divergence(vfieldk)
        complex(dpc), intent(in)  :: vfieldk(:, :, :, :)
        complex(dpc) :: div_vfieldk(nx_perproc, ny_half, nz)
        real(dp)     :: div_sfieldxx(nyy, nzz_perproc, nxx)
        real(dp)     :: div, my_div
    
        call diffops_div(vfieldk, div_vfieldk)
        call fftw_sk2x(div_vfieldk, div_sfieldxx)
    
        my_div = maxval(abs(div_sfieldxx(:, :, :)))
        call MPI_REDUCE(my_div, div, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, mpi_err)
    
        if (my_id == 0 .and. div >  divergence_th) then
            write(out,*) 'stats: Time = ', time
            write(out,*) 'stats: Worst divergence: ', div
        end if
            
    end subroutine stats_worst_divergence    

end module stats
!==============================================================================!
!
!  Fast Fourier Transform that uses FFTW3 library,
!  pseudospectral DNS code
!  Copyright (C) 2006 Sergei Chumakov, Natalia Vladimirova, Misha Stepanov
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT any WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the
!    Free Software Foundation, Inc.,
!    51 Franklin Street, Fifth Floor,
!    Boston, MA 02110-1301, USA
!
!==============================================================================!

!                     GNU GENERAL PUBLIC LICENSE
!                     Version 2, June 1991

! Copyright (C) 1989, 1991 Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
! Everyone is permitted to copy and distribute verbatim copies
! of this license document, but changing it is not allowed.

!                          Preamble

! The licenses for most software are designed to take away your
! freedom to share and change it.  By contrast, the GNU General Public
! License is intended to guarantee your freedom to share and change free
! software--to make sure the software is free for all its users.  This
! General Public License applies to most of the Free Software
! Foundation's software and to any other program whose authors commit to
! using it.  (Some other Free Software Foundation software is covered by
! the GNU Lesser General Public License instead.)  You can apply it to
! your programs, too.

! When we speak of free software, we are referring to freedom, not
! price.  Our General Public Licenses are designed to make sure that you
! have the freedom to distribute copies of free software (and charge for
! this service if you wish), that you receive source code or can get it
! if you want it, that you can change the software or use pieces of it
! in new free programs; and that you know you can do these things.

! To protect your rights, we need to make restrictions that forbid
! anyone to deny you these rights or to ask you to surrender the rights.
! These restrictions translate to certain responsibilities for you if you
! distribute copies of the software, or if you modify it.

! For example, if you distribute copies of such a program, whether
! gratis or for a fee, you must give the recipients all the rights that
! you have.  You must make sure that they, too, receive or can get the
! source code.  And you must show them these terms so they know their
! rights.

! We protect your rights with two steps: (1) copyright the software, and
! (2) offer you this license which gives you legal permission to copy,
! distribute and/or modify the software.

! Also, for each author's protection and ours, we want to make certain
! that everyone understands that there is no warranty for this free
! software.  If the software is modified by someone else and passed on, we
! want its recipients to know that what they have is not the original, so
! that any problems introduced by others will not reflect on the original
! authors' reputations.

! Finally, any free program is threatened constantly by software
! patents.  We wish to avoid the danger that redistributors of a free
! program will individually obtain patent licenses, in effect making the
! program proprietary.  To prevent this, we have made it clear that any
! patent must be licensed for everyone's free use or not licensed at all.

! The precise terms and conditions for copying, distribution and
! modification follow.

!                  GNU GENERAL PUBLIC LICENSE
! TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

! 0. This License applies to any program or other work which contains
! a notice placed by the copyright holder saying it may be distributed
! under the terms of this General Public License.  The "Program", below,
! refers to any such program or work, and a "work based on the Program"
! means either the Program or any derivative work under copyright law:
! that is to say, a work containing the Program or a portion of it,
! either verbatim or with modifications and/or translated into another
! language.  (Hereinafter, translation is included without limitation in
! the term "modification".)  Each licensee is addressed as "you".

! Activities other than copying, distribution and modification are not
! covered by this License; they are outside its scope.  The act of
! running the Program is not restricted, and the output from the Program
! is covered only if its contents constitute a work based on the
! Program (independent of having been made by running the Program).
! Whether that is true depends on what the Program does.

! 1. You may copy and distribute verbatim copies of the Program's
! source code as you receive it, in any medium, provided that you
! conspicuously and appropriately publish on each copy an appropriate
! copyright notice and disclaimer of warranty; keep intact all the
! notices that refer to this License and to the absence of any warranty;
! and give any other recipients of the Program a copy of this License
! along with the Program.

! You may charge a fee for the physical act of transferring a copy, and
! you may at your option offer warranty protection in exchange for a fee.

! 2. You may modify your copy or copies of the Program or any portion
! of it, thus forming a work based on the Program, and copy and
! distribute such modifications or work under the terms of Section 1
! above, provided that you also meet all of these conditions:

!  a) You must cause the modified files to carry prominent notices
!  stating that you changed the files and the date of any change.

!  b) You must cause any work that you distribute or publish, that in
!  whole or in part contains or is derived from the Program or any
!  part thereof, to be licensed as a whole at no charge to all third
!  parties under the terms of this License.

!  c) If the modified program normally reads commands interactively
!  when run, you must cause it, when started running for such
!  interactive use in the most ordinary way, to print or display an
!  announcement including an appropriate copyright notice and a
!  notice that there is no warranty (or else, saying that you provide
!  a warranty) and that users may redistribute the program under
!  these conditions, and telling the user how to view a copy of this
!  License.  (Exception: if the Program itself is interactive but
!  does not normally print such an announcement, your work based on
!  the Program is not required to print an announcement.)

! These requirements apply to the modified work as a whole.  If
! identifiable sections of that work are not derived from the Program,
! and can be reasonably considered independent and separate works in
! themselves, then this License, and its terms, do not apply to those
! sections when you distribute them as separate works.  But when you
! distribute the same sections as part of a whole which is a work based
! on the Program, the distribution of the whole must be on the terms of
! this License, whose permissions for other licensees extend to the
! entire whole, and thus to each and every part regardless of who wrote it.

! Thus, it is not the intent of this section to claim rights or contest
! your rights to work written entirely by you; rather, the intent is to
! exercise the right to control the distribution of derivative or
! collective works based on the Program.

! In addition, mere aggregation of another work not based on the Program
! with the Program (or with a work based on the Program) on a volume of
! a storage or distribution medium does not bring the other work under
! the scope of this License.

! 3. You may copy and distribute the Program (or a work based on it,
! under Section 2) in object code or executable form under the terms of
! Sections 1 and 2 above provided that you also do one of the following:

!  a) Accompany it with the complete corresponding machine-readable
!  source code, which must be distributed under the terms of Sections
!  1 and 2 above on a medium customarily used for software interchange; or,

!  b) Accompany it with a written offer, valid for at least three
!  years, to give any third party, for a charge no more than your
!  cost of physically performing source distribution, a complete
!  machine-readable copy of the corresponding source code, to be
!  distributed under the terms of Sections 1 and 2 above on a medium
!  customarily used for software interchange; or,

!  c) Accompany it with the information you received as to the offer
!  to distribute corresponding source code.  (This alternative is
!  allowed only for noncommercial distribution and only if you
!  received the program in object code or executable form with such
!  an offer, in accord with Subsection b above.)

! The source code for a work means the preferred form of the work for
! making modifications to it.  For an executable work, complete source
! code means all the source code for all modules it contains, plus any
! associated interface definition files, plus the scripts used to
! control compilation and installation of the executable.  However, as a
! special exception, the source code distributed need not include
! anything that is normally distributed (in either source or binary
! form) with the major components (compiler, kernel, and so on) of the
! operating system on which the executable runs, unless that component
! itself accompanies the executable.

! If distribution of executable or object code is made by offering
! access to copy from a designated place, then offering equivalent
! access to copy the source code from the same place counts as
! distribution of the source code, even though third parties are not
! compelled to copy the source along with the object code.

! 4. You may not copy, modify, sublicense, or distribute the Program
! except as expressly provided under this License.  Any attempt
! otherwise to copy, modify, sublicense or distribute the Program is
! void, and will automatically terminate your rights under this License.
! However, parties who have received copies, or rights, from you under
! this License will not have their licenses terminated so long as such
! parties remain in full compliance.

! 5. You are not required to accept this License, since you have not
! signed it.  However, nothing else grants you permission to modify or
! distribute the Program or its derivative works.  These actions are
! prohibited by law if you do not accept this License.  Therefore, by
! modifying or distributing the Program (or any work based on the
! Program), you indicate your acceptance of this License to do so, and
! all its terms and conditions for copying, distributing or modifying
! the Program or works based on it.

! 6. Each time you redistribute the Program (or any work based on the
! Program), the recipient automatically receives a license from the
! original licensor to copy, distribute or modify the Program subject to
! these terms and conditions.  You may not impose any further
! restrictions on the recipients' exercise of the rights granted herein.
! You are not responsible for enforcing compliance by third parties to
! this License.

! 7. If, as a consequence of a court judgment or allegation of patent
! infringement or for any other reason (not limited to patent issues),
! conditions are imposed on you (whether by court order, agreement or
! otherwise) that contradict the conditions of this License, they do not
! excuse you from the conditions of this License.  If you cannot
! distribute so as to satisfy simultaneously your obligations under this
! License and any other pertinent obligations, then as a consequence you
! may not distribute the Program at all.  For example, if a patent
! license would not permit royalty-free redistribution of the Program by
! all those who receive copies directly or indirectly through you, then
! the only way you could satisfy both it and this License would be to
! refrain entirely from distribution of the Program.

! If any portion of this section is held invalid or unenforceable under
! any particular circumstance, the balance of the section is intended to
! apply and the section as a whole is intended to apply in other
! circumstances.

! It is not the purpose of this section to induce you to infringe any
! patents or other property right claims or to contest validity of any
! such claims; this section has the sole purpose of protecting the
! integrity of the free software distribution system, which is
! implemented by public license practices.  Many people have made
! generous contributions to the wide range of software distributed
! through that system in reliance on consistent application of that
! system; it is up to the author/donor to decide if he or she is willing
! to distribute software through any other system and a licensee cannot
! impose that choice.

! This section is intended to make thoroughly clear what is believed to
! be a consequence of the rest of this License.

! 8. If the distribution and/or use of the Program is restricted in
! certain countries either by patents or by copyrighted interfaces, the
! original copyright holder who places the Program under this License
! may add an explicit geographical distribution limitation excluding
! those countries, so that distribution is permitted only in or among
! countries not thus excluded.  In such case, this License incorporates
! the limitation as if written in the body of this License.

! 9. The Free Software Foundation may publish revised and/or new versions
! of the General Public License from time to time.  Such new versions will
! be similar in spirit to the present version, but may differ in detail to
! address new problems or concerns.

! Each version is given a distinguishing version number.  If the Program
! specifies a version number of this License which applies to it and "any
! later version", you have the option of following the terms and conditions
! either of that version or of any later version published by the Free
! Software Foundation.  If the Program does not specify a version number of
! this License, you may choose any version ever published by the Free Software
! Foundation.

! 10. If you wish to incorporate parts of the Program into other free
! programs whose distribution conditions are different, write to the author
! to ask for permission.  For software which is copyrighted by the Free
! Software Foundation, write to the Free Software Foundation; we sometimes
! make exceptions for this.  Our decision will be guided by the two goals
! of preserving the free status of all derivatives of our free software and
! of promoting the sharing and reuse of software generally.

!                          NO WARRANTY

! 11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
! FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
! OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
! PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
! OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
! TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
! PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
! REPAIR OR CORRECTION.

! 12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
! WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
! REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
! INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
! OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
! TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
! YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
! PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGES.

!                   END OF TERMS AND CONDITIONS

!          How to Apply These Terms to Your New Programs

! If you develop a new program, and you want it to be of the greatest
! possible use to the public, the best way to achieve this is to make it
! free software which everyone can redistribute and change under these terms.

! To do so, attach the following notices to the program.  It is safest
! to attach them to the start of each source file to most effectively
! convey the exclusion of warranty; and each file should have at least
! the "copyright" line and a pointer to where the full notice is found.

!  <one line to give the program's name and a brief idea of what it does.>
!  Copyright (C) <year>  <name of author>

!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.

!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.

!  You should have received a copy of the GNU General Public License along
!  with this program; if not, write to the Free Software Foundation, Inc.,
!  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

! Also add information on how to contact you by electronic and paper mail.

! If the program is interactive, make it output a short notice like this
! when it starts in an interactive mode:

!  Gnomovision version 69, Copyright (C) year name of author
!  Gnomovision comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
!  This is free software, and you are welcome to redistribute it
!  under certain conditions; type `show c' for details.

! The hypothetical commands `show w' and `show c' should show the appropriate
! parts of the General Public License.  Of course, the commands you use may
! be called something other than `show w' and `show c'; they could even be
! mouse-clicks or menu items--whatever suits your program.

! You should also get your employer (if you work as a programmer) or your
! school, if any, to sign a "copyright disclaimer" for the program, if
! necessary.  Here is a sample; alter the names:

! Yoyodyne, Inc., hereby disclaims all copyright interest in the program
! `Gnomovision' (which makes passes at compilers) written by James Hacker.

! <signature of Ty Coon>, 1 April 1989
! Ty Coon, President of Vice

! This General Public License does not permit incorporating your program into
! proprietary programs.  If your program is a subroutine library, you may
! consider it more useful to permit linking proprietary applications with the
! library.  If this is what you want to do, use the GNU Lesser General
! Public License instead of this License.

module x_fftw
    use m_numbers
    use m_parameters
    use m_io
    use m_work
    
    
    
    ! FFTW parameters that do not change
    integer(i8), parameter :: FFTW_ESTIMATE =  0
    integer(i8), parameter :: FFTW_FORWARD  = -1
    integer(i8), parameter :: FFTW_BACKWARD =  1
    
    ! dedicated arrays for parallel FFT
    real(dp), allocatable :: zy_sheet(:, :), buff(:, :, :), x_stick(:)
    
    real(dp), allocatable :: buff2(:, :, :)
    
    ! order of message passing between processors
    integer, allocatable :: order(:), order_matrix(:, :)
    
    ! the arrays to store FFTW plans for the 2D r2c or c2r steps
    ! plan_r2c(1...nx, 1..ubound(wrk,4)) ... 
    integer(i8), allocatable :: plan_r2c(:, :), plan_c2r(:, :)
    
    ! FFTW plans for the 1D c2c forward/backward steps
    integer(i8) :: plan_f_c2c, plan_b_c2c
    
    ! k-vectors (real arrays)
    real(dp), allocatable :: akx(:), aky(:), akz(:)
    
    ! auxiliary parameters
    integer     :: nz21
    real(dp)    :: normfac
    
    ! array that contains indicator of aliasing when products are taken
    integer, allocatable :: ialias(:,:,:)
    real(dp), allocatable :: dealias_matrix(:,:,:) 
    
    ! error indicator:
    integer :: ierr
    
    contains 
    !==========================================================================
    !   SUBROUTINES
    !==========================================================================
    
    subroutine x_fftw_allocate(flag)
        !======================================================================
        ! Allocate(flag = 1)/Deallocate(flag = -1) FFTW arrays
        !======================================================================
                
        
        integer :: flag
        
        ierr = 0 
        if (flag == 1) then
            ! Allocate FFTW arrays
            allocate(&
                plan_r2c(nx, LBOUND(wrk,4):UBOUND(wrk,4)), &
                plan_c2r(nx, LBOUND(wrk,4):UBOUND(wrk,4)), &
                zy_sheet(nz, ny), buff(nz + 2, nx, nx), &
                x_stick(2 * nx_all), akz(nz + 2), aky(nx), akx(ny), &
                order(numprocs - 1), &
                buff2(nz + 2, nx, nx), &
                ialias(nz + 2, ny, nx), dealias_matrix(nz+2, ny, nx), stat = ierr)
                
            if (ierr /= 0) then
                write(out,*) 'X_FFTW_ALLOCATE: cannot allocate'
                flush(out)
                stop
            end if
            
            write(out,*) 'X_FFTW allocated'
            
            
            ! assign temporary values to allocated arrays
            plan_r2c = 0
            plan_c2r = 0
            zy_sheet = zero
            buff     = zero
            buff2    = zero
            x_stick  = zero
            order    = 0
            
        elseif (flag == -1) then
            
            if (allocated(plan_r2c)) then
                deallocate(plan_r2c, plan_c2r, &
                           zy_sheet, buff, x_stick, order, akx, aky, akz, &
                           buff2, ialias)
            end if
            write(out, *) "X_FFTW deallocated"
        
        else
            
            write (out, *) 'X_FFTW_ALLOCATE : Wrong value of flag' ,flag
            flush(out)
            stop
            
        end if
    
        flush(out)
        return
        
    end subroutine x_fftw_allocate
    
    subroutine x_fftw_init
        !======================================================================
        ! Initialize auxiliary arrays for FFT
        !======================================================================
        
        
        integer :: itmp, ix, iy, iz, n, i, j, k
        
        ierr = 0
        !----------------------------------------------------------------------
        ! fill up order_matrix
        !----------------------------------------------------------------------
        allocate(order_matrix(numprocs, numprocs), stat = ierr)
        if (ierr /= 0 ) then
            write(out, *) 'x_fftw_init: cannot allocate order_matrix'
            flush(out)
            stop
        end if
        
        order_matrix(1, 1) = 0
        itmp = 1
        
        do while (itmp < numprocs)
            
            do ix = 1, itmp
                do iy = 1, itmp
                    
                    order_matrix(ix + itmp, iy) = order_matrix(ix, iy) + itmp
                    order_matrix(ix, iy + itmp) = order_matrix(ix, iy) + itmp
                    order_matrix(ix + itmp, iy + itmp) = order_matrix(ix, iy)
                    
                end do
            end do
            
            itmp = 2 * itmp
            
        end do
        
        write (out, *) "printing order_matrix"
        
        ! Print order_matrix
        do i = 1, ubound(order_matrix, 1)
            
            write (out, *) (order_matrix(i, j), j = 1, ubound(order_matrix, 2))
            
            
        end do
        
        
        !----------------------------------------------------------------------
        !   fill order and deallocate order_matrix
        !----------------------------------------------------------------------        
        ! order stores numprocs-1 elements:
        ! myid + 1 % numprocs, 
        
        do ix = 1, numprocs
            if (order_matrix(ix, myid + 1) /= 0) then
                order(order_matrix(ix, myid + 1)) = ix - 1
            end if            
        end do

        write (out, *) "printing order"
        
        write (out, *) (order(j), j = 1, ubound(order, 1))
        
        
        deallocate(order_matrix)
        
        !----------------------------------------------------------------------
        !   initializing FFTW plans for 1st step in FFT  : 2D r2c
        !----------------------------------------------------------------------
        do n=lbound(wrk,4), ubound(wrk, 4)
            do iz = 1, nx
                call dfftw_plan_dft_r2c_2d (plan_r2c(iz, n), nz, ny, &
                                            zy_sheet, wrk(1, 1, iz, n), &
                                            FFTW_ESTIMATE)
            end do
        end do
        
        !----------------------------------------------------------------------
        !   initializing FFTW plans for 2nd step in FFT  : 1D c2c
        !----------------------------------------------------------------------
        
        call dfftw_plan_dft_1d (plan_f_c2c, nx_all, x_stick, x_stick, &
                                FFTW_FORWARD, FFTW_ESTIMATE)
                                
        !----------------------------------------------------------------------
        !   initializing FFTW plans for 1st step in IFFT : 1D c2c
        !----------------------------------------------------------------------        
        
        call dfftw_plan_dft_1d (plan_b_c2c, nx_all, x_stick, x_stick, &
                                FFTW_BACKWARD, FFTW_ESTIMATE)
                                
        !----------------------------------------------------------------------
        !   initializing FFTW plans for 2nd step in IFFT : 2D c2r
        !----------------------------------------------------------------------        
        
        do n = lbound(wrk, 4), ubound(wrk, 4)
            do iz = 1, nx
                
                call dfftw_plan_dft_c2r_2d (plan_c2r(iz, n), nz, ny, &
                                            wrk(1, 1, iz, n), zy_sheet, &
                                            FFTW_ESTIMATE)
                                
            end do          
        end do

        !---------------------------------------------------------------------!
        !  initialize constants
        !---------------------------------------------------------------------!
        
        nz21 = nz / 2 + 1
        normfac = one / real(nz * ny * nx_all, 8)
        
        !---------------------------------------------------------------------!
        !  fill up wavenumber arrays
        !---------------------------------------------------------------------!
        
        ! in Fourier space it is (nz / 2 + 1) complex numbers along kz-axis
        do ix = 1, nz + 1, 2
            
            akz (ix)       = real((ix - 1) / 2, 8) * (2.0d0 / Lz)
            akz (ix + 1)   = akz(ix)
            
        end do
                
        ! in Fourier space ky-axis is distributed among the processors
        do iy = 1, nx
            aky(iy)    = real(myid * nx + iy - 1, 8)
            if (aky(iy) > (0.5D0 * real(nx_all, 8))) aky(iy) = aky(iy) &
                                                         - real(nx_all, 8)
            
        end do 
        
        ! in Fourier space x wave numbers are aligned along the second index
        do iz = 1, ny
        
            akx (iz)   = real(iz - 1, 8)
            if (akx(iz) > (0.5d0 * real(ny, 8))) akx(iz) = akx(iz) &
                                                             - real(ny, 8)
        
        end do

        ! Define ialias
        ialias(:, :, :) = 0
        dealias_matrix(:, :, :) = one
        do k = 1, nx
            do j = 1, ny
                do i = 1, nz + 2                                                            
                    if  (( abs(akx(j)) > rnx &
                      .or. abs(aky(k)) > rny &
                      .or. abs(akz(i)) > rnz ) &
                      .or. ( abs(akx(j)) == zero &
                      .and. abs(aky(k)) == zero &
                      .and. abs(akz(i)) == zero )) then
                        ialias(i, j, k) = 1
                        dealias_matrix(i,j,k) = zero
                    end if
                end do                
            end do
        end do
            
        write(out, *) "x_fftw arrays are initialized."
        flush(out)
            
        return
    end subroutine x_fftw_init
    
    
    subroutine xFFT3d(flag, n)
        !=====================================================================!
        !   Subroutine that performs the FFT of a 3-D variable.  The variable 
        !   is contained within the array "wrk(:, :, :, n)".  Note that the
        !   result of FFT has different coordinate arrangement: in physical
        !   space it is (x, y, z), and in Fourier space it is (kx, kz, ky).
        !   Details can be extracted from very graphic comments in the body of
        !   the subroutine.
        !=====================================================================!
        
        use m_openmpi
         
        
        integer :: flag, n
        integer :: ix, iy, iz
        
        real(dp) :: rtmp
        
        integer :: iproc
        
        if (flag == 1) then
        ! Forward FFT
        
!------------------------------------------------------------------------------!
!  Direct FFT, step 1: 2-D real-to-complex transform
!------------------------------------------------------------------------------!
!  
!   R2C           (# = ny) A y          A y                 (# = ny) A k_y
!                          |            |                            |
!          +---+---+---+---+            +            +---+---+---+---+
!         /   /   /   /   /|           /|           /   /   /   /   /|
!        /   /   /   / wrk         zy_sheet        /   /   /   / wrk
!       /   /   /   /   /  |         /  |         /   /   /   /   /  |
!      +---+---+---+---+   |        +   |        +---+---+---+---+   |
!      |   |   |   |   |   |        |   |  R2C   |   |   |   |   |   |
!      |   |   |   |   |   | -----> |   | -----> |   |   |   |   |   |
!      |   |   |   |   |   |        |   |        |   |   |   |   |   |
!  <---|   |   |   |   |   +        |   +    <---|   |   |   |   |   +
!  x   |   |   |   |   |  /         |  /     x   |   |   |   |   |  /
!      |   |   |   |   | /          | /          |   |   |   |   | /
!      |   |   |   |   |/           |/           |   |   |   |   |/
!      +---+---+---+---+            +            +---+---+---+---+
!                     /            /                            /
!           (# = nz) V z          V z         (# = nz / 2 + 1) V k_z
!
!        3   2   1   0  --- myid                   3   2   1   0  --- myid
!
!------------------------------------------------------------------------------!        
            do iz = 1, nx
                do iy = 1, ny
                    do ix = 1, nz
                        zy_sheet(ix, iy) = wrk(ix, iy, iz, n)
                    end do
                end do
                
                call dfftw_execute(plan_r2c(iz, n))
            end do
        
!------------------------------------------------------------------------------!
!  Direct FFT, step 2:  transposing the variable via MPI messaging
!------------------------------------------------------------------------------!
!
!   MPI           (# = ny) A k_y   +---+          +---+          (# = nx_all) A x
!                          |      /buff|         /buff|                   |
!          +---+---+---+---+     /   / +        /   / +   +---+---+---+---+
!         /   /+++/   /   /|    /   / /  MPI   /   / /   /   /   /   /   /|
!        /   /+++/   /  wrk    +---+ /  ----> +---+ /   / ../.  /   / wrk
!       /   /+++/   /   /  |   |   |/         |   |/   / ../.. /   /   /  |
!      +---+---+---+---+   |   +---+          +---+   +---+---+---+---+   |
!      |   |+++|   |   |   |    A                \    |...|.  |   |   |   |
!      |   |   \____   |   ____/                  `--->+++|   |   |   |   |
!      |   |   |   |\_____/|                          |+++|   |   |   |   |
!  <---|   |   |   |   |   +                      <---|   |   |   |   |   +
!  x   |   |   |   |   |  /                       k_y |   |   |   |   |  /
!  (nx_all / numprocs) | /                        (# = ny / numprocs) | /
!      |   |   |   |   |/                             |   |   |   |   |/
!      +---+---+---+---+                              +---+---+---+---+
!                     /                                              /
!                    V k_z                                          V k_z
!
!        3   2   1   0  --- myid                        3   2   1   0  --- myid
!
!------------------------------------------------------------------------------!
            ! - "diagonal" messages, no need to use MPI

!!$          write(out,*) 'before diagonal message'
!!$                  
            do iz = 1, nx - 1
                do iy = iz + 1, nx
                    do ix = 1, nz + 2
                        
                        rtmp = wrk(ix, myid * nx + iy, iz, n)
                        wrk(ix, myid * nx + iy, iz, n) = &
                            & wrk(ix, myid * nx + iz, iy, n)
                        wrk(ix, myid * nx + iz, iy, n) = rtmp
                    
                    end do
                end do
            end do
            
            ! - sending and receiving MPI messages (off-diagonal elements)
!------------------------------------------------------------------------------
! Elements of the order matrix encodes which blocks are going to be exchanged
! at each iteration of the following loop. 
! For 4 processors, order matrix looks like this:
!            
!                          y     
!       +---+---+---+---+               
!       | 0 | 1 | 2 | 3 |  1
!       | 1 | 0 | 3 | 2 |  2
!       | 2 | 3 | 0 | 1 |  3
!       | 3 | 2 | 1 | 0 |  4
!       +---+---+---+---+  
!     x   1   2   3   4
! procid  0   1   2   3
!
! elements of the order array encodes communication destination at each step
! order(iproc) = process to be communicated at step iproc
! 
! In total, transposition is complete numprocs - 1 communications
!
!------------------------------------------------------------------------------

            do iproc = 1, numprocs - 1
                ! Fill the buffer with blocks of data to be exchanged:
                do iy = 1, nx
                    do iz = 1, nx
                        buff(:, iz, iy) = wrk(:, order(iproc) * nx + iy, iz, n)
                    end do
                end do
                
                ! Communication:         !DATA  !# of elements      !DTYPE
                call MPI_SENDRECV_REPLACE(buff, (nz+2) * nx * nx, MPI_REAL8,&
                        & order(iproc), &                         ! destination
                        & myid * numprocs + order(iproc), &       ! send tag
                        & order(iproc), &                         ! source
                        & order(iproc) * numprocs + myid, &       ! receive tag
                        & MPI_COMM_WORLD, mpi_status_var, mpi_err)
                
                ! Fill the wrk with transposed data:
                do iy = 1, nx
                    do iz = 1, nx
                        wrk(:, order(iproc) * nx + iz, iy, n) = buff(:, iz, iy)
                    end do
                end do
            end do
!            
! Note: hit3d's xFFT3d uses mpi_sendrecv rather than mpi_sendrecv_replace
!
!------------------------------------------------------------------------------!
!  Direct FFT, step 3: one-dimensional complex-to-complex FFT
!------------------------------------------------------------------------------!
!
!   C2C->                  A x                                            A k_x
!                          |        A x      A k_x                        |
!          +---+---+---+---+        |        |            +---+---+---+---+
!         /   /   /   /   /|        +        +           /   /   /   /   /|
!        /   /   /   /  wrk         |        |          /   /   /   / wrk
!       /   /   /   /   /  |     x_stick  x_stick      /   /   /   /   /  |
!      +---+---+---+---+   |        |        |        +---+---+---+---+   |
!      |   |   |   |   |   |        |  C2C   |        |   |   |   |   |   |
!      |   |   |   |   |   | -----> | -----> | -----> |   |   |   |   |   |
!      |   |   |   |   |   |        |        |        |   |   |   |   |   |
!  <---|   |   |   |   |   +        |        |    <---|   |   |   |   |   +
!  k_y |   |   |   |   |  /         +        +    k_y |   |   |   |   |  /
!      |   |   |   |   | /                            |   |   |   |   | /
!      |   |   |   |   |/                             |   |   |   |   |/
!      +---+---+---+---+                              +---+---+---+---+
!                     /                                              /
!                    V k_z                                          V k_z
!  
!        3   2   1   0  --- myid                        3   2   1   0  --- myid
!  
!------------------------------------------------------------------------------!

            do iy = 1, nx
                do ix = 1, nz21
                    do iz = 1, nx_all
                        x_stick(2 * iz - 1) = wrk(2 * ix - 1, iz, iy, n) ! Real
                        x_stick(2 * iz) = wrk(2 * ix, iz, iy, n)         ! Imag
                    end do
                    
                    call dfftw_execute(plan_f_c2c)
                    
                    do iz = 1, nx_all
                        wrk(2 * ix - 1, iz, iy, n) = x_stick(2 * iz - 1) ! Real
                        wrk(2 * ix, iz, iy, n)     = x_stick(2 * iz)     ! Imag
                    end do
                end do
            end do
            
!------------------------------------------------------------------------------
! End of forward fft
!------------------------------------------------------------------------------            
        
        elseif (flag == -1) then

            ! Do dealising when transforming Fourier to real / g 200805
            call x_dealias(n)

!-----------------------------------------------------------------------------!
!  Inverse FFT, step 1: one-dimensionsal complex-to-complex transform
!-----------------------------------------------------------------------------!
!   
!   <-C2C                  A k_x                                          A x
!                          |        A k_x    A x                          |
!          +---+---+---+---+        |        |            +---+---+---+---+
!         /   /   /   /   /|        +        +           /   /   /   /   /|
!        /   /   /   / wrk          |        |          /   /   /   / wrk
!       /   /   /   /   /  |     x_stick  x_stick      /   /   /   /   /  |
!      +---+---+---+---+   |        |        |        +---+---+---+---+   |
!      |   |   |   |   |   |        |  C2C   |        |   |   |   |   |   |
!      |   |   |   |   |   | -----> | -----> | -----> |   |   |   |   |   |
!      |   |   |   |   |   |        |        |        |   |   |   |   |   |
!  <---|   |   |   |   |   +        |        |    <---|   |   |   |   |   +
!  k_y |   |   |   |   |  /         +        +    k_y |   |   |   |   |  /
!      |   |   |   |   | /                            |   |   |   |   | /
!      |   |   |   |   |/                             |   |   |   |   |/
!      +---+---+---+---+                              +---+---+---+---+
!                     /                                              /
!                    V k_z                                          V k_z
!  
!        3   2   1   0  --- myid                        3   2   1   0  --- myid
!
!------------------------------------------------------------------------------!
            
            do iy = 1, nx
                do ix = 1, nz21
                    do iz = 1, nx_all
                        x_stick(2 * iz - 1) = wrk(2 * ix - 1, iz, iy, n) ! Real
                        x_stick(2 * iz)     = wrk(2 * ix, iz, iy, n)     ! Imag    
                    end do
                    
                    call dfftw_execute (plan_b_c2c)
                    
                    do iz = 1, nx_all
                        wrk(2 * ix - 1, iz, iy, n) = x_stick(2 * iz - 1)
                        wrk(2 * ix, iz, iy, n) = x_stick(2 * iz)
                    end do
                end do
            end do            
        
!-----------------------------------------------------------------------------!
!  Inverse FFT, step 2: transposing the variable via MPI messaging
!-----------------------------------------------------------------------------!
!
!   MPI       (# = nx_all) A x     +---+          +---+          (# = ny) A k_y
!                          |      /buff|         /buff|                   |
!          +---+---+---+---+     /   / +        /   / +   +---+---+---+---+
!         /   /+++/   /   /|    /   / /  MPI   /   / /   /   /   /   /   /|
!        /   /+++/   / wrk |   +---+ /  ----> +---+ /   / ../.  /   / wrk
!       /   /+++/   /   /  |   |   |/         |   |/   / ../.. /   /   /  |
!      +---+---+---+---+   |   +---+          +---+   +---+---+---+---+   |
!      |   |+++|   |   |   |    A                \    |...|.  |   |   |   |
!      |   |   \____   |   ____/                  `--->+++|   |   |   |   |
!      |   |   |   |\_____/|                          |+++|   |   |   |   |
!  <---|   |   |   |   |   +                      <---|   |   |   |   |   +
!  k_y |   |   |   |   |  /                       x   |   |   |   |   |  /
!  (# = ny / numprocs) | /                    (# = nx_all / numprocs) | /
!      |   |   |   |   |/                             |   |   |   |   |/
!      +---+---+---+---+                              +---+---+---+---+
!                     /                                              /
!                    V k_z                                          V k_z
!
!        3   2   1   0  --- myid                        3   2   1   0  --- myid
!
!-----------------------------------------------------------------------------!

            ! - "diagonal" messages, no need to use MPI
            do iy = 1, nx - 1
                do iz = iy + 1, nx
                    do ix = 1, nz + 2
                        rtmp = wrk(ix, myid * nx + iz, iy, n)
                        wrk(ix, myid * nx + iz , iy, n) = &
                                              & wrk(ix, myid * nx + iy, iz, n)
                        wrk(ix, myid * nx + iy, iz, n) = rtmp
                    end do
                end do
            end do
       
            ! - sending and receiving MPI messages
            do iproc = 1, numprocs - 1
                do iz = 1, nx
                    do iy = 1, nx
                        buff(:, iy, iz) = wrk(:, order(iproc) * nx + iz, iy, n)
                    end do
                end do
          
                call MPI_SENDRECV_REPLACE(buff, &
                   & (nz + 2) * nx * nx, MPI_REAL8, &
                   & order(iproc), myid * numprocs + order(iproc), &
                   & order(iproc), order(iproc) * numprocs + myid, &
                   & MPI_COMM_WORLD, mpi_status_var, mpi_err) 
                
                do iz = 1, nx
                    do iy = 1, nx
                        wrk(:, order(iproc) * nx + iy, iz, n) = buff(:, iy, iz)
                    end do
                end do
            end do

!-----------------------------------------------------------------------------!
!  Inverse FFT, step 3: 2-D complex-to-real transform
!-----------------------------------------------------------------------------!
!
!   C2R           (# = ny) A k_y        A y                 (# = ny) A y
!                          |            |                            |
!          +---+---+---+---+            +            +---+---+---+---+
!         /   /   /   /   /|           /|           /   /   /   /   /|
!        /   /   /   /    wrk      zy_sheet        /   /   /   /   wrk
!       /   /   /   /   /  |         /  |         /   /   /   /   /  |
!      +---+---+---+---+   |        +   |        +---+---+---+---+   |
!      |   |   |   |   |   |  C2R   |   |        |   |   |   |   |   |
!      |   |   |   |   |   | -----> |   | -----> |   |   |   |   |   |
!      |   |   |   |   |   |        |   |        |   |   |   |   |   |
!  <---|   |   |   |   |   +        |   +    <---|   |   |   |   |   +
!  x   |   |   |   |   |  /         |  /     x   |   |   |   |   |  /
!      |   |   |   |   | /          | /          |   |   |   |   | /
!      |   |   |   |   |/           |/           |   |   |   |   |/
!      +---+---+---+---+            +            +---+---+---+---+
!                     /            /                            /
!   (# = nz / 2 + 1) V k_z        V z                 (# = nz) V z
!
!        3   2   1   0  --- myid                   3   2   1   0  --- myid
!
!-----------------------------------------------------------------------------!

            do iz = 1, nx
                call DFFTW_EXECUTE(plan_c2r(iz, n))
                do iy = 1, ny
                    do ix = 1, nz
                    wrk(ix, iy, iz, n) = zy_sheet(ix, iy) ! unnormalized
                    end do
                end do
            end do

!-----------------------------------------------------------------------------!
            
        end if
        
        return
        
    end subroutine xFFT3d


    subroutine x_derivative(n, axis, nto)
        !=====================================================================!
        !    
        !   Compute derivative of wrk(:, :, :, n) with respect to axis and 
        !   store the result in wrk(:, :, :, nto). All in Fourier space.
        !   
        !=====================================================================!
    
    
    
    integer   :: ix, iy, iz, n, nto
    real(dp)    :: rtmp
    character :: axis
    
    ! In Fourier space, indices are ordered as (kz, kx, ky)
    
    ! First multiply by k:
    select case (axis)
    case('z')
        do iy = 1, nx
            do iz = 1, nx_all
                do ix = 1, nz + 2
                    wrk(ix, iz, iy, nto) = akz(ix) * wrk(ix, iz, iy, n)
                end do
            end do
        end do
    case('y')
        do iy = 1, nx
            wrk(:, :, iy, nto) = aky(iy) * wrk(:, :, iy, n)
        end do
    case('x')
        do iy = 1, nx
            do iz = 1, nx_all
                wrk(:, iz, iy, nto) = akx(iz) * wrk(:, iz, iy, n)
            end do
        end do
    case default 
        write (out, *) 'x_derivative: wrong value of axis: ', axis
        flush(out)
        stop        
    end select
    
    ! Multiply by i
    ! Real parts: wrk(ix, :, :, :), Imaginary  parts: wrk(ix + 1, :, :, :)
    ! for ix = 1, nz + 1, 2
    do iy = 1, nx
        do iz = 1, nx_all
            do ix = 1, nz + 1, 2
                
                rtmp = - wrk(ix + 1, iz, iy, nto)
                wrk(ix + 1, iz, iy, nto) = wrk(ix, iz, iy, nto)
                wrk(ix, iz, iy, nto) = rtmp
                
            end do
        end do
    end do
        
    end subroutine x_derivative

    subroutine x_dealias(i)
        
        integer, intent(in) :: i
        ! Dealias wrk(:, :, :, i)

        wrk(:, :, :, i) = dealias_matrix(:,:,:) * wrk(:, :, :, i)
    end subroutine x_dealias
    
end module x_fftw

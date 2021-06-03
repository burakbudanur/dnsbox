module m_numbers
    use iso_fortran_env, only : REAL32, REAL64, INT8, INT16, INT32, INT64

    ! Floating point numbers
    integer, parameter :: sp = REAL32
    integer, parameter :: dp = REAL64

    ! Integer numbers
    integer, parameter :: i1 = INT8
    integer, parameter :: i2 = INT16
    integer, parameter :: i4 = INT32
    integer, parameter :: i8 = INT64

    ! Corresponding formats

    ! (sign=1) + (leading=1) + (dot=1) + (eps=8)  + (exp=4) = 15
    character(*), parameter :: sp_f = 'es16.8'
    character(*), parameter :: sp_len = 'A16'
    ! (sign=1) + (leading=1) + (dot=1) + (eps=16)  + (exp=5) = 24
    character(*), parameter :: dp_f = 'es25.16'
    character(*), parameter :: dp_len = 'A25'

    ! (sign=1) + (digits=5) = 6
    character(*), parameter :: i2_f = 'i7'
    character(*), parameter :: i2_len = 'A7'
    ! (sign=1) + (digits=10) = 11
    character(*), parameter :: i4_f = 'i12'
    character(*), parameter :: i4_len = 'A12'
    ! (sign=1) + (digits=19) = 20
    character(*), parameter :: i8_f = 'i21'
    character(*), parameter :: i8_len = 'A21'

    ! Some constants

    real(dp), parameter :: PI     = 4.0_dp * DATAN(1.0_dp)
    real(dp), parameter :: TWO_PI = 2.0_dp * PI
    real(dp), parameter :: HALF   = 0.5_dp
    real(dp), parameter :: ZERO   = 0.0_dp
    real(dp), parameter :: ONE    = 1.0_dp
    real(dp), parameter :: TWO    = 2.0_dp
    real(dp), parameter :: THREE  = 3.0_dp
    real(dp), parameter :: FOUR   = 4.0_dp
    real(dp), parameter :: FIVE   = 5.0_dp
    real(dp), parameter :: SIX    = 6.0_dp

end module m_numbers

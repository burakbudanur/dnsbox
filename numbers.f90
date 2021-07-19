module numbers
    use, intrinsic :: iso_c_binding, only : C_FLOAT, C_DOUBLE, C_DOUBLE_COMPLEX, &
                                            C_INT32_T, C_INT64_T

    ! Floating point numbers
    integer, parameter :: sp  = C_FLOAT
    integer, parameter :: dp  = C_DOUBLE
    integer, parameter :: dpc = C_DOUBLE_COMPLEX

    ! Integer numbers
    integer, parameter :: i4 = C_INT32_T
    integer, parameter :: i8 = C_INT64_T

    ! Formats

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

    ! Constants
    complex(dpc), parameter :: imag_1 = (0, 1.0_dp)
    real(dp), parameter :: PI = 4.0_dp * ATAN(1.0_dp)

    ! https://numpy.org/doc/stable/reference/generated/numpy.finfo.html
    ! "The difference between 1.0 and the next smallest representable float
    !  larger than 1.0."
    real(dp), parameter :: epsilon = 2.221e-16_dp
    ! "The smallest positive floating point number with full precision."
    real(dp), parameter :: tiny = 2.226e-308_dp

    ! threshold below which floats will be considered as zero
    real(dp), parameter :: small = 1.222e-77_dp

    contains

    logical function are_equal(a,b)
        ! this should not be used if either a or b can be zero!
        real(dp), intent(in) :: a,b
        real(dp) :: reldiff

        are_equal = .false.
        reldiff = abs(a - b) / max(abs(a),abs(b))
        if (reldiff < epsilon) are_equal = .true.

    end function are_equal

end module numbers

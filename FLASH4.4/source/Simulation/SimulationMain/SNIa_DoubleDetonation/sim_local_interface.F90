module sim_local_interface

  use, intrinsic :: iso_fortran_env, only: dp=>real64

  implicit none

  real (dp), parameter :: zero = 0.0_dp
  real (dp), parameter :: one = 1.0_dp

  interface sim_interpolate1dWd
     subroutine sim_interpolate1dWd( volume, r_inner, r_outer, dens, temp, x)
     implicit none
     real, intent(in) :: volume, r_inner, r_outer
     real, intent(inout):: dens, temp, x(:)
     end subroutine
  end interface

  interface interp1d_linear
    module procedure interp1d_linear_vec
    module procedure interp1d_linear_scalar
  end interface interp1d_linear

contains

  function locate(x, n, xs)

    ! input variables
    integer, intent(in) :: n
    real(dp), intent(in) :: x, xs(n)

    ! function variable
    integer :: locate

    ! local variables
    integer :: ilo, ihi, imid

    if (x <= xs(1)) then
       locate = 1
    else if (x >= xs(n)) then
       locate = n+1
    else
       ilo = 1
       ihi = n
       do while (ilo+1 /= ihi)
          imid = (ilo+ihi)/2
          if (x <= xs(imid)) then
             ihi = imid
          else
             ilo = imid
          end if
       end do
       locate = ihi
    end if

  end function locate

  subroutine interp1d_linear_scalar( x_in, f_in, x_out, f_out )

    ! input variables
    real (dp), intent(in) :: x_in(:)
    real (dp), intent(in) :: f_in(:)
    real (dp), intent(in) :: x_out

    ! output variables
    real (dp), intent(out) :: f_out

    ! local variables
    integer :: i, ix, ix0, ix1, ix_max
    real (dp) :: dx
    real (dp) :: c(2), d(2)

    ix_max = size(x_in)
    ix = locate( x_out, ix_max, x_in ) - 1
    ix0 = min( max( ix  , 1 ), ix_max )
    ix1 = max( min( ix+1, ix_max ), 1 )
    if ( ix1 > ix0 ) then
      dx = ( x_out - x_in(ix0) ) / ( x_in(ix1) - x_in(ix0) )
    else
      dx = zero
    end if

    d(1) = one
    d(2) = dx

    c(1) = f_in(ix0)
    c(2) = f_in(ix1) - f_in(ix0)

    f_out = sum( d*c )

    return
  end subroutine interp1d_linear_scalar

  subroutine interp1d_linear_vec( x_in, f_in, x_out, f_out )

    ! input variables
    real (dp), intent(in) :: x_in(:)
    real (dp), intent(in) :: f_in(:)
    real (dp), intent(in) :: x_out(:)

    ! output variables
    real (dp), intent(out) :: f_out(size(x_out))

    ! local variables
    integer :: i, ix, ix0, ix1, ix_max
    real (dp) :: dx
    real (dp) :: c(2), d(2)

    ix_max = size(x_in)

    do i = 1, size(x_out,1)

      ix = locate( x_out(i), ix_max, x_in ) - 1
      ix0 = min( max( ix  , 1 ), ix_max )
      ix1 = max( min( ix+1, ix_max ), 1 )
      if ( ix1 > ix0 ) then
        dx = ( x_out(i) - x_in(ix0) ) / ( x_in(ix1) - x_in(ix0) )
      else
        dx = zero
      end if

      d(1) = one
      d(2) = dx

      c(1) = f_in(ix0)
      c(2) = f_in(ix1) - f_in(ix0)

      f_out(i) = sum( d*c )

    end do

    return
  end subroutine interp1d_linear_vec

end module sim_local_interface


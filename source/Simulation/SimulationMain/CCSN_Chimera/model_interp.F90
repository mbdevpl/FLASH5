module model_interp_module

  use, intrinsic :: iso_fortran_env, only: dp=>real64

  implicit none

  real (dp), parameter :: zero = 0.0_dp
  real (dp), parameter :: one = 1.0_dp

  interface interp1d_linear
    module procedure interp1d_linear_vec
    module procedure interp1d_linear_scalar
  end interface interp1d_linear

  interface interp2d_linear
    module procedure interp2d_linear_vec
    module procedure interp2d_linear_scalar
  end interface interp2d_linear

  interface interp3d_linear
    module procedure interp3d_linear_vec
    module procedure interp3d_linear_scalar
  end interface interp3d_linear

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
       locate = 0
    else if (x >= xs(n)) then
       locate = n
    else
       ilo = 1
       ihi = n-1
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

  function locate_descending(x, n, xs)

    ! input variables
    integer, intent(in) :: n
    real(dp), intent(in) :: x, xs(n)

    ! function variable
    integer :: locate_descending

    ! local variables
    integer :: ilo, ihi, imid

    if (x >= xs(1)) then
       locate_descending = 0
    else if (x <= xs(n)) then
       locate_descending = n
    else
       ilo = 1
       ihi = n-1
       do while (ilo+1 /= ihi)
          imid = (ilo+ihi)/2
          if (x >= xs(imid)) then
             ihi = imid
          else
             ilo = imid
          end if
       end do
       locate_descending = ihi
    end if

  end function locate_descending

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

  subroutine interp2d_linear_scalar( x_in, y_in, f_in, x_out, y_out, f_out )

    ! input variables
    real (dp), intent(in) :: x_in(:)
    real (dp), intent(in) :: y_in(:)
    real (dp), intent(in) :: f_in(:,:)
    real (dp), intent(in) :: x_out
    real (dp), intent(in) :: y_out

    ! output variables
    real (dp), intent(out) :: f_out

    ! local variables
    integer :: i, ix, ix0, ix1, ix_max
    integer :: j, iy, iy0, iy1, iy_max
    real (dp) :: dx, dy
    real (dp) :: c(4), d(4)

    ix_max = size(x_in)
    iy_max = size(y_in)

    ix = locate( x_out, ix_max, x_in ) - 1
    ix0 = min( max( ix  , 1 ), ix_max )
    ix1 = max( min( ix+1, ix_max ), 1 )
    if ( ix1 > ix0 ) then
      dx = ( x_out - x_in(ix0) ) / ( x_in(ix1) - x_in(ix0) )
    else
      dx = zero
    end if

    iy = locate( y_out, iy_max, y_in ) - 1
    iy0 = min( max( iy  , 1 ), iy_max )
    iy1 = max( min( iy+1, iy_max ), 1 )
    if ( iy1 > iy0 ) then
      dy = ( y_out - y_in(iy0) ) / ( y_in(iy1) - y_in(iy0) )
    else
      dy = zero
    end if

    d(1) = one
    d(2) = dx
    d(3) = dy
    d(4) = dx*dy

    c(1) = f_in(ix0,iy0)
    c(2) = f_in(ix1,iy0) - f_in(ix0,iy0)
    c(3) = f_in(ix0,iy1) - f_in(ix0,iy0)
    c(4) = f_in(ix1,iy1) - f_in(ix0,iy1) + &
           f_in(ix0,iy0) - f_in(ix1,iy0)

    f_out = sum( d*c )

    return
  end subroutine interp2d_linear_scalar

  subroutine interp2d_linear_vec( x_in, y_in, f_in, x_out, y_out, f_out )

    ! input variables
    real (dp), intent(in) :: x_in(:)
    real (dp), intent(in) :: y_in(:)
    real (dp), intent(in) :: f_in(:,:)
    real (dp), intent(in) :: x_out(:,:)
    real (dp), intent(in) :: y_out(:,:)

    ! output variables
    real (dp), intent(out) :: f_out(size(x_out,1),size(x_out,2))

    ! local variables
    integer :: i, ix, ix0, ix1, ix_max
    integer :: j, iy, iy0, iy1, iy_max
    real (dp) :: dx, dy
    real (dp) :: c(4), d(4)

    ix_max = size(x_in)
    iy_max = size(y_in)

    do j = 1, size(x_out,2)
      do i = 1, size(x_out,1)

        ix = locate( x_out(i,j), ix_max, x_in ) - 1
        ix0 = min( max( ix  , 1 ), ix_max )
        ix1 = max( min( ix+1, ix_max ), 1 )
        if ( ix1 > ix0 ) then
          dx = ( x_out(i,j) - x_in(ix0) ) / ( x_in(ix1) - x_in(ix0) )
        else
          dx = zero
        end if

        iy = locate( y_out(i,j), iy_max, y_in ) - 1
        iy0 = min( max( iy  , 1 ), iy_max )
        iy1 = max( min( iy+1, iy_max ), 1 )
        if ( iy1 > iy0 ) then
          dy = ( y_out(i,j) - y_in(iy0) ) / ( y_in(iy1) - y_in(iy0) )
        else
          dy = zero
        end if

        d(1) = one
        d(2) = dx
        d(3) = dy
        d(4) = dx*dy

        c(1) = f_in(ix0,iy0)
        c(2) = f_in(ix1,iy0) - f_in(ix0,iy0)
        c(3) = f_in(ix0,iy1) - f_in(ix0,iy0)
        c(4) = f_in(ix1,iy1) - f_in(ix0,iy1) + &
               f_in(ix0,iy0) - f_in(ix1,iy0)

        f_out(i,j) = sum( d*c )

      end do
    end do

    return
  end subroutine interp2d_linear_vec

  subroutine interp3d_linear_scalar( x_in, y_in, z_in, f_in, x_out, y_out, z_out, f_out )

    ! input variables
    real (dp), intent(in) :: x_in(:)
    real (dp), intent(in) :: y_in(:)
    real (dp), intent(in) :: z_in(:)
    real (dp), intent(in) :: f_in(:,:,:)
    real (dp), intent(in) :: x_out
    real (dp), intent(in) :: y_out
    real (dp), intent(in) :: z_out

    ! output variables
    real (dp), intent(out) :: f_out

    ! local variables
    integer :: i, ix, ix0, ix1, ix_max
    integer :: j, iy, iy0, iy1, iy_max
    integer :: k, iz, iz0, iz1, iz_max
    real (dp) :: dx, dy, dz
    real (dp) :: c(8), d(8)

    ix_max = size(x_in)
    iy_max = size(y_in)
    iz_max = size(z_in)

    ix = locate( x_out, ix_max, x_in ) - 1
    ix0 = min( max( ix  , 1 ), ix_max )
    ix1 = max( min( ix+1, ix_max ), 1 )
    if ( ix1 > ix0 ) then
      dx = ( x_out - x_in(ix0) ) / ( x_in(ix1) - x_in(ix0) )
    else
      dx = zero
    end if

    iy = locate( y_out, iy_max, y_in ) - 1
    iy0 = min( max( iy  , 1 ), iy_max )
    iy1 = max( min( iy+1, iy_max ), 1 )
    if ( iy1 > iy0 ) then
      dy = ( y_out - y_in(iy0) ) / ( y_in(iy1) - y_in(iy0) )
    else
      dy = zero
    end if

    iz = locate( z_out, iz_max, z_in ) - 1
    iz0 = min( max( iz  , 1 ), iz_max )
    iz1 = max( min( iz+1, iz_max ), 1 )
    if ( iz1 > iz0 ) then
      dz = ( z_out - z_in(iz0) ) / ( z_in(iz1) - z_in(iz0) )
    else
      dz = zero
    end if

    d(1) = one
    d(2) = dx
    d(3) = dy
    d(4) = dz
    d(5) = dx*dy
    d(6) = dx*dz
    d(7) = dy*dz
    d(8) = dx*dy*dz

    c(1) = f_in(ix0,iy0,iz0)
    c(2) = f_in(ix1,iy0,iz0) - f_in(ix0,iy0,iz0)
    c(3) = f_in(ix0,iy1,iz0) - f_in(ix0,iy0,iz0)
    c(4) = f_in(ix0,iy0,iz1) - f_in(ix0,iy0,iz0)
    c(5) = f_in(ix1,iy1,iz0) - f_in(ix0,iy1,iz0) + &
           f_in(ix0,iy0,iz0) - f_in(ix1,iy0,iz0)
    c(6) = f_in(ix1,iy0,iz1) - f_in(ix0,iy0,iz1) + &
           f_in(ix0,iy0,iz0) - f_in(ix1,iy0,iz0)
    c(7) = f_in(ix0,iy1,iz1) - f_in(ix0,iy0,iz1) + &
           f_in(ix0,iy0,iz0) - f_in(ix0,iy1,iz0)
    c(8) = f_in(ix1,iy1,iz1) - f_in(ix0,iy1,iz1) + &
           f_in(ix0,iy0,iz1) - f_in(ix1,iy0,iz1) + &
           f_in(ix0,iy1,iz0) - f_in(ix1,iy1,iz0) + &
           f_in(ix1,iy0,iz0) - f_in(ix0,iy0,iz0)

    f_out = sum( d*c )

    return
  end subroutine interp3d_linear_scalar

  subroutine interp3d_linear_vec( x_in, y_in, z_in, f_in, x_out, y_out, z_out, f_out )

    ! input variables
    real (dp), intent(in) :: x_in(:)
    real (dp), intent(in) :: y_in(:)
    real (dp), intent(in) :: z_in(:)
    real (dp), intent(in) :: f_in(:,:,:)
    real (dp), intent(in) :: x_out(:,:,:)
    real (dp), intent(in) :: y_out(:,:,:)
    real (dp), intent(in) :: z_out(:,:,:)

    ! output variables
    real (dp), intent(out) :: f_out(size(x_out,1),size(x_out,2),size(x_out,3))

    ! local variables
    integer :: i, ix, ix0, ix1, ix_max
    integer :: j, iy, iy0, iy1, iy_max
    integer :: k, iz, iz0, iz1, iz_max
    real (dp) :: dx, dy, dz
    real (dp) :: c(8), d(8)

    ix_max = size(x_in)
    iy_max = size(y_in)
    iz_max = size(z_in)

    do k = 1, size(x_out,3)
      do j = 1, size(x_out,2)
        do i = 1, size(x_out,1)

          ix = locate( x_out(i,j,k), ix_max, x_in ) - 1
          ix0 = min( max( ix  , 1 ), ix_max )
          ix1 = max( min( ix+1, ix_max ), 1 )
          if ( ix1 > ix0 ) then
            dx = ( x_out(i,j,k) - x_in(ix0) ) / ( x_in(ix1) - x_in(ix0) )
          else
            dx = zero
          end if

          iy = locate( y_out(i,j,k), iy_max, y_in ) - 1
          iy0 = min( max( iy  , 1 ), iy_max )
          iy1 = max( min( iy+1, iy_max ), 1 )
          if ( iy1 > iy0 ) then
            dy = ( y_out(i,j,k) - y_in(iy0) ) / ( y_in(iy1) - y_in(iy0) )
          else
            dy = zero
          end if

          iz = locate( z_out(i,j,k), iz_max, z_in ) - 1
          iz0 = min( max( iz  , 1 ), iz_max )
          iz1 = max( min( iz+1, iz_max ), 1 )
          if ( iz1 > iz0 ) then
            dz = ( z_out(i,j,k) - z_in(iz0) ) / ( z_in(iz1) - z_in(iz0) )
          else
            dz = zero
          end if

          d(1) = one
          d(2) = dx
          d(3) = dy
          d(4) = dz
          d(5) = dx*dy
          d(6) = dx*dz
          d(7) = dy*dz
          d(8) = dx*dy*dz

          c(1) = f_in(ix0,iy0,iz0)
          c(2) = f_in(ix1,iy0,iz0) - f_in(ix0,iy0,iz0)
          c(3) = f_in(ix0,iy1,iz0) - f_in(ix0,iy0,iz0)
          c(4) = f_in(ix0,iy0,iz1) - f_in(ix0,iy0,iz0)
          c(5) = f_in(ix1,iy1,iz0) - f_in(ix0,iy1,iz0) + &
                 f_in(ix0,iy0,iz0) - f_in(ix1,iy0,iz0)
          c(6) = f_in(ix1,iy0,iz1) - f_in(ix0,iy0,iz1) + &
                 f_in(ix0,iy0,iz0) - f_in(ix1,iy0,iz0)
          c(7) = f_in(ix0,iy1,iz1) - f_in(ix0,iy0,iz1) + &
                 f_in(ix0,iy0,iz0) - f_in(ix0,iy1,iz0)
          c(8) = f_in(ix1,iy1,iz1) - f_in(ix0,iy1,iz1) + &
                 f_in(ix0,iy0,iz1) - f_in(ix1,iy0,iz1) + &
                 f_in(ix0,iy1,iz0) - f_in(ix1,iy1,iz0) + &
                 f_in(ix1,iy0,iz0) - f_in(ix0,iy0,iz0)

          f_out(i,j,k) = sum( d*c )

        end do
      end do
    end do

    return
  end subroutine interp3d_linear_vec

end module model_interp_module

! 
!! Routines for Euler PC integration.


! Order 2: Adams-Bashforth (explicit)
subroutine sm_pc_AB2(n,X0,DX0,DXm1,dt,m,X1)
  implicit none

  integer, intent(in) :: n, m
  real, intent(in) :: dt(-m:0)
  real, dimension(n), intent(in) :: X0,DX0,DXm1
  real, dimension(n), intent(out):: X1

  ! integ parameters
  real, dimension(2) :: alpha

  ! Compute the alpha parameters
  alpha = (/ 1 + dt(0)/(2.*dt(-1)), &
             -dt(0)/(2.*dt(-1)) /)

  ! Integrate
  !X1(1:n) = X0(1:n) + dt*0.5*( 3*DX0(1:n) - DXm1(1:n) ) 
  X1(1:n) = X0(1:n) + dt(0)*( alpha(1)*DX0(1:n) + alpha(2)*DXm1(1:n) ) 

  return
end subroutine sm_pc_AB2

! Order 2: Adams-Moulton (implicit correction)
subroutine sm_pc_AM2(n,X0,DX0,DXm1,DX1,dt,m,X1)
  implicit none

  integer, intent(in) :: n,m
  real, intent(in) :: dt(-m:0)
  real, dimension(n), intent(in) :: DXm1,X0,DX0,DX1
  real, dimension(n), intent(out):: X1
  ! integ parameters
  real, dimension(3) :: alpha

  alpha = (/(2 + dt(-1)/(dt(0) + dt(-1)))/6.,          &
            (3 + dt(0)/dt(-1))/6.,                     &
            -dt(0)**2/(6.*dt(-1)*(dt(0) + dt(-1))) /)

  !X1(1:n) = X0(1:n) + dt*( 5./12.*DX1(1:n) + 2./3.*DX0(1:n) - 1./12.*DXm1(1:n) )
  X1(1:n) = X0(1:n) + dt(0)*( alpha(1)*DX1(1:n) + alpha(2)*DX0(1:n) + alpha(3)*DXm1(1:n) )

  return
end subroutine sm_pc_AM2

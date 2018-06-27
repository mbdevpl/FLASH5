! 
!! Routines for Euler PC integration.


! Order 3: Adams-Bashforth (explicit)
subroutine sm_pc_AB3(n,X0,DX0,DXm1,DXm2,dt,m,X1)
  implicit none

  integer, intent(in) :: n,m
  real, intent(in) :: dt(-m:0)
  real, dimension(n), intent(in) :: X0,DX0,DXm1,DXm2
  real, dimension(n), intent(out):: X1
  ! integ parameters
  real, dimension(3) :: alpha

  alpha = (/ (2*dt(0)**2 + 6*dt(-1)*(dt(-1) + dt(-2)) + &
          3*dt(0)*(2*dt(-1) + dt(-2)))/                 &
          (6.*dt(-1)*(dt(-1) + dt(-2))),                &
          -(dt(0)*(2*dt(0) + 3*(dt(-1) + dt(-2))))/     &
          (6.*dt(-1)*dt(-2)),                           &
          (dt(0)*(2*dt(0) + 3*dt(-1)))/                 &
          (6.*dt(-2)*(dt(-1) + dt(-2))) /)

  !X1(1:n) = X0(1:n) + dt*( 23./12.*DX0(1:n) - 4./3.*DXm1(1:n) + 5./12.*DXm2(1:n) ) 
  X1(1:n) = X0(1:n) + dt(0)*( alpha(1)*DX0(1:n) + alpha(2)*DXm1(1:n) + alpha(3)*DXm2(1:n) ) 

  return
end subroutine sm_pc_AB3

! Order 3: Adams-Moulton (implicit correction)
subroutine sm_pc_AM3(n,X0,DX0,DXm1,DXm2,DX1,dt,m,X1)
  implicit none

  integer, intent(in) :: n, m
  real, intent(in) :: dt(-m:0)
  real, dimension(n), intent(in) :: DXm1,DXm2,X0,DX0,DX1
  real, dimension(n), intent(out):: X1
  ! integ parameters
  real, dimension(4) :: alpha

  alpha = (/ (3*dt(0)**2 +                              &
             6*dt(-1)*(dt(-1) + dt(-2)) +               &
             4*dt(0)*(2*dt(-1) + dt(-2)))/              &
             (12.*(dt(0) + dt(-1))*                     &
             (dt(0) + dt(-1) + dt(-2))),                &
             (dt(0)**2 + 6*dt(-1)*(dt(-1) + dt(-2)) +   &
             2*dt(0)*(2*dt(-1) + dt(-2)))/              &
             (12.*dt(-1)*(dt(-1) + dt(-2))),            &
             -(dt(0)**2*(dt(0) + 2*(dt(-1) + dt(-2))))/ &
             (12.*dt(-1)*(dt(0) + dt(-1))*dt(-2)),      &
             (dt(0)**2*(dt(0) + 2*dt(-1)))/             &
             (12.*dt(-2)*(dt(-1) + dt(-2))*             &
             (dt(0) + dt(-1) + dt(-2))) /)

  !X1(1:n) = X0(1:n) + dt*( 3./8.*DX1(1:n) + 19./24.*DX0(1:n) &
  !                        - 5./24.*DXm1(1:n) + 1./24.*DXm2(1:n) )
  X1(1:n) = X0(1:n) + dt(0)*( alpha(1)*DX1(1:n)  + alpha(2)*DX0(1:n) &
                             +alpha(3)*DXm1(1:n) + alpha(4)*DXm2(1:n) )


  return
end subroutine sm_pc_AM3

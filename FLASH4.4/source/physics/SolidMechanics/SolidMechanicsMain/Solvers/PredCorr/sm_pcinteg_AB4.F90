! 
!! Routines for Euler PC integration.


! Order 4: Adams-Bashforth (explicit)
subroutine sm_pc_AB4(n,X0,DX0,DXm1,DXm2,DXm3,dt,m,X1)
  implicit none

  integer, intent(in) :: n, m
  real, intent(in) :: dt(-m:0)
  real, dimension(n), intent(in) :: X0,DX0,DXm1,DXm2,DXm3
  real, dimension(n), intent(out):: X1
  ! integ parameters
  real, dimension(4) :: alpha

  ! build alpha's (generated in Mathematica)
  alpha = (/ (3*dt(0)**3 +                                &
       12*dt(-1)*(dt(-1) + dt(-2))*                       &
       (dt(-1) + dt(-2) + dt(-3)) +                       &
       4*dt(0)**2*(3*dt(-1) + 2*dt(-2) + dt(-3)) +        &
       6*dt(0)*((dt(-1) + dt(-2))*(3*dt(-1) + dt(-2)) +   &
       (2*dt(-1) + dt(-2))*dt(-3)))/                      &
       (12.*dt(-1)*(dt(-1) + dt(-2))*                     &
       (dt(-1) + dt(-2) + dt(-3))),                       &
       (dt(0)*(-3*dt(0)**2 -                              &
       6*(dt(-1) + dt(-2))*                               &
       (dt(-1) + dt(-2) + dt(-3)) -                       &
       4*dt(0)*(2*(dt(-1) + dt(-2)) + dt(-3))))/          &
       (12.*dt(-1)*dt(-2)*(dt(-2) + dt(-3))),             &
       (dt(0)*(3*dt(0)**2 +                               &
       6*dt(-1)*(dt(-1) + dt(-2) + dt(-3)) +              &
       4*dt(0)*(2*dt(-1) + dt(-2) + dt(-3))))/            &
       (12.*dt(-2)*(dt(-1) + dt(-2))*dt(-3)),             &
       (dt(0)*(-3*dt(0)**2 - 6*dt(-1)*(dt(-1) + dt(-2)) - &
       4*dt(0)*(2*dt(-1) + dt(-2))))/                     &
       (12.*dt(-3)*(dt(-2) + dt(-3))*                     &
       (dt(-1) + dt(-2) + dt(-3))) /)

  !X1(1:n) = X0(1:n) + dt*(   55./24.*DX0(1:n)  - 59./24.*DXm1(1:n)  &
  !                         + 37./24.*DXm2(1:n) -   3./8.*DXm3(1:n) ) 
  X1(1:n) = X0(1:n) + dt(0)*( alpha(1)*DX0(1:n)  + alpha(2)*DXm1(1:n)  &
                             +alpha(3)*DXm2(1:n) + alpha(4)*DXm3(1:n) )

  return
end subroutine sm_pc_AB4

! Order 4: Adams-Moulton (implicit correction)
subroutine sm_pc_AM4(n,X0,DX0,DXm1,DXm2,DXm3,DX1,dt,m,X1)
  implicit none

  integer, intent(in) :: n, m
  real, intent(in) :: dt(-m:0)
  real, dimension(n), intent(in) :: X0,DX0,DX1,DXm1,DXm2,DXm3
  real, dimension(n), intent(out):: X1
  ! integ parameters
  real, dimension(5) :: alpha

  ! build alpha's (generated in Mathematica)
  alpha = (/ (12*dt(0)**3 +                        &
       30*dt(-1)*(dt(-1) + dt(-2))*                &
       (dt(-1) + dt(-2) + dt(-3)) +                &
       15*dt(0)**2*                                &
       (3*dt(-1) + 2*dt(-2) + dt(-3)) +            &
       20*dt(0)*((dt(-1) + dt(-2))*                &
       (3*dt(-1) + dt(-2)) +                       &
       (2*dt(-1) + dt(-2))*dt(-3)))/               &
       (60.*(dt(0) + dt(-1))*                      &
       (dt(0) + dt(-1) + dt(-2))*                  &
       (dt(0) + dt(-1) + dt(-2) + dt(-3))),        &
       (3*dt(0)**3 + 30*dt(-1)*(dt(-1) + dt(-2))*  &
       (dt(-1) + dt(-2) + dt(-3)) +                &
       5*dt(0)**2*                                 &
       (3*dt(-1) + 2*dt(-2) + dt(-3)) +            &
       10*dt(0)*((dt(-1) + dt(-2))*                &
       (3*dt(-1) + dt(-2)) +                       &
       (2*dt(-1) + dt(-2))*dt(-3)))/               &
       (60.*dt(-1)*(dt(-1) + dt(-2))*              &
       (dt(-1) + dt(-2) + dt(-3))),                &
       (dt(0)**2*(-3*dt(0)**2 -                    &
       10*(dt(-1) + dt(-2))*                       &
       (dt(-1) + dt(-2) + dt(-3)) -                &
       5*dt(0)*(2*(dt(-1) + dt(-2)) + dt(-3)))     &
       )/                                          &
       (60.*dt(-1)*(dt(0) + dt(-1))*dt(-2)*        &
       (dt(-2) + dt(-3))),                         &
       (dt(0)**2*(3*dt(0)**2 +                     &
       10*dt(-1)*(dt(-1) + dt(-2) + dt(-3)) +      &
       5*dt(0)*(2*dt(-1) + dt(-2) + dt(-3))))/     &
       (60.*dt(-2)*(dt(-1) + dt(-2))*              &
       (dt(0) + dt(-1) + dt(-2))*dt(-3)),          &
       (dt(0)**2*(-3*dt(0)**2 -                    &
       10*dt(-1)*(dt(-1) + dt(-2)) -               &
       5*dt(0)*(2*dt(-1) + dt(-2))))/              &
       (60.*dt(-3)*(dt(-2) + dt(-3))*              &
       (dt(-1) + dt(-2) + dt(-3))*                 &
       (dt(0) + dt(-1) + dt(-2) + dt(-3))) /)

  !X1(1:n) = X0(1:n) + dt/720.*( 251.*DX1(1:n)  + 646.*DX0(1:n)  &
  !                             -264.*DXm1(1:n) + 106.*DXm2(1:n) &
  !                             - 19.*DXm3(1:n) )

  X1(1:n) = X0(1:n) + dt(0)*( alpha(1)*DX1(1:n)  + alpha(2)*DX0(1:n)  &
                             +alpha(3)*DXm1(1:n) + alpha(4)*DXm2(1:n) &
                             +alpha(5)*DXm3(1:n) )

  return
end subroutine sm_pc_AM4

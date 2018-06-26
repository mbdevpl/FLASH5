! 
!! Routines for Euler PC integration.


! Explicit Euler step:
subroutine sm_pcem(n,X0,DX0,dt,X1)
  implicit none

  integer, intent(in) :: n
  real, intent(in) :: dt
  real, dimension(n), intent(in) :: X0,DX0
  real, dimension(n), intent(out):: X1

  X1(1:n) = X0(1:n) + dt*DX0(1:n)

  return
end subroutine sm_pcem

! Modified Euler Correction:
subroutine sm_pcmem(n,X0,DX0,DX1,dt,X1)
  implicit none

  integer, intent(in) :: n
  real, intent(in) :: dt
  real, dimension(n), intent(in) :: X0,DX0,DX1
  real, dimension(n), intent(out):: X1

  real :: halfdt

  halfdt = 0.5*dt
 
  X1(1:n) = X0(1:n) + halfdt*(DX0(1:n)+DX1(1:n))

  return
end subroutine sm_pcmem

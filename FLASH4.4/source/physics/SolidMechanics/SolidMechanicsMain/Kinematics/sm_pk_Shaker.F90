!     
! File:   sm_pk_Shaker.F90
! Author: tim
!
! 
#include "Flash.h"
#include "constants.h"

subroutine sm_pk_Shaker(time,n,Xe,v,vd,vdd,parameters)
  ! This subroutine computes the prescribed coordinates {v, vd, vdd}
  ! for a shaker 
  !
  ! The inputs:
  !     n: number of points to computer
  !    Xe: (3,n) grid nodes (ref configuration) to move
  !     v: (3,n) global displacement of the nodes under the kinematic transformation (not position)
  !    vd: (3,n) global velocities of the prescribed nodes
  !   vdd: (3,n) global accelerations of the prescribed nodes
  !   
  
  use sm_pk_interface, only: sm_pk_Shaker_harmonic
  implicit none
    
  !
  ! IO Variables
  !
  real, intent(in)    :: time
  integer, intent(in) :: n
  real, intent(in),  dimension(NDIM,n) :: Xe
  real, intent(in),  dimension(:)      :: parameters
  real, intent(out), dimension(NDIM,n) :: v, vd, vdd
  
  !
  ! Local Variables
  !
  real :: Ax, Ay, Az, omega_x, omega_y, omega_z, phi_x, phi_y, phi_z, tau
  real :: qx, qdx, qddx, qy, qdy, qddy, qz, qdz, qddz
  integer :: j

  !
  ! Get the parameters
  !
  tau = parameters(1)
  ! x:
  Ax      = parameters(2)
  omega_x = parameters(3)
  phi_x   = parameters(4)
  ! y:
  Ay      = parameters(5)
  omega_y = parameters(6)
  phi_y   = parameters(7)
  ! z:
#if NDIM == 3
  Az      = parameters(8)
  omega_z = parameters(9)
  phi_z   = parameters(10)
#endif  

  !
  ! Compute the displacements for each direction
  !
  call sm_pk_Shaker_harmonic(time, omega_x, phi_x, Ax, tau, qx, qdx, qddx)
  call sm_pk_Shaker_harmonic(time, omega_y, phi_y, Ay, tau, qy, qdy, qddy)
#if NDIM == 3
  call sm_pk_Shaker_harmonic(time, omega_z, phi_z, Az, tau, qz, qdz, qddz)
#endif
    
  !
  ! loop over each node of RestSurface and set (v, vd, vdd) -> (qn, qdn, qddn)
#if NDIM == 2    
  do j = 1,n    
     V(1:NDIM,j)   = (/   qx,   qy /)
     Vd(1:NDIM,j)  = (/  qdx,  qdy /)
     Vdd(1:NDIM,j) = (/ qddx, qddy /)
  enddo
#elif NDIM == 3
  do j = 1,n    
     V(1:NDIM,j)   = (/   qx,   qy,   qz /)
     Vd(1:NDIM,j)  = (/  qdx,  qdy,  qdz /)
     Vdd(1:NDIM,j) = (/ qddx, qddy, qddz /)
  enddo
#endif
    
end subroutine sm_pk_Shaker

subroutine sm_pk_Shaker_harmonic(time, omega, phi, Amp, tau, qe, qde, qdde)

  implicit none
  
  ! IO
  real, intent(in) :: time, omega, phi, Amp, tau
  real, intent(out) :: qe, qde, qdde

  ! Internal
  real :: damp, s, c

  ! Damped start
  damp = exp(-time/tau)
  s = Amp*sin(omega*time + phi)
  c = Amp*cos(omega*time + phi)

  ! Compute Displacement
  qe = (1.-damp)*s

  ! Compute the velocity
  qde = damp*s/tau + (1.-damp)*omega*c

  ! Compute the acceleration
  qdde = -damp*s/tau**2 + 2*c*damp*omega/tau - (1. - damp)*s*omega**2

  return

end subroutine sm_pk_Shaker_harmonic

!     
! File:   sm_pk_Flapping_BermanWang.F90
! Author: tim
!
! routine is based on the kine
#include "Flash.h"
#include "constants.h"

subroutine sm_pk_Flapping_BermanWang(time,n,Xe,v,vd,vdd,parameters)
  ! This subroutine computes the prescribed coordinates {v, vd, vdd} for a wing.
  ! It assumes that the base wing is more or less defined in the x-y plane with y being forward, and x is the span of the wing.
  ! The inputs:
  !     n: number of points to computer
  !    Xe: (3,n) grid nodes (ref configuration) to move
  !     v: (3,n) global displacement of the nodes under the kinematic transformation (not position)
  !    vd: (3,n) global velocities of the prescribed nodes
  !   vdd: (3,n) global accelerations of the prescribed nodes
  !   

  use sm_pk_interface, only: sm_pk_flapping_BermanWangFunctions,sm_pk_flapping_BermanWang_RotMat 
  use Driver_interface, only :  Driver_abortFlash

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
  real, dimension(NDIM,NDIM) :: R, Rd, Rdd
  real, dimension(NDIM)      :: phi, theta, eta, q, qd, qdd
  real    :: tau, damp, damp_d, damp_dd
  integer :: a,i,j,idx

  v  =0.0
  vd =0.0
  vdd=0.0

#if NDIM == MDIM
  ! set the damping parameter
  tau = parameters(12)
  if( tau < 1.e-12 ) then
     damp    = 0.
     damp_d  = 0.
     damp_dd = 0.
  else
     damp    = exp(-time/tau)
     damp_d  = -damp/tau
     damp_dd = -damp_d/tau
  end if
    
  ! get current angles:
  call sm_pk_flapping_BermanWangFunctions(time, phi, theta, eta, parameters)
    
  ! get rotation matrix
  call sm_pk_flapping_BermanWang_RotMat(phi, theta, eta, R, Rd, Rdd)
    
  !
  ! loop over each node of RestSurface and set (v, vd, vdd) -> (qn, qdn, qddn)
  !
  ! Recall that x+q = R*x, so:
  !      q   = R*x - x = (R-I)*x
  !      qd  = Rd*x
  !      qdd = Rdd*x    
  !
  !      With a damped start, this becomes:
  !      Q   = (1-damp)*q
  !      Qd  = -damp_d*q + (1-damp)*qd
  !      Qdd = -damp_dd*qdd - 2*damp_d*qd + (1-damp)*qdd
  do j = 1,n
        
     q   = matmul(R  , Xe(1:NDIM,j) ) - Xe(1:NDIM,j)
     qd  = matmul(Rd , Xe(1:NDIM,j) )
     qdd = matmul(Rdd, Xe(1:NDIM,j) )

     V(1:NDIM,j)   = (1.-damp)*q
     Vd(1:NDIM,j)  = -damp_d*q  + (1.-damp)*qd
     Vdd(1:NDIM,j) = -damp_dd*q - 2.*damp_d*qd + (1.-damp)*qdd

  enddo

#else

  call Driver_abortFlash("sm_pk_Flapping_BermanWang: Error, Funtion not defined for 2D")

#endif

  return
    
end subroutine sm_pk_Flapping_BermanWang
    
subroutine sm_pk_flapping_BermanWangFunctions(t, phi, theta, eta, parameters)
  ! This subroutine computes the 3-2-1 (ish) angles defined the Berman and Wang 2007 JFM paper.

  implicit none
    
  !
  ! IO variables
  !
  real, intent(in)                :: t
  real, intent(out), dimension(MDIM) :: phi, theta, eta
  real, intent(in), dimension(:)  :: parameters

  !
  ! local variables
  !
  real :: omegaf, phi_m, Cphi, &
          theta_m, theta_0, PHI_theta, N, &
          eta_m, eta_0, PHI_eta, Ceta, tau

  ! set parameters:
  omegaf    = parameters(1)
  phi_m     = parameters(2)
  Cphi      = parameters(3)
  ! theta:
  theta_m   = parameters(4)
  theta_0   = parameters(5)
  PHI_theta = parameters(6)
  N         = parameters(7)
  ! eta:
  eta_m     = parameters(8)
  eta_0     = parameters(9)
  PHI_eta   = parameters(10)
  Ceta      = parameters(11)
  ! damped start (not currently used inside R)
  tau       = parameters(12)

  ! compute phi:
  phi(1) = phi_m/asin(Cphi)*asin(Cphi*sin(omegaf*t))
  phi(2) = phi_m*Cphi*omegaf*cos(omegaf*t)/asin(Cphi)/sqrt(1.0 - (Cphi*sin(omegaf*t))**2)
  phi(3) = 2.0*sqrt(2.0)*phi_m*Cphi*(Cphi**2-1.0)*omegaf**2*sin(omegaf*t)/asin(Cphi) &
       /(2.d0-Cphi**2 + Cphi**2 * cos(2.d0*omegaf*t))**1.5

  ! compute theta
  theta(1) = theta_0 + theta_m*cos(N*omegaf*t + PHI_theta)
  theta(2) = -theta_m*N*omegaf*sin(N*omegaf*t + PHI_theta)
  theta(3) = -theta_m*N**2*omegaf**2*cos(N*omegaf*t + PHI_theta)

  ! compute eta
  eta(1) = eta_0 + eta_m/tanh(Ceta)*tanh(Ceta*sin(omegaf*t + PHI_eta))
  eta(2) = Ceta*eta_m*omegaf*cos(PHI_eta + omegaf*t)/tanh(Ceta)/cosh(Ceta*sin(PHI_eta + omegaf*t))**2
  eta(3) = -Ceta*eta_m*omegaf**2/tanh(Ceta)/cosh(Ceta*sin(PHI_eta + omegaf*t))**2 &
       *( sin(PHI_eta + omegaf*t) + 2*Ceta*cos(PHI_eta + omegaf*t)**2*tanh(Ceta*sin(PHI_eta + omegaf*t)))
              
end subroutine sm_pk_flapping_BermanWangFunctions

subroutine sm_pk_flapping_BermanWang_RotMat(phi, theta, eta, R, Rd, Rdd)
  implicit none

  !
  ! IO variables
  !
  real, intent(in),  dimension(MDIM)   :: phi, theta, eta
  real, intent(out), dimension(MDIM,MDIM) :: R, Rd, Rdd

  !
  ! Internal Variables
  !
  real :: c1,s1, c2, s2, c3, s3
  real, dimension(MDIM,MDIM) :: R1,R2,R3
  !real, parameter :: pi = 4.0*atan(1.0)

  !
  ! Compute the rotation matrix from x->X
  !
  ! 3-rotation in phi
  c3 = cos(phi(1))
  s3 = sin(phi(1))
  R3 = reshape( (/ c3,s3,0.0,  -s3,c3,0.0,  0.0,0.0,1.0 /), (/3,3/), ORDER=(/2,1/) )

  ! 2-rotation in theta
  c2 = cos(theta(1))
  s2 = sin(theta(1))
  R2 = reshape( (/ c2,0.0,-s2,  0.0,1.0,0.0,   s2,0.0,c2 /), (/3,3/), ORDER=(/2,1/) )

  ! 1-rotation in eta
  c1 = -cos(eta(1)) ! = cos(eta+pi)
  s1 = -sin(eta(1)) ! = sin(eta+pi)
  R1 = reshape( (/ 1.0,0.0,0.0,  0.0,c1,s1,   0.0,-s1,c1 /), (/3,3/), ORDER=(/2,1/) )

  ! compute the rotation matrix
  ! R = (R1*R2*R3)^T
  R = transpose( matmul(matmul(R1,R2),R3) )

  !
  ! Compute \dot{R}
  !
  Rd =  reshape( (/ &
       -(c2*phi(2)*s3) - c3*s2*theta(2), &
       c1*c3*(-phi(2) + eta(2)*s2) + s1*(eta(2)*s3 - phi(2)*s2*s3 + c2*c3*theta(2)), &
       c1*(eta(2) - phi(2)*s2)*s3 + c3*(phi(2)*s1 - eta(2)*s1*s2 + c1*c2*theta(2)), &
       c2*c3*phi(2) - s2*s3*theta(2), &
       c3*s1*(-eta(2) + phi(2)*s2) + s3*(-(c1*phi(2)) + c1*eta(2)*s2 + c2*s1*theta(2)), &
       s1*(phi(2) - eta(2)*s2)*s3 + c1*(-(c3*eta(2)) + c3*phi(2)*s2 + c2*s3*theta(2)), &
       -(c2*theta(2)), &
       c1*c2*eta(2) - s1*s2*theta(2),&
       -(c2*eta(2)*s1) - c1*s2*theta(2) &
       /), (/3,3/), ORDER=(/2,1/) )

  !
  ! Compute \ddot{R}
  !
  Rdd = reshape( (/ &
       -(c2*(phi(3)*s3 + c3*(phi(2)**2 + theta(2)**2))) + s2*(2*phi(2)*s3*theta(2) - c3*theta(3)), &
       c1*((eta(2)**2 + phi(2)**2 - 2*eta(2)*phi(2)*s2)*s3 + c3*(-phi(3) + eta(3)*s2 + 2*c2*eta(2)*theta(2))) &
       - s1*(s3*(-eta(3) + phi(3)*s2 + 2*c2*phi(2)*theta(2)) + c3*(-2*eta(2)*phi(2) + eta(2)**2*s2 + s2*(phi(2)**2 & 
       + theta(2)**2) - c2*theta(3))), &
       s1*(-((eta(2)**2 + phi(2)**2 - 2*eta(2)*phi(2)*s2)*s3) + c3*(phi(3) - eta(3)*s2 - 2*c2*eta(2)*theta(2))) & 
       - c1*(s3*(-eta(3) + phi(3)*s2 + 2*c2*phi(2)*theta(2)) + c3*(-2*eta(2)*phi(2) + eta(2)**2*s2 + s2*(phi(2)**2 & 
       + theta(2)**2) - c2*theta(3))), &
       c3*(c2*phi(3) - 2*phi(2)*s2*theta(2)) - s3*(c2*(phi(2)**2 + theta(2)**2) + s2*theta(3)), &
       -(c1*(c3*(eta(2)**2 + phi(2)**2 - 2*eta(2)*phi(2)*s2) + s3*(phi(3) - eta(3)*s2 - 2*c2*eta(2)*theta(2)))) &
       - s1*(-2*eta(2)*phi(2)*s3 + c3*(eta(3) - phi(3)*s2 - 2*c2*phi(2)*theta(2)) + s2*s3*(eta(2)**2 + phi(2)**2 & 
       + theta(2)**2) - c2*s3*theta(3)), &
       s1*(c3*(eta(2)**2 + phi(2)**2 - 2*eta(2)*phi(2)*s2) + s3*(phi(3) - eta(3)*s2 - 2*c2*eta(2)*theta(2))) &
       - c1*(-2*eta(2)*phi(2)*s3 + c3*(eta(3) - phi(3)*s2 - 2*c2*phi(2)*theta(2)) + s2*s3*(eta(2)**2 + phi(2)**2 &
       + theta(2)**2) - c2*s3*theta(3)), &
       s2*theta(2)**2 - c2*theta(3), &
       c1*c2*eta(3) - 2*c1*eta(2)*s2*theta(2) - s1*(c2*(eta(2)**2 + theta(2)**2) + s2*theta(3)), &
       -(c2*eta(3)*s1) + 2*eta(2)*s1*s2*theta(2) - c1*(c2*(eta(2)**2 + theta(2)**2) + s2*theta(3)) &
       /), (/3,3/), ORDER=(/2,1/) )

end subroutine sm_pk_flapping_BermanWang_RotMat

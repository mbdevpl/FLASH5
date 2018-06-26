
#include "Flash.h"
#include "constants.h"

subroutine sm_pk_Whirl(time,n,Xe,v,vd,vdd,parameters)
  ! This subroutine computes the prescribed coordinates {v, vd, vdd} for a wing.
  ! It assumes that the base wing is more or less defined in the x-y plane with y being forward, and x is the span of the wing.
  ! The inputs:
  !     n: number of points to computer
  !    Xe: (3,n) grid nodes (ref configuration) to move
  !     v: (3,n) global displacement of the nodes under the kinematic transformation (not position)
  !    vd: (3,n) global velocities of the prescribed nodes
  !   vdd: (3,n) global accelerations of the prescribed nodes
  !   
  ! The parameters(:) array is as follows
  ! 1-3: pivot point
  ! 4-6: unit vector of the axis of rotation
  ! 7:  angular velocity
  ! 8:  damped time

  use sm_pk_interface, only: sm_pk_Whirl_functions, sm_pk_Whirl_RotMat
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
  real, dimension(NDIM,NDIM) :: R, Rd, Rdd, eye, Ri
  real, dimension(NDIM)      :: theta, nhat, a, pivot
  real    :: omega, tau
  integer :: i,j

  v   = 0.0
  vd  = 0.0
  vdd = 0.0

#if NDIM == MDIM

  ! Identity Matrix
  eye = reshape( (/ 1.,0.,0.,  0.,1.,0.,  0.,0.,1. /), (/3,3/) , ORDER=(/2,1/) )

  ! Unload needed parameters
  pivot = (/ parameters(1), parameters(2), parameters(3) /)
  nhat  = (/ parameters(4), parameters(5), parameters(6) /)

  ! get current angles:
  call sm_pk_Whirl_functions(time, theta, parameters)

  ! get rotation matrix
  call sm_pk_Whirl_RotMat(theta, nhat, R, Rd, Rdd)

  !
  ! loop over each node of RestSurface and set (v, vd, vdd) -> (qn, qdn, qddn)
  !
  ! if a = X - pivot
  ! then displacement q = R*a - a = (R-I)*a
  !      q   = (R-I)*a
  !      qd  = Rd*a
  !      qdd = Rdd*a

  Ri = R-eye

  do j = 1,n
     
     do i = 1,NDIM
        a(i) = Xe(i,j) - pivot(i)
     end do

     V(1:NDIM,j)   = matmul(Ri , a )
     Vd(1:NDIM,j)  = matmul(Rd , a )
     Vdd(1:NDIM,j) = matmul(Rdd, a )

  end do
#endif

  return

end subroutine sm_pk_Whirl

subroutine sm_pk_Whirl_functions(time, theta, parameters)
  implicit none
  real, intent(in)                     :: time
  real, intent(out), dimension(MDIM)   :: theta
  real, intent(in),  dimension(:)      :: parameters


  ! Internal Variables
  real :: pivot(MDIM), nhat(MDIM), omega, tau
  real :: damp, dampd, dampdd
  
  ! Unload parameters
  pivot = (/ parameters(1), parameters(2), parameters(3) /)
  nhat  = (/ parameters(4), parameters(5), parameters(6) /)
  omega = parameters(7)
  tau   = parameters(8)

  ! Compute the 'damped start' stuff
  damp   = 1. - exp(-time/tau)
  dampd  = (1. - damp)/tau       ! = exp(-time/tau)/tau
  dampdd = -dampd/tau            ! = -exp(-time/tau)/tau**2

  ! compute the values of 
  theta(1) = damp*omega*time
  theta(2) = damp*omega + dampd*omega*time
  theta(3) = 2.*omega*dampd + dampdd*omega*time

  return
end subroutine sm_pk_Whirl_functions

subroutine sm_pk_Whirl_RotMat(theta, nhat, R, Rd, Rdd)
  implicit none
  real, intent(in),  dimension(MDIM)      :: theta  ! angle of rotation (and its derivatives in time)
  real, intent(in),  dimension(MDIM)      :: nhat   ! unit vector
  real, intent(out), dimension(MDIM,MDIM) :: R, Rd, Rdd

  ! Local variables
  real, dimension(MDIM, MDIM) :: eye, A, B, Rtheta, Rtheta2
  real :: s,c

  ! Define some local vars
  s = sin(theta(1))
  c = cos(theta(1))
  eye = reshape( (/ 1.,0.,0.,  0.,1.,0.,  0.,0.,1. /), (/3,3/) , ORDER=(/2,1/) )

  ! spin(nhat)
  A   = reshape( (/       0.,-nhat(3), nhat(2),     &
                     nhat(3),      0.,-nhat(1),     &
                    -nhat(2), nhat(1),       0. /), &
                    (/3,3/) , ORDER=(/2,1/) )

  ! nhat*nhat^T - I
  B   = reshape( (/ nhat(1)**2 - 1., nhat(1)*nhat(2), nhat(1)*nhat(3),     &
                    nhat(1)*nhat(2), nhat(2)**2 - 1., nhat(2)*nhat(3),     &
                    nhat(1)*nhat(3), nhat(2)*nhat(3), nhat(3)**2 - 1. /) , &
                    (/3,3/), ORDER=(/2,1/) )
  
  ! Build Rotation matrix R (using Rodrigues formula)
  R = eye + s*A + (1.-c)*B

  ! dR/dtheta
  Rtheta = c*A + s*B

  ! dR/dt = dR/dtheta * dtheta/dt
  Rd = Rtheta*theta(2)

  ! d^2R/dtheta^2
  Rtheta2 = -s*A + c*B

  ! d^2R/dt^2 = d^2theta/dt^2 * dR/dtheta + (dtheta/dt)^2 * d^2R/dtheta^2 
  Rdd = theta(3)*Rtheta + theta(2)**2 * Rtheta2

  return
end subroutine sm_pk_Whirl_RotMat

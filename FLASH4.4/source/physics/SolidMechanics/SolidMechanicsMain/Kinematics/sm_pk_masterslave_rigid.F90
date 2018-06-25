!
!
!



subroutine sm_pk_masterslave_rigid(maxdofs,nfix,xm,qmstr,qdmstr,qddmstr, &
                                   TNB,NwB_N,NaB_N,Xi,XiB,v,vd,vdd)

  use sm_Misc_interface, only : sm_crossProd

  implicit none
#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

  integer, intent(in) :: nfix,maxdofs
  real, intent(in)  :: xm(MDIM,1)
  real, intent(in)  :: qmstr(maxdofs),qdmstr(maxdofs),qddmstr(maxdofs)
  real, intent(in)  :: TNB(MDIM,MDIM),NwB_N(MDIM,1),NaB_N(MDIM,1)
  real, intent(in)  :: Xi(MDIM,nfix),XiB(MDIM,nfix)
  real, intent(out) :: v(maxdofs,nfix),vd(maxdofs,nfix),vdd(maxdofs,nfix)

  ! Local Vars
  real :: u(MDIM,1),ud(MDIM,1),udd(MDIM,1)
  real :: WP(MDIM,1),WPP(MDIM,1),AP(MDIM,1)

  real ::  xPB_B(MDIM,1),xPB_N(MDIM,1)

  integer :: i

  !!! Initialize vars:
  v=0.;   vd=0.;    vdd=0.;
  u=0.;   ud=0.;    udd=0.;

  ! Linear Displacements, velocities and accelerations of master point
  u(1:NDIM,1)  =  qmstr(1:NDIM)
  ud(1:NDIM,1) = qdmstr(1:NDIM)
  udd(1:NDIM,1)=qddmstr(1:NDIM) 

  ! Obtain Linear displacements, velocities and accelerations of slave points:
  xPB_B = 0.
  do i=1,nfix
     xPB_B(1:NDIM,1) = XiB(1:NDIM,i)
     xPB_N(1:NDIM,1) = matmul(TNB(1:NDIM,1:NDIM),xPB_B(1:NDIM,1)) ! Vector from B tp Point P in the N frame.

     ! Displacement of point P uP = (xB0 - xP0) + uB + xpb_N  
     v(1:NDIM,i)  = (xm(1:NDIM,1) - Xi(1:NDIM,i)) + u(1:NDIM,1) + xPB_N(1:NDIM,1)

#if NDIM == MDIM
     ! Velocity of Point P uPd = uBd + NwB_N x xpb_N
     WP(1:NDIM,1) = sm_crossProd(NwB_N(1:NDIM,1),xPB_N(1:NDIM,1))
     vd(1:NDIM,i) = ud(1:NDIM,1) + WP(1:NDIM,1)

     ! Acceleration of Point P uPdd = uBdd + NwB_N x NwB_N x xpb_N + NaB_N x xpb_N
     WPP(1:NDIM,1) = sm_crossProd(NwB_N(1:NDIM,1),WP(1:NDIM,1))
     AP(1:NDIM,1)  = sm_crossProd(NaB_N(1:NDIM,1),xPB_N(1:NDIM,1))
     vdd(1:NDIM,i) = udd(1:NDIM,1) + WPP(1:NDIM,1) + AP(1:NDIM,1)

#else /* 2d */

     ! Velocity of Point P uPd = uBd + NwB_N x xpb_N
     ! ux = ubx - w*ypb
     vd(IAXIS,i) = ud(IAXIS,1) - NwB_N(1,1)*xPB_N(JAXIS,1)
     ! uy = uby + w*xpb
     vd(JAXIS,i) = ud(JAXIS,1) + NwB_N(1,1)*xPB_N(IAXIS,1)

     ! Acceleration of Point P uPdd = uBdd + NwB_N x NwB_N x xpb_N + NaB_N x xpb_N
     ! udx = udbx - w^2*xpb - a*ypb
     vdd(IAXIS,i) = udd(IAXIS,1) - NwB_N(1,1)**2. * xPB_N(IAXIS,1) - NaB_N(1,1)*xPB_N(JAXIS,1)
     ! udy = udby - w^2*ypb + a*xpb 
     vdd(JAXIS,i) = udd(JAXIS,1) - NwB_N(1,1)**2. * xPB_N(JAXIS,1) + NaB_N(1,1)*xPB_N(IAXIS,1)

#endif

  enddo

  return

end subroutine sm_pk_masterslave_rigid

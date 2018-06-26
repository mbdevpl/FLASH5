!!****if* source/physics/IncompNS/IncompNSMain/vardens/IncompNS_data
!!
!! NAME
!!
!!  IncompNS_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE IncompNS_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the Incmp_Navier_Stokes module.
!!
!!***
 
 
module IncompNS_data

#include "Flash.h"
#include "constants.h"

  integer, save :: ins_isgs
  integer, save :: ins_cflflg
  integer, save :: ins_intschm
  integer, save :: ins_nstep
  integer, save :: ins_meshMe
  integer, save :: ins_meshNumProcs
  integer, save :: ins_meshComm
  integer, save :: ins_iConvU

  real, save :: ins_cfl
  real, save :: ins_sigma
  real, save :: ins_invRe
  real, save :: ins_dtspec

  real, save :: ins_gravX, ins_gravY, ins_gravZ, ins_dampC
  real, save :: ins_xDampL, ins_xDampR, ins_yDampL, ins_zDampL

  integer, parameter :: AB2_SCHM = 2
  integer, parameter :: RK3_SCHM = 3

  integer, parameter :: rkstep = 3

  real, save :: ins_alf(rkstep)
  real, save :: ins_gam(rkstep)
  real, save :: ins_rho(rkstep)
  real, save :: ins_alfa, ins_gama, ins_rhoa

  logical, save :: ins_prescorr
  real, save :: ins_prescoeff

  logical, save :: ins_restart

  real, save :: ins_tlevel

  real, save :: ins_Qin
  real, save :: ins_Qout

  logical, save :: ins_predcorrflg

  logical, save :: ins_outflowgridChanged

  real, save :: ins_convvel(LOW:HIGH,MDIM)

  integer, save :: ins_prol_method

  ! Arrays for convective outflow BC
!  real, save :: uvel_x(NGUARD+1,NYB+2*NGUARD*K2D,NZB+2*NGUARD*K3D,LOW:HIGH,MAXBLOCKS)
!  real, save :: vvel_x(NGUARD+1,NYB+2*NGUARD*K2D+K2D,NZB+2*NGUARD*K3D,LOW:HIGH,MAXBLOCKS)
!  real, save :: wvel_x(NGUARD+1,NYB+2*NGUARD*K2D,NZB+2*NGUARD*K3D+K3D,LOW:HIGH,MAXBLOCKS)
  real, save :: uvel_x(NGUARD+1,NYB+2*NGUARD*K2D,NZB+2*NGUARD*K3D,LOW:HIGH,MAXBLOCKS*10)
  real, save :: vvel_x(NGUARD+1,NYB+2*NGUARD*K2D+K2D,NZB+2*NGUARD*K3D,LOW:HIGH,MAXBLOCKS*10)
  real, save :: wvel_x(NGUARD+1,NYB+2*NGUARD*K2D,NZB+2*NGUARD*K3D+K3D,LOW:HIGH,MAXBLOCKS*10)

!  real, save :: uvel_y(NGUARD+1,NXB+2*NGUARD+1,NZB+2*NGUARD*K3D,LOW:HIGH,MAXBLOCKS)
!  real, save :: vvel_y(NGUARD+1,NXB+2*NGUARD,NZB+2*NGUARD*K3D,LOW:HIGH,MAXBLOCKS)
!  real, save :: wvel_y(NGUARD+1,NXB+2*NGUARD,NZB+2*NGUARD*K3D+K3D,LOW:HIGH,MAXBLOCKS)
  real, save :: uvel_y(NGUARD+1,NXB+2*NGUARD+1,NZB+2*NGUARD*K3D,LOW:HIGH,MAXBLOCKS*10)
  real, save :: vvel_y(NGUARD+1,NXB+2*NGUARD,NZB+2*NGUARD*K3D,LOW:HIGH,MAXBLOCKS*10)
  real, save :: wvel_y(NGUARD+1,NXB+2*NGUARD,NZB+2*NGUARD*K3D+K3D,LOW:HIGH,MAXBLOCKS*10)

#if NDIM == 3
!  real, save :: uvel_z(NGUARD+1,NXB+2*NGUARD+1,NYB+2*NGUARD*K2D,LOW:HIGH,MAXBLOCKS)
!  real, save :: vvel_z(NGUARD+1,NXB+2*NGUARD,NYB+2*NGUARD*K2D+1,LOW:HIGH,MAXBLOCKS)
!  real, save :: wvel_z(NGUARD+1,NXB+2*NGUARD,NYB+2*NGUARD*K2D,LOW:HIGH,MAXBLOCKS)
  real, save :: uvel_z(NGUARD+1,NXB+2*NGUARD+1,NYB+2*NGUARD*K2D,LOW:HIGH,MAXBLOCKS*10)
  real, save :: vvel_z(NGUARD+1,NXB+2*NGUARD,NYB+2*NGUARD*K2D+1,LOW:HIGH,MAXBLOCKS*10)
  real, save :: wvel_z(NGUARD+1,NXB+2*NGUARD,NYB+2*NGUARD*K2D,LOW:HIGH,MAXBLOCKS*10)
#endif

end module IncompNS_data

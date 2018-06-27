!     
! File:   sm_pk_fixed.F90
! Author: tim
!
! routine is based on the kine
#include "Flash.h"
#include "constants.h"

subroutine sm_pk_fixed(time,n,Xe,v,vd,vdd,parameters)
  ! This subroutine computes the prescribed coordinates {v, vd, vdd} for a surface constrainted to be zero
  ! The inputs:
  !     n: number of points to computer
  !    Xe: (3,n) grid nodes (ref configuration) to move
  !     v: (3,n) global displacement of the nodes under the kinematic transformation (not position)
  !    vd: (3,n) global velocities of the prescribed nodes
  !   vdd: (3,n) global accelerations of the prescribed nodes
  !   

  implicit none

  !
  ! IO Variables
  !
  real, intent(in)    :: time
  integer, intent(in) :: n
  real, intent(in),  dimension(NDIM,n) :: Xe
  real, intent(out), dimension(NDIM,n) :: v, vd, vdd
  real, intent(in),  dimension(:)      :: parameters

  V   = parameters(1)
  Vd  = 0.
  Vdd = 0.

  return

end subroutine sm_pk_fixed

!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Misc/sm_crossProd_mat
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!
!!***

#include "Flash.h"
#include "constants.h"

function sm_crossProd_mat(a)

#if NDIM != MDIM
  use Driver_interface,            only: Driver_abortFlash
#endif

  implicit none

  real, intent(in) :: a(NDIM)
  real :: sm_crossProd_mat(NDIM,NDIM)

  ! Compute the cross product 3x3 matrix
  sm_crossProd_mat = 0.

#if NDIM == MDIM
  sm_crossProd_mat(1,2) = -a(3)
  sm_crossProd_mat(1,3) =  a(2)
  
  sm_crossProd_mat(2,1) =  a(3)
  sm_crossProd_mat(2,3) = -a(1)

  sm_crossProd_mat(3,1) = -a(2)
  sm_crossProd_mat(3,2) =  a(1)

#else
  call Driver_abortFlash('Function is only defined in 3D <sm_crossProd_mat.F90>')
#endif

  return

end function sm_crossProd_mat

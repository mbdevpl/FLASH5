!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Misc/sm_crossProd
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

function sm_crossProd(a,b)

  implicit none
  real, intent(in) :: a(NDIM),b(NDIM)

#if NDIM == 2
  real :: sm_crossProd
#elif NDIM == 3
  real :: sm_crossProd(NDIM)
#endif

  ! Compute the cross product

#if NDIM == 2

  sm_crossProd = a(1)*b(2) - a(2)*b(1)

#elif NDIM == 3

  sm_crossProd(1) = a(2)*b(3) - a(3)*b(2)

  sm_crossProd(2) = a(3)*b(1) - a(1)*b(3)

  sm_crossProd(3) = a(1)*b(2) - a(2)*b(1)

#endif

  return

end function sm_crossProd

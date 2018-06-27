!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Testing/DryRun/sm_surf_assembleFluidForce_toPoint
!!
!! NAME
!! 
!!
!! SYNOPSIS
!!    Stub function
!!  
!! DESCRIPTION 
!! 
!!
!! ARGUMENTS 
!!
!!***

#include "Flash.h"

subroutine sm_surf_assembleFluidForce_toPoint(ibd, point, force_p, force_v, moment_p, moment_v)
  implicit none

  ! IO variables 
  integer, intent(in) :: ibd
  real, dimension(NDIM), intent(in) :: point  ! location to compute moments about
  real, dimension(NDIM), intent(out):: force_p, force_v, moment_p, moment_v

  force_v  = 0.
  force_p  = 0.
  moment_v = 0.
  moment_p = 0.

  return

end subroutine sm_surf_assembleFluidForce_toPoint

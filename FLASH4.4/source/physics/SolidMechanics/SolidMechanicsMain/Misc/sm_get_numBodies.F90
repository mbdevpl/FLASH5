!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Misc/sm_get_numBodies
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

subroutine sm_get_numBodies(numBodies)
  use SolidMechanics_data, only: sm_NumBodies
  implicit none
  integer, intent(out) :: numBodies
  
  numBodies = sm_NumBodies

  return

end subroutine sm_get_numBodies

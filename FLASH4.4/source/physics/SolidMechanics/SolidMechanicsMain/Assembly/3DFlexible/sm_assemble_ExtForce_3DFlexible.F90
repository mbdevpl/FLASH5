!     
! File:  
! Author: tim

#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

subroutine sm_assemble_ExtForce_3DFlexible(ibd, time)

  use SolidMechanics_data, only: sm_BodyInfo, sm_structure
  use Driver_interface, only : Driver_abortFlash

  implicit none
    
  !IO Variables
  integer, intent(in) :: ibd ! body number
  real, intent(in)    :: time
    
  ! Define internal variables
  type(sm_structure), pointer :: body
    
  body => sm_BodyInfo(ibd)
    
  ! If there are body forces defined on the reference configuration, put them here

  ! Gravity Forces here
  if( body%gravity_flag == SM_TRUE ) then
     call Driver_abortFlash("ToDo: implement gravity on the 3DFlexible body.")
  end if

  
  return
    
end subroutine sm_assemble_ExtForce_3DFlexible

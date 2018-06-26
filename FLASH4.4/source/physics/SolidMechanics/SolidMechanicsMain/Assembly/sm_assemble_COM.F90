!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Assembly/sm_assemble_COM
!!
!! NAME
!!  
!!
!! SYNOPSIS
!!
!!  
!! DESCRIPTION 
!!   subroutine to compute the center of mass of a body
!!  
!!
!! ARGUMENTS 
!!
!!***

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_assemble_COM(ibd)
  
  use SolidMechanics_data,   only : sm_BodyInfo, sm_structure
  use Driver_interface,      only : Driver_abortFlash
  use sm_assemble_interface, only : sm_assemble_COM_3DFlexible, &
                                    sm_assemble_COM_rigid, sm_assemble_COM_rbc
  implicit none
    
  !IO Variables
  integer, intent(in) :: ibd ! body number
    
  ! Define internal variables
  type(sm_structure), pointer :: body

  ! Point to body
  body => sm_BodyInfo(ibd)


  ! Initialize the Center of Mass
  body%COM = 0.
    
  ! Load and allocate based on BodyType
   select case( body%BodyType )

   case( BODYTYPE_RIGID ) 
      call sm_assemble_COM_rigid(ibd)

   case( BODYTYPE_2DFLEXIBLE )
      call Driver_abortFlash('BodyType not yet implemented')

   case( BODYTYPE_3DFLEXIBLE )
      call sm_assemble_COM_3DFlexible(ibd, SM_IOPT_QN )

   case( BODYTYPE_RBC )
      call sm_assemble_COM_rbc(ibd)
     
   case default
      call Driver_abortFlash("BodyType not yet implemented <sm_assemble_COM.F90>")

   end select

  return
    
end subroutine sm_assemble_COM

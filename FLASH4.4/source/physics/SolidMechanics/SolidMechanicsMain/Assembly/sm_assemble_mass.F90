!     
! File:   physics/SolidMechanics/SolidMechanicsMain/Assembly/sm_assemble_mass.F90
! Author: tim
! Description: assembles mass matrix
!     
!

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_assemble_mass(ibd, iopt)
  
  use SolidMechanics_data, only: sm_nen,  sm_BodyInfo, sm_structure
  use Driver_interface, only : Driver_abortFlash
  use sm_assemble_interface, only: sm_assemble_mass_3DFlexible, &
                                   sm_assemble_mass_rigid, sm_assemble_mass_rbc

  implicit none
    
  !IO Variables
  integer, intent(in)    :: ibd ! body number
  integer, optional, intent(in) :: iopt
    
  ! Define internal variables
  type(sm_structure), pointer :: body
  integer :: flag = SM_FALSE

  if( PRESENT( iopt ) ) flag = iopt
    
  body => sm_BodyInfo(ibd)
    
  ! Load and allocate based on BodyType
   select case( body%BodyType )

   case( BODYTYPE_RIGID ) 
      call sm_assemble_mass_rigid(ibd,flag)

   case( BODYTYPE_2DFLEXIBLE )
      call Driver_abortFlash('BodyType not yet implemented')

   case( BODYTYPE_3DFLEXIBLE )
      call sm_assemble_mass_3DFlexible(ibd,flag)

   case( BODYTYPE_RBC )
      call sm_assemble_mass_RBC(ibd,flag)

   case default
      call Driver_abortFlash("BodyType not yet implemented <sm_assemble_mass.F90>")

   end select

  return
    
end subroutine sm_assemble_mass

!     
! File:   physics/SolidMechanics/SolidMechanicsMain/Assembly/sm_assemble_ExtForce.F90
! Author: tim
! Description: assembles External force vector
!     
!

#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

subroutine sm_assemble_ExtForce(ibd, time)
  
  use SolidMechanics_data, only:  sm_BodyInfo, sm_structure
  use Driver_interface, only : Driver_abortFlash
  use sm_assemble_interface, only: sm_assemble_ExtForce_3DFlexible, &
                                   sm_assemble_ExtForce_rigid,      &
                                   sm_assemble_ExtForce_rbc

  use sm_surf_interface, only: sm_surf_assembleFluidForce, &
                               sm_surf_assembleFluidForce_rigid 

  implicit none
    
  !IO Variables
  integer, intent(in) :: ibd ! body number
  real,    intent(in) :: time
    
  ! Define internal variables
  type(sm_structure), pointer :: body
  
  body => sm_BodyInfo(ibd)

  ! Initalize Hs, always add to Hs inside subroutines below
  body%Hs(:) = 0.
  
  ! Load and allocate based on BodyType
  select case( body%BodyType )

  case( BODYTYPE_RIGID ) 
     call sm_assemble_ExtForce_rigid(ibd,time)

     ! Assemble the Fluid Forces for rigid
     call sm_surf_assembleFluidForce_rigid(ibd)

#if NDIM == MDIM

  case( BODYTYPE_3DFLEXIBLE )
     call sm_assemble_ExtForce_3DFlexible(ibd,time)

     ! Assemble the Fluid Forces
     ! if you don't need to do this, then stub the function (see SolidMechanicsMain/Testing/DryRun)
     call sm_surf_assembleFluidForce(ibd)

  case( BODYTYPE_RBC )
     call sm_assemble_ExtForce_RBC(ibd,time)

     ! Assemble the Fluid Forces
     ! if you don't need to do this, then stub the function (see SolidMechanicsMain/Testing/DryRun)
     call sm_surf_assembleFluidForce(ibd)

#else /* 2D */

  case( BODYTYPE_2DFLEXIBLE )
     call Driver_abortFlash('BodyType 2DFLEXIBLE not yet implemented <sm_assemble_ExtForce.F90>')


#endif

  case default
     call Driver_abortFlash("BodyType not yet implemented  <sm_assemble_ExtForce.F90>")

  end select



  return
    
end subroutine sm_assemble_ExtForce

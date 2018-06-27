!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/sm_surf_currentFluidForce
!!
!! NAME
!!  
!!
!! SYNOPSIS
!!
!!  
!! DESCRIPTION 
!!   subroutine to compute the current fluid forces on body ibd
!!  
!!
!! ARGUMENTS 
!!
!!***

#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

subroutine sm_surf_currentFluidForce(ibd, ndofs, Hs_pres, Hs_visc )
  
  use SolidMechanics_data, only : sm_BodyInfo, sm_structure
  use Driver_interface,    only : Driver_abortFlash
  use sm_surf_interface,   only : sm_surf_currentFluidForce_3DFlex, &
                                  sm_surf_currentFluidForce_rigid
  implicit none
    
  !IO Variables
  integer, intent(in)  :: ibd ! body number
  integer, intent(in)  :: ndofs
  real,    intent(out) :: Hs_pres(ndofs), Hs_visc(ndofs)
  
    
  ! Define internal variables
  type(sm_structure), pointer :: body

  ! Point to body
  body => sm_BodyInfo(ibd)


  ! Initialize the force containers
  Hs_pres = 0.
  Hs_visc = 0.
    
  ! Load and allocate based on BodyType
   select case( body%BodyType )

   case( BODYTYPE_RIGID ) 
      call sm_surf_currentFluidForce_rigid(ibd, ndofs, Hs_pres, Hs_visc )

   case( BODYTYPE_2DFLEXIBLE )
      call Driver_abortFlash('BodyType 2DFLEXIBLE not yet implemented <sm_surf_currentFluidForce.F90>')

   case( BODYTYPE_3DFLEXIBLE )
      call sm_surf_currentFluidForce_3DFlex(ibd, ndofs, Hs_pres, Hs_visc )

   case( BODYTYPE_RBC )
      call sm_surf_currentFluidForce_3DFlex(ibd, ndofs, Hs_pres, Hs_visc )
     
   case default
      call Driver_abortFlash("BodyType not yet implemented <sm_surf_currentFluidForce.F90>")

   end select

  return
  
end subroutine sm_surf_currentFluidForce

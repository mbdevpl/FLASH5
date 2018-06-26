
module sm_surf_interface
#include "Flash.h"
#include "SolidMechanics.h"
#include "constants.h"

  implicit none
  interface 
     subroutine sm_surf_MapParticles(ibd)
       implicit none
       integer, intent(in) :: ibd
     end subroutine sm_surf_MapParticles
  end interface

  interface
     subroutine sm_surf_CountParticles_wsIEN(ibd, numPart)
       implicit none
       integer, intent(in) :: ibd
       integer, intent(out) :: numPart
     end subroutine sm_surf_CountParticles_wsIEN
  end interface

  interface
     subroutine sm_surf_assembleFluidForce(ibd)
       implicit none
       integer, intent(in) :: ibd
     end subroutine sm_surf_assembleFluidForce
  end interface
         
  interface
     subroutine sm_surf_assembleFluidForce_toPoint(ibd, point, &
          force_pres, force_visc, moment_pres, moment_visc)
       implicit none
       integer, intent(in) :: ibd
       real, dimension(NDIM), intent(in) :: point  ! location to compute moments about
       real, dimension(NDIM), intent(out):: force_pres, force_visc, moment_pres, moment_visc
     end subroutine sm_surf_assembleFluidForce_toPoint
  end interface

  interface
     subroutine sm_surf_assembleFluidForce_rigid(ibd)
       implicit none
       integer, intent(in) :: ibd
     end subroutine sm_surf_assembleFluidForce_rigid
  end interface

  interface
     subroutine sm_surf_repoParticlesBC(ibd)
       implicit none
       integer, intent(in) :: ibd
     end subroutine sm_surf_repoParticlesBC
  end interface

  interface
     subroutine sm_surf_currentFluidForce(ibd, ndofs, Hs_pres, Hs_visc)
       implicit none
       integer, intent(in)  :: ibd, ndofs
       real,    intent(out) :: Hs_pres(ndofs), Hs_visc(ndofs)
     end subroutine sm_surf_currentFluidForce
  end interface

  interface
     subroutine sm_surf_currentFluidForce_3DFlex(ibd, ndofs, Hs_pres, Hs_visc)
       implicit none
       integer, intent(in)  :: ibd, ndofs
       real,    intent(out) :: Hs_pres(ndofs), Hs_visc(ndofs)
     end subroutine sm_surf_currentFluidForce_3DFlex
  end interface

  interface
     subroutine sm_surf_currentFluidForce_rigid(ibd, ndofs, Hs_pres, Hs_visc)
       implicit none
       integer, intent(in)  :: ibd, ndofs
       real,    intent(out) :: Hs_pres(ndofs), Hs_visc(ndofs)
     end subroutine sm_surf_currentFluidForce_rigid
  end interface


end module sm_surf_interface

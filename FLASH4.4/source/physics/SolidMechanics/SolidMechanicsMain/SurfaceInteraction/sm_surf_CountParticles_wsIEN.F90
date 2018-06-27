!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SurfaceInteraction/sm_surf_CountParticles_wsIEN
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!! 
!! Counts the number of marker particles for one solid body.
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_surf_CountParticles_wsIEN(ibd, numPart)

  use SolidMechanics_data, only : sm_BodyInfo
  use Grid_data, only : gr_minCellSizes
  use Grid_interface, ONLY : Grid_getMinCellSize
  use Driver_interface, only : Driver_abortFlash
  use sm_element_interface, only: sm_el01_countParticles, &
          sm_el02_countParticles, sm_el10_countParticles
  implicit none

  ! IO Variables
  integer, intent(in)  :: ibd
  integer, intent(out) :: numPart

  ! Internal Variables
  integer :: e, nXi, nEta, ptelem
  real :: Dmin

#ifdef FLASH_GRID_PARAMESH 
  ! Get the Eulerian cell size
  call Grid_getMinCellSize(Dmin)
#else /* Uniform Grid */
  ! map_use_dmin is true if the max surface element
  ! side length is divided the Eulerian grid Dmin.
  ! In case surface elements are are skewed in the
  ! same direction as Eulerian grid Dmax can be used.
  if (sm_bodyinfo(ibd)%map_use_dmin .eq. SM_TRUE) then
    call Grid_getMinCellSize(Dmin)
  else
    Dmin = maxval(gr_minCellSizes(IAXIS:NDIM))
  endif
#endif

  ! Loop over wet-surface elements
  numPart = 0

  do e = 1,sm_bodyinfo(ibd)%ws_nel

     select case( sm_bodyInfo(ibd)%ws_eltype(e) )
        
#if NDIM == MDIM

        case( THREE_NODE_TRIANGLE )
           call sm_el02_countParticles(sm_bodyInfo(ibd), e, Dmin, nXi, nEta, ptelem, SM_IOPT_QN)
           
        case( NINE_NODE_QUADRILATERAL  )
           call sm_el10_countParticles(sm_bodyInfo(ibd), e, Dmin, nXi, nEta, ptelem, SM_IOPT_QN)

#else /* 2D */

        case( TWO_NODE_LINE )
           call sm_el01_countParticles(sm_bodyInfo(ibd), e, Dmin, nXi, nEta, ptelem, SM_IOPT_QN)

#endif

        case default
           call Driver_abortFlash("Element type not recognized.")

     end select

     sm_bodyInfo(ibd)%ws_nXi(e)   = nXi
     sm_bodyInfo(ibd)%ws_nEta(e)  = nEta
     sm_bodyInfo(ibd)%ws_ptelem(e)= ptelem

     numPart = numPart + ptelem

  end do
  
  return

end subroutine sm_surf_CountParticles_wsIEN

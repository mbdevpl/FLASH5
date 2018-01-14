!!****if* source/Grid/GridMain/AMR/Amrex/Grid_updateRefinement
!!
!! NAME
!!
!!  Grid_updateRefinement
!!
!! SYNOPSIS
!!  
!!  call Grid_updateRefinement(integer(IN)            :: nstep,
!!                             real(IN)               :: time,
!!                             logical(OUT), optional :: gridChanged)
!!
!! DESCRIPTION
!!
!!  If the indicated step qualifies as a refinement step, then this routine
!!    (1) fills all guardcells at all levels with EoS run on the guardcells and
!!    (2) triggers AMReX to execute grid refinement.
!!
!!  Note that all EoS runs are done in the mode specified by the eosMode
!!  runtime parameter.
!!
!!  Note also that a step qualifies as a refinement step if the given step
!!  number is a multiple of the nrefs runtime parameter.
!!
!!  AMReX has FLASH identify blocks requiring refinement via the 
!!  gr_markRefineDerefineCallback routine.
!!
!!  After the refinement, AMReX fills cell-centered data in newly-created child
!!  blocks via conservative linear interpolation from the coarse parents.
!!  This step also includes filling guardcells.  Finally, the
!!  EoS is called on the block interiors to make them thermodynamically
!!  consistent. 
!!
!!  Presently, this routine does not alter any particle information.
!!
!! ARGUMENTS
!!
!!  nstep - current step number
!!  time  - current evolution time
!!  gridChanged - returns TRUE if grid may actually have changed.
!!
!! SEE ALSO
!!  Grid_fillGuardCells
!!  gr_markRefineDerefineCallback
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine Grid_updateRefinement(nstep, time, gridChanged)
  use amrex_amrcore_module, ONLY : amrex_regrid

  use Grid_interface,       ONLY : Grid_fillGuardCells
  use Grid_data,            ONLY : gr_nrefs, &
!                                   gr_maxRefine, &
                                   gr_refine_var, &
                                   gr_numRefineVars, &
                                   gr_eosMode
  use Timers_interface,     ONLY : Timers_start, Timers_stop
 
  implicit none

  integer, intent(in)            :: nstep
  real,    intent(in)            :: time
  logical, intent(out), OPTIONAL :: gridChanged

  integer, parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  
  logical, save :: gcMaskArgsLogged = .FALSE.
  
  logical :: gcMask(maskSize)
  integer :: iref
  integer :: i

  ! We only consider refinements every nrefs timesteps.
  if (mod(nstep, gr_nrefs) == 0) then
#ifdef DEBUG_GRID
     write(*,'(A,I4,A,E9.3)') "[Grid_updateRefinement] AMReX Regridding @ step=", & 
                              nstep, " / time = ", time
#endif
 
     ! DEV: TODO Allow for updates to gr_maxRefine here?
     ! Does this work for initialization and all callbacks?
!     if(gr_lrefineMaxRedDoByTime) then
!        call gr_markDerefineByTime()
!     end if
!     
!     if(gr_lrefineMaxByTime) then
!        call gr_setMaxRefineByTime()
!     end if
!
!    if(gr_lrefineMaxRedDoByLogR) then
!       call gr_unmarkRefineByLogRadius(gr_lrefineCenterI, &
!                                       gr_lrefineCenterJ, &
!                                       gr_lrefineCenterK)
!    end if

     gcMask = .FALSE.
     do i = 1, gr_numRefineVars
        iref = gr_refine_var(i)
        gcMask(iref) = (iref > 0)
     end do
     gcMask(NUNK_VARS+1:min(maskSize,NUNK_VARS+NDIM*NFACE_VARS)) = .TRUE.

     call Timers_start("Grid_updateRefinement")
     call Grid_fillGuardCells(CENTER, ALLDIR, &
                              eosMode=gr_eosMode, doEos=.TRUE., &
                              maskSize=maskSize, mask=gcMask, &
                              makeMaskConsistent=.TRUE., &
                              doLogMask=(.NOT. gcMaskArgsLogged))
     call amrex_regrid(0, time)
     call Timers_stop("Grid_updateRefinement")

     ! Only log on the first call
     gcMaskArgsLogged = .TRUE.

    ! DEV: TODO What happens with particles here?

     if (present(gridChanged)) then
        ! DEV: FIXME: Shouldn't this actually check if AMReX
        ! decided to regrid based on the contents of the physical quantities?
        gridChanged = .TRUE.
     end if
  else
     if (present(gridChanged)) then
        gridChanged = .FALSE.
     end if
  end if

end subroutine Grid_updateRefinement


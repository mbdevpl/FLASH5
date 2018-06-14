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
!!    (1) restricts data from leaf blocks down to all ancestors,
!!    (2) fills all guardcells at all levels,
!!    (3) runs EoS on all interiors and all guardcells, and
!!    (4) triggers AMReX to execute grid refinement.
!!
!!  It is assumed that the data in all leaf block interiors is correct and that
!!  EoS has been run for these.  No assumptions are made about the quality of
!!  guardcell data nor ancestor blocks.
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
  use amrex_amrcore_module,      ONLY : amrex_regrid, &
                                        amrex_get_finest_level, &
                                        amrex_geom, &
                                        amrex_ref_ratio
  use amrex_fillpatch_module,    ONLY : amrex_fillpatch
  use amrex_interpolater_module, ONLY : amrex_interp_cell_cons
  
  use Grid_interface,            ONLY : Grid_fillGuardCells, &
                                        Grid_getBlkPtr, Grid_releaseBlkPtr, &
                                        Grid_getLeafIterator, &
                                        Grid_releaseLeafIterator
  use Grid_data,                 ONLY : gr_nrefs, &
!                                        gr_maxRefine, &
                                        gr_refine_var, &
                                        gr_numRefineVars, &
                                        gr_eosMode, &
                                        gr_amrexDidRefinement, &
                                        lo_bc_amrex, hi_bc_amrex
  use gr_interface,              ONLY : gr_getBlkIterator, &
                                        gr_releaseBlkIterator
  use gr_amrexInterface,         ONLY : gr_primitiveToConserve, &
                                        gr_conserveToPrimitiveLevel, &
                                        gr_fillPhysicalBC, &
                                        gr_averageDownLevels
  use gr_physicalMultifabs,      ONLY : unk
  use leaf_iterator,             ONLY : leaf_iterator_t
  use block_metadata,            ONLY : block_metadata_t
  use Eos_interface,             ONLY : Eos_wrapped
  use Timers_interface,          ONLY : Timers_start, Timers_stop
 
  implicit none

  integer, intent(in)            :: nstep
  real,    intent(in)            :: time
  logical, intent(out), OPTIONAL :: gridChanged

  integer, parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS

  logical, save :: gcMaskArgsLogged = .FALSE.

  integer :: finest_level
!  logical :: gcMask(maskSize)
!  integer :: iref
  integer :: lev
  integer :: i

  type(leaf_iterator_t)          :: itor
  type(block_metadata_t)         :: blockDesc
  real,                  pointer :: solnData(:, :, :, :) => null()

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

!     gcMask = .FALSE.
!     do i = 1, gr_numRefineVars
!        iref = gr_refine_var(i)
!        DEV: Does this mean that iref could be below the array's lbound?
!        gcMask(iref) = (iref > 0)
!     end do
!     gcMask(NUNK_VARS+1:min(maskSize,NUNK_VARS+NDIM*NFACE_VARS)) = .TRUE.

     call Timers_start("Grid_updateRefinement")
     
     !!!!! POPULATE ALL BLOCKS AT ALL LEVELS WITH CONSERVED-FORM DATA
     ! We are only concerned with data on interior at this point
     call Grid_getLeafIterator(itor, tiling=.FALSE.)
     do while (itor%is_valid())
       call itor%blkMetaData(blockDesc)
       call gr_primitiveToConserve(blockDesc)

       call itor%next()
     end do
     call Grid_releaseLeafIterator(itor)

     ! Restrict data from leaves to coarser blocks
     call gr_averageDownLevels

     !!!!! POPULATE GUARDCELLS IN ALL BLOCKS
     ! DEV: TODO Confirm with AMReX team if non-parent ancestor blocks can
     ! influence refinement decisions.  If not, we just need to GC fill
     ! on all levels with leaf blocks.
     lev = 0
     call amrex_fillpatch(unk(lev), 1.0d0, unk(lev), &
                                    0.0d0, unk(lev), &
                                    amrex_geom(lev), gr_fillPhysicalBC, &
                                    0.0d0, UNK_VARS_BEGIN, &
                                    UNK_VARS_BEGIN, NUNK_VARS)

     finest_level = amrex_get_finest_level()
     do lev=1, finest_level
        call amrex_fillpatch(unk(lev), 1.0d0, unk(lev-1), &
                                       0.0d0, unk(lev-1), &
                                       amrex_geom(lev-1), gr_fillPhysicalBC, &
                                       1.0e0, unk(lev  ), &
                                       0.0d0, unk(lev  ), &
                                       amrex_geom(lev  ), gr_fillPhysicalBC, &
                                       0.0d0, UNK_VARS_BEGIN, &
                                       UNK_VARS_BEGIN, NUNK_VARS, &
                                       amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                       lo_bc_amrex, hi_bc_amrex) 
     end do

     ! Let AMReX regrid with callbacks
     gr_amrexDidRefinement = .FALSE.
     call amrex_regrid(0, time)

     ! Revert variables to primitive form if necessary
     do lev = 0, finest_level
       call gr_conserveToPrimitiveLevel(lev, .TRUE.)
     end do

     if (gr_amrexDidRefinement) then
       ! We don't know which leaf blocks were created, so run EoS 
       ! on all blocks
       ! DEV: TODO: make gr_amrexDidRefinement into an array so that
       ! we only run EoS on those levels that were really changed.
       call Grid_getLeafIterator(itor)
       do while (itor%is_valid())
          call itor%blkMetaData(blockDesc)

          call Grid_getBlkPtr(blockDesc, solnData, CENTER)
          call Eos_wrapped(gr_eosMode, blockDesc%limitsGC, solnData)
          call Grid_releaseBlkPtr(blockDesc, solnData, CENTER)

          call itor%next()
       end do
       call Grid_releaseLeafIterator(itor)
     end if

     call Timers_stop("Grid_updateRefinement")

     ! Only log on the first call
     gcMaskArgsLogged = .TRUE.

     ! DEV: TODO What happens with particles here?

     if (present(gridChanged)) then
        gridChanged = gr_amrexDidRefinement
     end if
  else
     if (present(gridChanged)) then
        gridChanged = .FALSE.
     end if
  end if

end subroutine Grid_updateRefinement


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
  use amrex_amrcore_module, ONLY : amrex_regrid, &
                                   amrex_get_finest_level

  use Grid_interface,       ONLY : Grid_fillGuardCells, &
                                   Grid_getBlkPtr, Grid_releaseBlkPtr
  use Grid_data,            ONLY : gr_nrefs, &
!                                   gr_maxRefine, &
                                   gr_refine_var, &
                                   gr_numRefineVars, &
                                   gr_eosMode, &
                                   gr_amrexDidRefinement
  use gr_interface,         ONLY : gr_getBlkIterator, &
                                   gr_releaseBlkIterator
  use gr_amrexInterface,    ONLY : gr_primitiveToConserveLevel, &
                                   gr_conserveToPrimitiveLevel
  use gr_iterator,          ONLY : gr_iterator_t
  use block_metadata,       ONLY : block_metadata_t
  use Eos_interface,        ONLY : Eos_wrapped
  use Timers_interface,     ONLY : Timers_start, Timers_stop
 
  implicit none

  integer, intent(in)            :: nstep
  real,    intent(in)            :: time
  logical, intent(out), OPTIONAL :: gridChanged

  integer, parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS

  logical, save :: gcMaskArgsLogged = .FALSE.

  integer :: finest_level
  logical :: gcMask(maskSize)
  integer :: iref
  integer :: lev
  integer :: i

  type(gr_iterator_t)            :: itor
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

     gcMask = .FALSE.
     do i = 1, gr_numRefineVars
        iref = gr_refine_var(i)
        gcMask(iref) = (iref > 0)
     end do
     gcMask(NUNK_VARS+1:min(maskSize,NUNK_VARS+NDIM*NFACE_VARS)) = .TRUE.

     call Timers_start("Grid_updateRefinement")
     call Grid_fillGuardCells(CENTER, ALLDIR, doEos=.FALSE.)

     ! Run EoS on non-leaf interiors/GC to get best possible 
     ! refinement decisions
     ! DEV: TODO This should be over all ANCESTOR blocks as the leaves were
     ! given to us with EoS run.
     ! DEV: TODO Confirm with AMReX team if non-parent ancestor blocks can
     ! influence refinement decisions.  If no, then we need only apply EoS to
     ! parent blocks.
     call gr_getBlkIterator(itor)
     do while (itor%is_valid())
        call itor%blkMetaData(blockDesc)

        call Grid_getBlkPtr(blockDesc, solnData, CENTER)
        ! DEV: TODO Add masking as in Grid_fillGuardCell?
        call Eos_wrapped(gr_eosMode, blockDesc%limitsGC, solnData)
        call Grid_releaseBlkPtr(blockDesc, solnData, CENTER)

        call itor%next()
     end do
     call gr_releaseBlkIterator(itor)

     ! Regridding requires interpolation when leaf blocks are created
     ! => data must be in conserved form
     finest_level = amrex_get_finest_level()
     do lev = 0, finest_level
       call gr_primitiveToConserveLevel(lev)
     end do

     gr_amrexDidRefinement = .FALSE.
     call amrex_regrid(0, time)

     do lev = 0, finest_level
       call gr_conserveToPrimitiveLevel(lev, .TRUE.)
     end do

     if (gr_amrexDidRefinement) then
       ! We don't know which leaf blocks were created, so run EoS 
       ! on all blocks
       ! DEV: TODO: make gr_amrexDidRefinement into an array so that
       ! we only run EoS on those levels that were really changed.
       ! Could we get away with just running EoS on leaf blocks on 
       ! these levels?
       call gr_getBlkIterator(itor)
       do while (itor%is_valid())
          call itor%blkMetaData(blockDesc)

          call Grid_getBlkPtr(blockDesc, solnData, CENTER)
          call Eos_wrapped(gr_eosMode, blockDesc%limitsGC, solnData)
          call Grid_releaseBlkPtr(blockDesc, solnData, CENTER)

          call itor%next()
       end do
       call gr_releaseBlkIterator(itor)
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


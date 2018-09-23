!!****if* source/Grid/GridMain/AMR/Amrex/Grid_updateRefinement
!!
!! NAME
!!  Grid_updateRefinement
!!
!! SYNOPSIS
!!  call Grid_updateRefinement(integer(IN)  :: nstep,
!!                             real(IN)     :: time,
!!                   optional, logical(OUT) :: gridChanged)
!!
!! DESCRIPTION
!!  If the indicated step qualifies as a refinement step, then this routine
!!    (1) converts all primitive form leaf data to conserved form,
!!    (2) restricts data from leaf blocks down to all ancestors,
!!    (3) fills all guardcells at all levels,
!!    (4) reverts conserved form leaf data to primitive form where 
!!        necessary,
!!    (5) runs EoS on interiors/GC of all blocks,
!!    (6) triggers AMReX to execute grid refinement, and
!!    (7) runs EoS on interiors/GC of all leaf blocks if the mesh was updated.
!!
!!  Note that steps (1) and (4) are skipped if the runtime parameters 
!!  convertToConsvdInMeshInterp and convertToConsvdForMeshCalls indicate that
!!  conversion is not desired.
!!
!!  It is assumed that the data in all leaf block interiors is correct.
!!  and that face variable data is not used to identify blocks for 
!!  refinement/derefinement.  Upon termination, this routine only guarantees
!!  correct data on the interiors and guardcells of leaf blocks.
!!
!!  Note that all EoS runs are done in the mode specified by the eosMode
!!  runtime parameter.  Note also that a step qualifies as a refinement step 
!!  if the given step number is a multiple of the nrefs runtime parameter.
!!
!!  AMReX has FLASH identify blocks requiring refinement via the 
!!  gr_markRefineDerefineCallback routine.
!!
!!  After the refinement, AMReX fills cell-centered data in newly-created child
!!  blocks via callbacks listed below.  Please refer to the documentation of
!!  these for more information on how the data in these blocks is set.
!!
!!  Presently, this routine does not alter any particle information.
!!
!! ARGUMENTS
!!  nstep - current step number
!!  time  - current evolution time
!!  gridChanged - returns TRUE if grid actually have changed.
!!
!! SEE ALSO
!!  Grid_fillGuardCells
!!  gr_markRefineDerefineCallback
!!  gr_makeFineLevelFromCoarseCallback
!!  gr_remakeLevelCallback
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
  
  use Grid_interface,            ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Grid_data,                 ONLY : gr_nrefs, &
!                                        gr_maxRefine, &
                                        gr_refine_var, &
                                        gr_numRefineVars, &
                                        gr_eosMode, &
                                        gr_convertToConsvdInMeshInterp, &
                                        gr_smallrho, &
                                        gr_smalle, &
                                        gr_amrexDidRefinement, &
                                        gr_interpolator, &
                                        lo_bc_amrex, hi_bc_amrex
  use gr_interface,              ONLY : gr_getBlkIterator, &
                                        gr_releaseBlkIterator
  use gr_amrexInterface,         ONLY : gr_conserveToPrimitive, &
                                        gr_cleanDensityData, &
                                        gr_cleanEnergyData, &
                                        gr_fillPhysicalBC, &
                                        gr_restrictAllLevels
  use gr_physicalMultifabs,      ONLY : unk
  use gr_iterator,               ONLY : gr_iterator_t
  use block_metadata,            ONLY : block_metadata_t
  use Eos_interface,             ONLY : Eos_wrapped
  use Timers_interface,          ONLY : Timers_start, Timers_stop
 
  implicit none

  integer, intent(in)            :: nstep
  real,    intent(in)            :: time
  logical, intent(out), OPTIONAL :: gridChanged

!  integer, parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS

!  logical, save :: gcMaskArgsLogged = .FALSE.

!  logical :: gcMask(maskSize)
!  integer :: iref
  integer :: lev
  integer :: i
  
  logical :: needConversion

  type(gr_iterator_t)            :: itor
  type(block_metadata_t)         :: blockDesc
  real,                  pointer :: solnData(:, :, :, :) => null()

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

     !!!!! POPULATE ALL BLOCKS AT ALL LEVELS WITH CONSERVATIVE FORM DATA
     ! Only convert if requested
     needConversion = gr_convertToConsvdInMeshInterp
     
     ! Restrict data from leaves to coarser blocks.  Leave in conservative
     ! form as this is potentially needed for interpolation with fillpatch
     call gr_restrictAllLevels(CENTER, convertPtoC=needConversion, &
                                       convertCtoP=.FALSE.)

     !!!!! POPULATE GUARDCELLS IN ALL BLOCKS
     ! DEV: TODO Confirm with AMReX team if non-parent ancestor blocks can
     ! influence refinement decisions.  If not, we just need to GC fill
     ! on all levels with leaf blocks.
     ! DEV: TODO Should we restrict GC fill to only refinement variables?
     lev = 0
     call amrex_fillpatch(unk(lev), 1.0d0, unk(lev), &
                                    0.0d0, unk(lev), &
                                    amrex_geom(lev), gr_fillPhysicalBC, &
                                    0.0d0, UNK_VARS_BEGIN, &
                                    UNK_VARS_BEGIN, NUNK_VARS)

     do lev=1, amrex_get_finest_level()
        call amrex_fillpatch(unk(lev), 1.0d0, unk(lev-1), &
                                       0.0d0, unk(lev-1), &
                                       amrex_geom(lev-1), gr_fillPhysicalBC, &
                                       1.0e0, unk(lev  ), &
                                       0.0d0, unk(lev  ), &
                                       amrex_geom(lev  ), gr_fillPhysicalBC, &
                                       0.0d0, UNK_VARS_BEGIN, &
                                       UNK_VARS_BEGIN, NUNK_VARS, &
                                       amrex_ref_ratio(lev-1), gr_interpolator, &
                                       lo_bc_amrex, hi_bc_amrex)
     end do

     ! Clean data to account for possible unphysical values caused by
     ! interpolation, revert to primitive form if needed, and
     ! run EoS on all interiors/GC to get best possible refinement
     ! DEV: TODO Confirm with AMReX team if non-parent ancestor blocks can
     ! influence refinement decisions.  If no, then we need only apply EoS to
     ! parent blocks.
     ! DEV: TODO Add masking as in Grid_fillGuardCell?
     if (needConversion)then
       call gr_getBlkIterator(itor)
       do while (itor%is_valid())
          call itor%blkMetaData(blockDesc)
          call Grid_getBlkPtr(blockDesc, solnData, CENTER)

          call gr_cleanDensityData(gr_smallrho, &
                                   blockDesc%limitsGC(LOW,  :), &
                                   blockDesc%limitsGC(HIGH, :), &
                                   solnData, &
                                   blockDesc%limitsGC(LOW,  :), &
                                   blockDesc%limitsGC(HIGH, :), &
                                   NUNK_VARS)
          call gr_conserveToPrimitive(blockDesc%limitsGC(LOW,  :), &
                                      blockDesc%limitsGC(HIGH, :), &
                                      solnData, &
                                      blockDesc%limitsGC(LOW,  :), &
                                      blockDesc%limitsGC(HIGH, :), &
                                      NUNK_VARS, &
                                      UNK_VARS_BEGIN, NUNK_VARS)
          call gr_cleanEnergyData(gr_smalle, &
                                  blockDesc%limitsGC(LOW,  :), &
                                  blockDesc%limitsGC(HIGH, :), &
                                  solnData, &
                                  blockDesc%limitsGC(LOW,  :), &
                                  blockDesc%limitsGC(HIGH, :), &
                                  NUNK_VARS)

          call Eos_wrapped(gr_eosMode, blockDesc%limitsGC, solnData)
 
          call Grid_releaseBlkPtr(blockDesc, solnData, CENTER)
          call itor%next()
       end do
       call gr_releaseBlkIterator(itor)
     else
       call gr_getBlkIterator(itor)
       do while (itor%is_valid())
          call itor%blkMetaData(blockDesc)
          call Grid_getBlkPtr(blockDesc, solnData, CENTER)

          call gr_cleanDensityData(gr_smallrho, &
                                   blockDesc%limitsGC(LOW,  :), &
                                   blockDesc%limitsGC(HIGH, :), &
                                   solnData, &
                                   blockDesc%limitsGC(LOW,  :), &
                                   blockDesc%limitsGC(HIGH, :), &
                                   NUNK_VARS)
          call gr_cleanEnergyData(gr_smalle, &
                                  blockDesc%limitsGC(LOW,  :), &
                                  blockDesc%limitsGC(HIGH, :), &
                                  solnData, &
                                  blockDesc%limitsGC(LOW,  :), &
                                  blockDesc%limitsGC(HIGH, :), &
                                  NUNK_VARS)

          call Eos_wrapped(gr_eosMode, blockDesc%limitsGC, solnData)

          call Grid_releaseBlkPtr(blockDesc, solnData, CENTER)
          call itor%next()
       end do
       call gr_releaseBlkIterator(itor)
     end if

     ! Let AMReX regrid with callbacks
     gr_amrexDidRefinement = .FALSE.
     call amrex_regrid(0, time)

     ! Rerun EoS on interiors+GC of all leaf blocks if callbacks report that
     !   1) a new level was created,
     !   2) an existing level was remade, or
     !   3) a level was removed completely
     ! as there will be new leaf blocks.  We iterate over all leaf blocks
     ! as we do not know which leaf blocks are new.
     !
     ! We run EoS here rather than in the callbacks as these would not
     ! run EoS on blocks that were parents but that are now leaves due
     ! to the removal of a level.
     if (gr_amrexDidRefinement) then
       call gr_getBlkIterator(itor, LEAF, tiling=.FALSE.)
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
!     gcMaskArgsLogged = .TRUE.

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


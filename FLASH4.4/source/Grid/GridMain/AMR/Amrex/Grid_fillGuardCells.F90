!!****if* source/Grid/GridMain/AMR/Amrex/Grid_fillGuardCells
!!
!! NAME
!!  Grid_fillGuardCells
!!
!! SYNOPSIS
!!
!!  call Grid_fillGuardCells(integer(IN) :: gridDataStruct,
!!                           integer(IN) :: idir,
!!                  optional,integer(IN) :: minLayers,
!!                  optional,integer(IN) :: eosMode,
!!                  optional,logical(IN) :: doEos,
!!                  optional,integer(IN) :: maskSize,
!!                  optional,logical(IN) :: mask(maskSize),
!!                  optional,logical(IN) :: makeMaskConsistent,
!!                  optional,integer(IN) :: selectBlockType,
!!                  optional,logical(IN) :: unitReadsMeshDataOnly)
!!
!!
!! DESCRIPTION 
!!  
!!  The argument "gridDataStruct" can take on one of many valid 
!!  values to determine a specific grid data structure on which to apply
!!  the guardcell fill operation. The currently available options are listed with
!!  the arguments. Most users will use CENTER as the option, 
!!  since applications typically use the cell centered grid data, and they want
!!  guardcells to be filled for all the variables.
!!  More specialized applications, such as the unsplit methods, may want to use
!!  other options. 
!!  The user can also choose to fill guard cells either in a single direction,
!!  or all of them. 
!!
!!
!! ARGUMENTS 
!!  
!!
!!  gridDataStruct - integer constant, defined in "constants.h", 
!!                   indicating which grid data structure 
!!                   variable's guardcells to fill.
!!                   UG has 4 data structures for grid variables that
!!                   can have their guardcells filled. 
!!
!!                   unk                all cell centered variables in the grid
!!                   facex,facey,facez  all face centered variables along i,j,k 
!!                                      direction respectively
!!                   
!!                   valid values of gridDataStruct are  
!!                   CENTER             unk only
!!                   WORK               has no meaning in UG
!!                   The face variables are not yet implemented in UG
!!                   FACES              facex,facey, and facez
!!                   FACEX              facex
!!                   FACEY              facey
!!                   FACEZ              facez
!!                   CENTER_FACES     unk,facex,facey,facez
!!
!!  idir - direction of guardcell fill.  User can specify ALLDIR for all (x,y,z)
!!         directions, or if for example the algorithm only does one directional
!!         sweep at a time then time can be saved by filling only the guardcell
!!         direction that is needed.  A user would pass in the constants defined
!!         in constants.h IAXIS, JAXIS or KAXIS to fill guardcells in only one 
!!         direction.        
!!         All layers of guardcells in the given direction(s) are filled.
!!         In the current UG implementation, idir is ignored and a full
!!         guardcell fill in all directions is always performed.
!!
!!         THE REMAINING ARGUMENTS HAVE NO MEANING IN UG
!!
!!  minLayers - number of guardcell layers requested for all directions.
!!
!!   eosMode  - The mode in which eos is to be applied
!!   doEos    - the UG implementation does not act upon this argument
!!   maskSize - the size of the mask array. 
!! 
!!  mask -  It is a one-dimensional logical array 
!!          with indices corresponding to variables in the grid data
!!          structures. If a variable should have its guardcells filled,
!!          the corresponding element in "mask" is true, otherwise it is
!!          false.
!!          The mask is always ignored if the runtime parameter
!!          enableMaskedGCFill is set .FALSE.
!!  
!! makeMaskConsistent - If true when mask is applied, it is made sure that for
!!          all the selected variables in the mask, the ones they are dependent
!!          on are true too. It is also determined whether there is a need to 
!!          apply Eos if doEos argument is true.
!!
!! selectBlockType - IGNORED
!!
!! unitReadsMeshDataOnly - specifies that the unit calling Grid_fillGuardCells
!!                         does not update any internal grid data.  This
!!                         allows us to skip the next guard cell fill because
!!                         the guard cells already contain up to date data.
!!
!! EXAMPLE
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells(CENTER, IAXIS)
!!
!!     This call will fill all guardcells for all cell-centered 
!!     variables in the x direction.
!!     
!! EXAMPLE 2
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells(CENTER_FACES, ALLDIR)
!!     
!!     This call fills guardcells along all directions in both
!!     cell centered and face centered data structures.
!!
!! NOTES
!!
!!   The masking functionality is not yet included in UG
!!  
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_fillGuardCells(gridDataStruct, idir, &
                               minLayers, &
                               eosMode, doEos, &
                               maskSize, mask, makeMaskConsistent, doLogMask, &
                               selectBlockType, &
                               unitReadsMeshDataOnly)
#include "constants.h"
#include "Flash.h"

  use, INTRINSIC :: iso_c_binding
  use amrex_amrcore_module,   ONLY : amrex_get_finest_level
  use amrex_fillpatch_module, ONLY : amrex_fillpatch
  use Grid_data, ONLY : gr_justExchangedGC, gr_eosModeNow, &
                        gr_convertToConsvdInMeshInterp
!  use Grid_data, ONLY : gr_bndOrder, gr_allPeriodic, gr_justExchangedGC, &
!       gr_eosModeNow, gr_convertToConsvdInMeshInterp, gr_maxRefineLevel
  use gr_bcInterface, ONLY : gr_bcApplyToAllBlks
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_guardCells
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

  integer, intent(in)           :: gridDataStruct
  integer, intent(in)           :: idir
  integer, intent(in), optional :: minLayers
  integer, intent(in), optional :: eosMode
  logical, intent(in), optional :: doEos
  integer, intent(in), optional :: maskSize
  logical, intent(in), optional :: mask(:)
  logical, intent(in), optional :: makeMaskConsistent
  logical, intent(in), optional :: doLogMask
  integer, intent(in), optional :: selectBlockType
  logical, intent(in), optional :: unitReadsMeshDataOnly

  integer :: axis
  logical :: isWork = .FALSE.
  integer :: listBlockType = ALL_BLKS
  integer :: finest_level = -1
  integer :: layersArray(MDIM)
  integer :: gcEosMode
  integer :: numBlocks
  logical :: needEos
  integer :: ierr = 0

! TODO: Implement once code finalized
!#ifdef DEBUG_GRID
!  logical:: validDataStructure,validMaskSize
!  validDataStructure = (gridDataStruct==CENTER).or.&
!                       (gridDataStruct==FACES).or.&
!                       (gridDataStruct==FACEX).or.&
!                       (gridDataStruct==FACEY).or.&
!                       (gridDataStruct==FACEZ).or.&
!                       (gridDataStruct==WORK).or.&
!                       (gridDataStruct==CENTER_FACES)
!  if(.not.validDataStructure)then
!     call Driver_abortFlash("GCfill: invalid data structure")
!  end if
!#endif

  if (gridDataStruct /= CENTER .and. gridDataStruct /= CENTER_FACES) then
     !DEV CD.  I am accepting CENTER_FACES for the time being because it
     !is passed by Grid_markRefineDerefine.  I do not support FACE variables
     !yet so CENTER_FACES is just CENTER for now.
     call Driver_abortFlash("[Grid_fillGuardCells]: Non-center not yet coded")
  end if

  ! TODO: Implement code for optional parameters.
  if (present(mask)) then
     call Driver_abortFlash("[Grid_fillGuardCells]: No guard cell masks yet")
  elseif (present(selectBlockType)) then
     call Driver_abortFlash("[Grid_fillGuardCells]: No block selection yet")
  elseif (present(unitReadsMeshDataOnly)) then
     call Driver_abortFlash("[Grid_fillGuardCells]: No support for unitReadsMeshData yet")
  end if

  ! TODO: From paramesh.  Do we need this?
!  !We can skip this guard cell fill if the guard cells are up to date.
!  if (gridDataStruct /= WORK) then
!     skipThisGcellFill = gr_gcellsUpToDate
!  else
!     skipThisGcellFill = .false.
!  end if

  if(present(eosMode)) then
     gcEosMode=eosMode
  else
     gcEosMode=gr_eosModeNow
  end if

  !! If masking is not done then Eos should be applied since it is not known
  !! which variables are of interest
  needEos=.true.
  layersArray = NGUARD
  ! TODO: Paramesh uses ACTIVE_BLKS, which means leaves and their parents.  
  ! By doing this, a parent block is ready for use if its child is removed.
  ! The AMReX iterator is just going over leaves.  Is that acceptable?
  listBlockType = ACTIVE_BLKS

  ! TODO: From paramesh version.  Necessary here?
!  if (ndim<2) gcell_on_fc(KAXIS,:) = .false.
 
  ! GC data could be managed by other processor.
  ! Wait for work on all data structures across full mesh to finish 
  ! before GC filling
!  call Timers_start("guardcell Barrier")
!  call MPI_BARRIER(gr_meshComm, ierr)
!  call Timers_stop("guardcell Barrier")

  call Timers_start("guardcell internal")

  ! Transform center-based data from primitive to conserved form
  ! TODO: Check if this is being done automatically by AMReX
!  if((gridDataStruct==CENTER_FACES).or.(gridDataStruct==CENTER)) then
!      ! TODO: Can we uncomment these?
!!     if (.NOT. skipThisGcellFill) then
!        itor = block_iterator_t(listBlockType) 
!        do while (itor%is_valid())
!           call itor%blkMetaData(block)
!           call gr_primitiveToConserve(block)
!
!           call itor%next()
!        end do
!!     end if
!
!      ! TODO: These are commented out in paramesh version.  Same here?
!      if (gr_convertToConsvdInMeshInterp) then
!        ! TODO: Update interface of gr_sanitizeDataAfterInterp
!        itor = block_iterator_t(listBlockType) 
!        do while (itor%is_valid())
!           call itor%blkMetaData(block)
!           call gr_sanitizeDataAfterInterp(block, 'after gc filling', layersArray)
!
!           call itor%next()
!        end do
!      end if
!  end if

#ifdef DEBUG_GRID
  print*, 'amrex_fillpatch(PE', iopt, guard, layers,')'
#endif
  call Timers_start("amr_guardcell")

!  finest_level = amrex_get_finest_level()
!  do lev=0, finest_level
!      ! Do not use interpolation in time
!      call amrex_fillpatch(mf(lev)%p, 0.0_wp, smf(lev), 0.0_wp, ns, &
!                           0, 0, NUNK_VARS, gr_fillPhysBc)
!  end do
  call Timers_stop("amr_guardcell")

  gr_justExchangedGC = .TRUE.

  if(present(doEos)) then
     if(doEos .AND. needEos) then
        call Timers_start("eos gc")

        ! TODO: Paramesh version disallows CENTER and CENTER_FACES.
        ! Implement same here?
!        itor = block_iterator_t(listBlockType)
!        do while (itor%is_valid())
!            call itor%blkMetaData(block)
! 
!            call Grid_getBlkPtr(block, solnData)
!            ! TODO: Has Eos_guardCells been impolemented to take data ptr?
!            ! TODO: Paramesh version specifies skipSrl=.TRUE.
!            call Eos_guardCells(gcEosMode, solnData, corners=.TRUE., &
!                                layers=layersArray)
!            call Grid_releaseBlkPtr(block, solnData)
!            nullify(solnData)
!
!            call itor%next()
!        end do
!
        call Timers_stop("eos gc")
     end if
  end if

  ! TODO: Paramesh has code to detect if we can skip the next GC fill.
  ! Implement same here?
  ! TODO: Is AMReX managing this logical?
!  gr_gcellsUpToDate = .FALSE.

  call Timers_stop("guardcell internal")
end subroutine Grid_fillGuardCells


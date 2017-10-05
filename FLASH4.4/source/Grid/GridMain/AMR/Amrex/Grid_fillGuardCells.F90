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

#include "constants.h"
#include "Flash.h"

subroutine Grid_fillGuardCells(gridDataStruct, idir, &
                               minLayers, &
                               eosMode, doEos, &
                               maskSize, mask, makeMaskConsistent, doLogMask, &
                               selectBlockType, &
                               unitReadsMeshDataOnly)
  use, INTRINSIC :: iso_c_binding
  use amrex_amrcore_module,      ONLY : amrex_get_finest_level, &
                                        amrex_geom, &
                                        amrex_ref_ratio
  use amrex_fillpatch_module,    ONLY : amrex_fillpatch
  use amrex_bc_types_module,     ONLY : amrex_bc_int_dir
  use amrex_interpolater_module, ONLY : amrex_interp_cell_cons
  
  use Grid_data,                 ONLY : gr_justExchangedGC, &
        gr_eosModeNow
  use gr_physicalMultifabs,      ONLY : unk
  use gr_amrexInterface,         ONLY : gr_fillPhysicalBC, &
                                        gr_averageDownLevels
  use Driver_interface,          ONLY : Driver_abortFlash
  use Timers_interface,          ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use gr_interface, ONLY : gr_setGcFillNLayers
  use Eos_interface, ONLY : Eos_guardCells
  use block_iterator, ONLY : block_iterator_t, destroy_iterator
  use block_metadata, ONLY : block_metadata_t

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

  integer :: guard, gcEosMode
  integer,dimension(MDIM) :: layers, returnLayers
  integer :: listBlockType
  real,dimension(:,:,:,:),pointer::solnData
  type(block_iterator_t) :: itor
  type(block_metadata_t) :: blockDesc

  logical,parameter :: needEos = .TRUE.

  integer :: lo_bc(NDIM, 1)
  integer :: hi_bc(NDIM, 1)

  integer :: lev
  integer :: finest_level

#ifdef DEBUG_GRID
  logical:: validDataStructure

  validDataStructure = (gridDataStruct==CENTER).or.&
                       (gridDataStruct==FACES).or.&
                       (gridDataStruct==FACEX).or.&
                       (gridDataStruct==FACEY).or.&
                       (gridDataStruct==FACEZ).or.&
                       (gridDataStruct==WORK).or.&
                       (gridDataStruct==CENTER_FACES)
  if (.not.validDataStructure) then
     call Driver_abortFlash("[Grid_fillGuardcell] invalid data structure")
  end if
#endif

  finest_level = -1

  ! DEVNOTE: TODO Implement this functionality
#ifdef DEBUG_GRID
  if      (present(minLayers)) then
    call Driver_abortFlash("[Grid_fillGuardCells] minLayers *not* implemented for yet AMReX") 
  end if
  if (present(eosMode)) then
    write(*,*) "eosMode = ", eosMode
    call Driver_abortFlash("[Grid_fillGuardCells] eosMode *not* implemented for yet AMReX") 
  end if
  if (present(doEos)) then
    write(*,*) "doEos = ", doEos
    call Driver_abortFlash("[Grid_fillGuardCells] doEos *not* implemented for yet AMReX") 
  end if
  if (present(maskSize)) then
    call Driver_abortFlash("[Grid_fillGuardCells] maskSize *not* implemented for yet AMReX") 
  end if
  if (present(mask)) then
    call Driver_abortFlash("[Grid_fillGuardCells] mask *not* implemented for yet AMReX") 
  end if
  if (present(makeMaskConsistent)) then
    call Driver_abortFlash("[Grid_fillGuardCells] makeMaskConsistent *not* implemented for yet AMReX") 
  end if
  if (present(doLogMask)) then
    call Driver_abortFlash("[Grid_fillGuardCells] doLogMask *not* implemented for yet AMReX") 
  end if
  if (present(unitReadsMeshDataOnly)) then
    call Driver_abortFlash("[Grid_fillGuardCells] unitReadsMeshDataOnly *not* implemented for yet AMReX") 
  end if
#endif


  if (gridDataStruct /= CENTER .and. gridDataStruct /= CENTER_FACES) then
     !DEV CD.  I am accepting CENTER_FACES for the time being because it
     !is passed by Grid_markRefineDerefine.  I do not support FACE variables
     !yet so CENTER_FACES is just CENTER for now.
     call Driver_abortFlash("[Grid_fillGuardCells]: Non-center not yet coded")
  end if

  ! DEV: Filling by direction is not needed any longer
#ifdef DEBUG_GRID
  if (idir /= ALLDIR) then
    call Driver_abortFlash("[Grid_fillGuardCells] idir must be ALLDIR with AMReX")
  end if
#endif

  if(present(eosMode)) then
     gcEosMode=eosMode
  else
     gcEosMode=gr_eosModeNow
  end if

  ! GC data could be managed by other processor.
  ! Wait for work on all data structures across full mesh to finish 
  ! before GC filling
  ! DEV: TODO Does AMReX handle synchronization?
!  call Timers_start("guardcell Barrier")
!  call MPI_BARRIER(gr_meshComm, ierr)
!  call Timers_stop("guardcell Barrier")

  ! DEV: TODO How to do guardcell fill by direction?
  call Timers_stop("guardcell internal")

  guard = NGUARD
  listBlockType = ACTIVE_BLKS

  !----------------------------------------------------------------
  ! Figure out nlayers arguments to amr_guardcell based on our arguments
  call gr_setGcFillNLayers(layers, idir, guard, minLayers, returnLayers)

  if (present(selectBlockType)) then
     listBlockType = selectBlockType
     select case (selectBlockType)
     case(LEAF)
        ! ok
     case(ACTIVE_BLKS)
        ! ok
     case(ALL_BLKS)
        call Driver_abortFlash('Grid_fillGuardCells: unsupported value ALL_BLKS for selectBlockType!')
#ifdef DEBUG_GRID
     case default
        call Driver_abortFlash("[Grid_fillGuardCells] selectBlockType *not* implemented for yet AMReX") 
#endif
     end select
  end if

  ! Restrict data from leaves to coarser blocks
  call gr_averageDownLevels

  ! DEVNOTE: FIXME Currently fixing BC to periodic here
  ! DEVNOTE: FIXME Currently fixing interpolation mode to cell conserved
  !                linear (AMReX_Interpolater.H)
  ! DEVNOTE: TODO Since we are not using subcycling, should we just use
  !               amrex_fi_fillpatch_two directly?
  lo_bc(:, :) = amrex_bc_int_dir
  hi_bc(:, :) = amrex_bc_int_dir

  ! DEV: Using fill_boundary didn't work on finest levels since the GC outside
  ! the domain were zero (no periodic BC).  AMReX recommended using fillpatch,
  ! which is copying *all* data, including the GC.
  lev = 0
  call Timers_start("amr_guardcell")
  call amrex_fillpatch(unk(lev), 1.0d0, unk(lev), &
                                 0.0d0, unk(lev), &
                       amrex_geom(lev), gr_fillPhysicalBC, &
                       0.0d0, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS)

  finest_level = amrex_get_finest_level()
  do lev=1, finest_level
    call amrex_fillpatch(unk(lev), 1.0d0, unk(lev-1), &
                                   0.0d0, unk(lev-1), &
                         amrex_geom(lev-1), gr_fillPhysicalBC, &
                               1.0e0, unk(lev  ), &
                               0.0d0, unk(lev  ), &
                         amrex_geom(lev  ), gr_fillPhysicalBC, &
                         0.0d0, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, &
                         amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                         lo_bc, hi_bc)
  end do
  call Timers_stop("amr_guardcell")

  gr_justExchangedGC = .TRUE.

  if(present(doEos)) then
     if(doEos.and.needEos) then
        call Timers_start("eos gc")
        itor = block_iterator_t(listBlockType)
        do while (itor%is_valid())
                call itor%blkMetaData(blockDesc)
                
                call Grid_getBlkPtr(blockDesc, solnData)
                call Eos_guardCells(gcEosMode, solnData, corners=.true., &
                                    layers=returnLayers)
                call Grid_releaseBlkPtr(blockDesc, solnData)
                nullify(solnData)

                call itor%next()
        end do
#if defined(__GFORTRAN__) && (__GNUC__ <= 4)
        call destroy_iterator(itor)
#endif
        call Timers_stop("eos gc")
     end if
  end if

  
  call Timers_stop("guardcell internal")

#ifdef DEBUG_GRID
  write(*,'(A,I3)') "[Grid_fillGuardcell] From level 1 to level ", &
                    finest_level+1
#endif

end subroutine Grid_fillGuardCells


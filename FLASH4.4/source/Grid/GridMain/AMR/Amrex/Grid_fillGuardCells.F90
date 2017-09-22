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
  use, INTRINSIC :: iso_c_binding
  use amrex_amrcore_module,   ONLY : amrex_get_finest_level, &
                                     amrex_geom
  
  use Grid_data,              ONLY : gr_justExchangedGC
  use gr_physicalMultifabs,   ONLY : unk
  use Driver_interface,       ONLY : Driver_abortFlash
  use Timers_interface,       ONLY : Timers_start, Timers_stop

  implicit none

#include "constants.h"
#include "Flash.h"

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
  if      (present(minLayers)) then
    call Driver_abortFlash("[Grid_fillGuardCells] minLayers *not* implemented for yet AMReX") 
  else if (present(eosMode)) then
    call Driver_abortFlash("[Grid_fillGuardCells] eosMode *not* implemented for yet AMReX") 
  else if (present(doEos)) then
    call Driver_abortFlash("[Grid_fillGuardCells] doEos *not* implemented for yet AMReX") 
  else if (present(maskSize)) then
    call Driver_abortFlash("[Grid_fillGuardCells] maskSize *not* implemented for yet AMReX") 
  else if (present(mask)) then
    call Driver_abortFlash("[Grid_fillGuardCells] mask *not* implemented for yet AMReX") 
  else if (present(makeMaskConsistent)) then
    call Driver_abortFlash("[Grid_fillGuardCells] makeMaskConsistent *not* implemented for yet AMReX") 
  else if (present(doLogMask)) then
    call Driver_abortFlash("[Grid_fillGuardCells] doLogMask *not* implemented for yet AMReX") 
  else if (present(selectBlockType)) then
    call Driver_abortFlash("[Grid_fillGuardCells] selectBlockType *not* implemented for yet AMReX") 
  else if (present(unitReadsMeshDataOnly)) then
    call Driver_abortFlash("[Grid_fillGuardCells] unitReadsMeshDataOnly *not* implemented for yet AMReX") 
  end if


  if (gridDataStruct /= CENTER .and. gridDataStruct /= CENTER_FACES) then
     !DEV CD.  I am accepting CENTER_FACES for the time being because it
     !is passed by Grid_markRefineDerefine.  I do not support FACE variables
     !yet so CENTER_FACES is just CENTER for now.
     call Driver_abortFlash("[Grid_fillGuardCells]: Non-center not yet coded")
  end if

  ! GC data could be managed by other processor.
  ! Wait for work on all data structures across full mesh to finish 
  ! before GC filling
  ! DEV: TODO Does AMReX handle synchronization?
!  call Timers_start("guardcell Barrier")
!  call MPI_BARRIER(gr_meshComm, ierr)
!  call Timers_stop("guardcell Barrier")

  ! DEV: TODO How to do guardcell fill by direction?
  call Timers_start("amr_guardcell")

  finest_level = amrex_get_finest_level()
  do lev=0, finest_level
    call unk(lev)%fill_boundary(amrex_geom(lev), cross=.FALSE.)
  end do
  
  gr_justExchangedGC = .TRUE.
  write(*,*) "[Grid_fillGuardcell] From level 0 to level ", finest_level

  call Timers_stop("amr_guardcell")
end subroutine Grid_fillGuardCells


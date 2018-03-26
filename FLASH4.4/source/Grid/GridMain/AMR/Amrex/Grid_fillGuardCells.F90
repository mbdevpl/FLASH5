!!****if* source/Grid/GridMain/AMR/Amrex/Grid_fillGuardCells
!!
!! NAME
!!
!!  Grid_fillGuardCells
!!
!! SYNOPSIS
!!
!!  call Grid_fillGuardCells(integer(IN)           :: gridDataStruct,
!!                           integer(IN)           :: idir,
!!                           integer(IN), optional :: minLayers,
!!                           integer(IN), optional :: eosMode,
!!                           logical(IN), optional :: doEos,
!!                           integer(IN), optional :: maskSize,
!!                           logical(IN), optional :: mask(maskSize),
!!                           logical(IN), optional :: makeMaskConsistent,
!!                           integer(IN), optional :: selectBlockType,
!!                           logical(IN), optional :: unitReadsMeshDataOnly)
!!
!! DESCRIPTION 
!!  
!!  Restrict data from leaf blocks down to all ancestors and fill the guardcells
!!  of the physical quantities indicated with the given mask parameters and that
!!  are of the indicated grid data structure type.  The fill can be restricted 
!!  to a given direction.
!!
!!  Note that the fill step might require that AMReX execute prolongation operations
!!  using the AMReX conservative linear interpolation algorithm for guardcells
!!  at fine/coarse boundaries.
!!
!!  It is assumed that the data in all leaf block interiors is correct and that EoS
!!  has been run for these.  No assumptions are made about the quality of
!!  guardcell data nor ancestor blocks.
!!
!!  Note that while the interior data of ancestor blocks are improved through
!!  restriction, EoS has not been run on them.
!!
!!  If EoS is run on the guardcells, it is done on all levels for all affected
!!  variables and all affected blocks.  If doEos is set to true but eosMode is
!!  not given, then EoS is run using the mode given by the eosMode runtime
!!  parameter.
!!
!! ARGUMENTS 
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
  use amrex_fort_module,         ONLY : wp => amrex_real
  use amrex_amrcore_module,      ONLY : amrex_get_finest_level, &
                                        amrex_geom, &
                                        amrex_ref_ratio
  use amrex_fillpatch_module,    ONLY : amrex_fillpatch
  use amrex_interpolater_module, ONLY : amrex_interp_cell_cons
  
  use Grid_interface,            ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Grid_data,                 ONLY : gr_justExchangedGC, &
                                        gr_eosMode, &
                                        gr_enableMaskedGCFill, &
                                        gr_meshMe, gr_meshComm, &
                                        gr_gcellsUpToDate, &
                                        lo_bc_amrex, hi_bc_amrex
  use Eos_interface,             ONLY : Eos_guardCells
  use Driver_interface,          ONLY : Driver_abortFlash
  use Timers_interface,          ONLY : Timers_start, Timers_stop
  use Logfile_interface,         ONLY : Logfile_stampMessage, &
                                        Logfile_stampVarMask, &
                                        Logfile_stamp
  use gr_amrexInterface,         ONLY : gr_fillPhysicalBC, &
                                        gr_averageDownLevels
  use gr_interface,              ONLY : gr_setGcFillNLayers, &
                                        gr_setMasks_gen, &
                                        gr_makeMaskConsistent_gen, &
                                        gr_getBlkIterator, &
                                        gr_releaseBlkIterator
  use gr_physicalMultifabs,      ONLY : unk
  use gr_iterator,               ONLY : gr_iterator_t
  use block_metadata,            ONLY : block_metadata_t

#include "Flash_mpi_implicitNone.fh"

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

  logical,dimension(NUNK_VARS) :: gcell_on_cc
  integer :: guard, gcEosMode
  integer,dimension(MDIM) :: layers, returnLayers
  real,dimension(:,:,:,:),pointer::solnData
  type(gr_iterator_t) :: itor
  type(block_metadata_t) :: blockDesc

  integer :: ierr

  logical :: needEos

  logical, save :: maskWarningDone = .FALSE.
  logical :: skipThisGcellFill, skipNextGcellFill
  character(len=10) :: tagext
  integer :: scompCC, ncompCC, lcompCC

  type(c_ptr) :: lo_bc_ptr(UNK_VARS_BEGIN:UNK_VARS_END)
  type(c_ptr) :: hi_bc_ptr(UNK_VARS_BEGIN:UNK_VARS_END)

  integer :: lev, j
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

  ! DEV: TODO Implement this functionality
#ifdef DEBUG_GRID
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

#ifdef DEBUG_GRID
  ! Filling by direction is not needed any longer
  if (idir /= ALLDIR) then
    call Driver_abortFlash("[Grid_fillGuardCells] idir must be ALLDIR with AMReX")
  end if
#endif

  skipThisGcellFill = .FALSE.   ! for now

  if(present(eosMode)) then
     gcEosMode=eosMode
  else
     gcEosMode=gr_eosMode
  end if

  needEos=.true.

  if (.NOT. gr_enableMaskedGCFill) then

     !! If masking is disabled then a warning is issued and all masking related
     !! processing is skipped

     if (.NOT. maskWarningDone) then
        call Logfile_stampMessage( 'INFO: Grid_fillGuardCells is ignoring masking.')
        if (gr_meshMe==MASTER_PE) print*,    'INFO: Grid_fillGuardCells is ignoring masking.'
        maskWarningDone = .TRUE.
     end if

  else
     
     !! if masking is not explicitly disabled then the presence of a mask allows 
     !! masking to proceed

     if(present(mask))then
        if(present(maskSize)) then

           !! If both mask and masksize are present, apply the mask
           call gr_setMasks_gen(gridDataStruct,maskSize,mask, &
                gcell_on_cc,                                  &
                enableMaskedGCFill=gr_enableMaskedGCFill)
           if(present(makeMaskConsistent))then
              if(makeMaskConsistent) then
                 !! if the caller routine is asking for a consistency check
                 !! then mask may be modified, and also determine if eos needs
                 !! to be applied based upon the mask consistency
                 call gr_makeMaskConsistent_gen(gridDataStruct,gcEosMode,needEos,gcell_on_cc)
              end if
           end if
        else  !! if mask is present without the maskSize, abort
           call Driver_abortFlash("gcfill :: maskSize must be present with mask")
        end if
     end if
  end if

  ! GC data could be managed by other processor.
  ! Wait for work on all data structures across full mesh to finish 
  ! before GC filling
  if (.not. skipThisGcellFill) then
     call Timers_start("guardcell Barrier")
     call MPI_BARRIER(gr_meshComm, ierr)
     call Timers_stop("guardcell Barrier")
  end if

  call Timers_start("guardcell internal")
  !! appropriately mask the data structures to ensure that only the correct data
  !! structure is filled.
  if((gridDataStruct/=CENTER_FACES).and.(gridDataStruct/=CENTER))gcell_on_cc = .false.

  scompCC = UNK_VARS_BEGIN
  ncompCC = NUNK_VARS

  if(present(mask))then
     if(present(maskSize)) then
        if (gr_enableMaskedGCFill) then
            scompCC = maxloc(merge(1.,0.,gcell_on_cc),dim=1) ! maxloc(gcell_on_cc,dim=1)
            lcompCC = UNK_VARS_END + 1 - &
                      maxloc(merge(1.,0.,gcell_on_cc(UNK_VARS_END:UNK_VARS_BEGIN:-1)),dim=1)
            ncompCC = lcompCC - scompCC + 1
            gcell_on_cc(scompCC:lcompCC) = .TRUE.
        end if

        if (present(doLogMask)) then
           if (doLogMask) then
              if (skipThisGcellFill) then
                 tagext = '(skipped)'
              else
                 tagext = ''
              end if
              if (present(doEos)) then
                 if (doEos) then
                    call Logfile_stampVarMask(gcell_on_cc, needEos, '[Grid_fillGuardCells]'//tagext, 'gcSet')
                 else
                    call Logfile_stampVarMask(gcell_on_cc, needEos, '[Grid_fillGuardCells]'//tagext, 'gcSet[no doEos]')
                 end if
              else
                 call Logfile_stampVarMask(gcell_on_cc, needEos, '[Grid_fillGuardCells]'//tagext, 'gcSet[nop doEos]')
              end if
           end if
        end if
     end if
  end if

  guard = NGUARD

  !----------------------------------------------------------------
  ! Figure out nlayers arguments to amr_guardcell based on our arguments
  call gr_setGcFillNLayers(layers, idir, guard, minLayers, returnLayers)

  if (present(selectBlockType)) then
     call Driver_abortFlash("[Grid_fillGuardCells] selectBlockType *not* implemented for yet AMReX") 
  end if

  ! Restrict data from leaves to coarser blocks
  call gr_averageDownLevels

  ! Using fill_boundary didn't work on finest levels since the GC outside
  ! the domain were zero (no periodic BC).  AMReX recommended using fillpatch,
  ! which is copying *all* data, including the GC.
  call Timers_start("amr_guardcell")

  ! GC Fill with GC EoS on coarsest level
  lev = 0
  call amrex_fillpatch(unk(lev), 1.0d0, unk(lev), &
                                 0.0d0, unk(lev), &
                                 amrex_geom(lev), gr_fillPhysicalBC, &
                                 0.0d0, scompCC, scompCC, ncompCC)

  if (present(doEos) .AND. needEos) then
     if (doEos) then
        call Timers_start("eos gc")
        call gr_getBlkIterator(itor, level=lev+1)
        do while (itor%is_valid())
           call itor%blkMetaData(blockDesc)
           
           call Grid_getBlkPtr(blockDesc, solnData)
           call Eos_guardCells(gcEosMode, solnData, corners=.true., &
                               layers=returnLayers)
           call Grid_releaseBlkPtr(blockDesc, solnData)

           call itor%next()
        end do
        call gr_releaseBlkIterator(itor)
        call Timers_stop("eos gc")
     end if
  end if

  ! GC Fill with GC EoS up to finest level
  finest_level = amrex_get_finest_level()
  do lev=1, finest_level
     call amrex_fillpatch(unk(lev), 1.0d0, unk(lev-1), &
                                    0.0d0, unk(lev-1), &
                                    amrex_geom(lev-1), gr_fillPhysicalBC, &
                                    1.0e0, unk(lev  ), &
                                    0.0d0, unk(lev  ), &
                                    amrex_geom(lev  ), gr_fillPhysicalBC, &
                                    0.0d0, scompCC, scompCC, ncompCC, &
                                    amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                    lo_bc_amrex, hi_bc_amrex) 
  
     if (present(doEos) .AND. needEos) then
        if (doEos) then
           call Timers_start("eos gc")
           call gr_getBlkIterator(itor, level=lev+1)
           do while (itor%is_valid())
              call itor%blkMetaData(blockDesc)
              
              call Grid_getBlkPtr(blockDesc, solnData)
              call Eos_guardCells(gcEosMode, solnData, corners=.true., &
                                  layers=returnLayers)
              call Grid_releaseBlkPtr(blockDesc, solnData)

              call itor%next()
           end do
           call gr_releaseBlkIterator(itor)
           call Timers_stop("eos gc")
        end if
     end if
  end do
  call Timers_stop("amr_guardcell")

  gr_justExchangedGC = .TRUE.

  call Logfile_stamp(finest_level+1, &
          '[Grid_fillGuardCells] GC fill/GC EoS up to level ')

  !We now test whether we can skip the next guard cell fill.
  skipNextGcellFill = .false.
  if(present(unitReadsMeshDataOnly)) then
     if (unitReadsMeshDataOnly) then
        if (gr_gcellsUpToDate) then
           !If *all* guard cells were up to date on entry to
           !Grid_fillGuardCells then they will continue to be up to date.
           skipNextGcellFill = .true.
        else
           !Check whether we filled guardcells for all layers, all
           !variables and all active blocks.  This ensures all
           !guard cells are up to date for the next unit.
           if ((gridDataStruct == CENTER_FACES .OR. &
                (gridDataStruct == CENTER .AND. (NFACE_VARS < 1))) &
                .and. idir == ALLDIR) then
              skipNextGcellFill = .true.
              if (present(minLayers)) then
                 if (minval(layers(1:NDIM)) < guard) then
                    skipNextGcellFill = .false.
                 end if
              end if
              if (present(mask)) then
                 if (.not.all(mask .eqv. .true.)) then
                    skipNextGcellFill = .false.
                 end if
              end if
              if (present(selectBlockType)) then
                 if (selectBlockType /= ACTIVE_BLKS) then
                    skipNextGcellFill = .false.
                 end if
              end if
           end if
        end if
     end if
  end if
  gr_gcellsUpToDate = skipNextGcellFill

  call Timers_stop("guardcell internal")

#ifdef DEBUG_GRID
  write(*,'(A,I3)') "[Grid_fillGuardcell] From level 1 to level ", &
                    finest_level+1
#endif

end subroutine Grid_fillGuardCells


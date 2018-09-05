!!****if* source/Grid/GridMain/AMR/Amrex/Grid_fillGuardCells
!!
!! NAME
!!  Grid_fillGuardCells
!!
!! SYNOPSIS
!!  call Grid_fillGuardCells(integer(IN) :: gridDataStruct,
!!                           integer(IN) :: idir,
!!                 optional, integer(IN) :: minLayers,
!!                 optional, integer(IN) :: eosMode,
!!                 optional, logical(IN) :: doEos,
!!                 optional, integer(IN) :: maskSize,
!!                 optional, logical(IN) :: mask(maskSize),
!!                 optional, logical(IN) :: makeMaskConsistent,
!!                 optional, logical(IN) :: doLogMask
!!                 optional, integer(IN) :: selectBlockType,
!!                 optional, logical(IN) :: unitReadsMeshDataOnly)
!!
!! DESCRIPTION 
!!  For all leaf blocks, fill the guardcells of the physical quantities 
!!  in accord with the given mask parameters and that are of the indicated grid
!!  data structure type.
!!
!!  Specifically, this routine
!!    (1) converts all primitive form leaf data to conserved form,
!!    (2) restricts data from leaf blocks down to all ancestors,
!!    (3) fills all guardcells at all levels,
!!    (4) reverts conserved form leaf data to primitive form where 
!!        necessary, and
!!    (5) runs EoS on cell-centered leaf block guardcells if so desired.
!!
!!  Note that steps (1) and (4) are skipped if the runtime parameters 
!!  convertToConsvdInMeshInterp and convertToConsvdForMeshCalls indicate that
!!  conversion is not desired.
!!
!!  The fill step might require prolongation operations, which use the AMReX
!!  conservative linear interpolation algorithm for guardcells at fine/coarse
!!  boundaries.
!!
!! ARGUMENTS 
!!  gridDataStruct - integer constant that indicates which grid data structure 
!!                   variable's guardcells to fill.  Valid values are  
!!                     CENTER             cell-centered data only
!!                     FACEX              X face-centered data only
!!                     FACEY              Y face-centered data only
!!                     FACEZ              Z face-centered data only
!!                     FACES              All face-centered data only
!!                     CENTER_FACES       cell-centered and all face-centered
!!  idir - For AMReX, the only valid value is ALLDIR, which does the fill along
!!         all directions.
!!  minLayers - number of guardcell layers requested for all directions.
!!  eosMode  - The mode in which eos is to be applied.  If not given, then
!!             EoS will be run in the mode given by the eosMode runtime parameter.
!!  doEos    - run EoS on leaf block guardcells
!!  maskSize - the size of the mask array. 
!!  mask - an array whose indices correspond to variables in the grid data
!!          structures.  The guardcell data for each variable will be filled if
!!          the corresponding element in mask is set to .TRUE.  If the runtime
!!          parameter enableMaskedGCFill is .FALSE., then the routine performs
!!          the fill as if all elements in mask are set to .TRUE.  Note
!!          that the presence of this variable requires that maskSize also be
!!          given.
!!  makeMaskConsistent - If .TRUE., then the mask is altered so that the fill
!!          is also done for all quantities on which unmasked quanitities
!!          depend.  If doEos is .TRUE., it is also determined if the contents
!!          of mask require that EoS be run.  Note that the presence of
!!          this variable requires that mask also be given.
!!  doLogMask - log masking information if given.
!!  selectBlockType - it is an error to give this parameter.
!!  unitReadsMeshDataOnly - it is an error to give this parameter.
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
  use amrex_interpolater_module, ONLY : amrex_interp_cell_cons
  
  use Grid_interface,            ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                                        Grid_getLeafIterator, &
                                        Grid_releaseLeafIterator
  use Grid_data,                 ONLY : gr_justExchangedGC, &
                                        gr_eosMode, &
                                        gr_enableMaskedGCFill, &
                                        gr_convertToConsvdInMeshInterp, &
                                        gr_smallrho, &
                                        gr_smalle, &
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
                                        gr_fillPhysicalFaceBC, &
                                        gr_restrictAllLevels, &
                                        gr_conserveToPrimitive, &
                                        gr_cleanDensityData, &
                                        gr_cleanEnergyData
  use gr_interface,              ONLY : gr_setGcFillNLayers, &
                                        gr_setMasks_gen, &
                                        gr_makeMaskConsistent_gen
  use gr_physicalMultifabs,      ONLY : unk, &
                                        facevarx, facevary, facevarz
  use leaf_iterator,             ONLY : leaf_iterator_t
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
  real,pointer :: solnData(:,:,:,:) => null()
  type(leaf_iterator_t)  :: itor
  type(block_metadata_t) :: blockDesc

  integer :: ierr

  logical :: needEos
  logical :: needConversion

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
                       (gridDataStruct==CENTER_FACES)
  if (.not.validDataStructure) then
     call Driver_abortFlash("[Grid_fillGuardcell] invalid data structure")
  end if
#endif

  ! DEV: TODO Implement this functionality?
  if (       (gridDataStruct /= CENTER) .AND. (gridDataStruct /= CENTER_FACES) &
       .AND. (gridDataStruct /= FACES)  .AND. (gridDataStruct /= FACEX) &
       .AND. (gridDataStruct /= FACEY)  .AND. (gridDataStruct /= FACEZ)) then
     write(*,*) "Unsupported gridDataStruct ", gridDataStruct 
     call Driver_abortFlash("[Grid_fillGuardCells]: Unsupported gridDataStruct")
  else if (idir /= ALLDIR) then
     call Driver_abortFlash("[Grid_fillGuardCells] idir must be ALLDIR with AMReX")
  else if (present(selectBlockType)) then
     call Driver_abortFlash("[Grid_fillGuardCells] selectBlockType *not* implemented for AMReX yet") 
  else if (present(unitReadsMeshDataOnly)) then
     call Driver_abortFlash("[Grid_fillGuardCells] unitReadsMeshDataOnly *not* implemented for AMReX yet") 
  end if

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
!  if((gridDataStruct/=CENTER_FACES).and.(gridDataStruct/=CENTER))gcell_on_cc = .false.

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

  !!!!! POPULATE ALL BLOCKS AT ALL LEVELS WITH CONSERVATIVE FORM DATA
  ! Only convert if requested
  needConversion = gr_convertToConsvdInMeshInterp

  ! Restrict data from leaves to coarser blocks.  Leave in conservative
  ! form as this is potentially needed for interpolation with fillpatch
  call gr_restrictAllLevels(gridDataStruct, convertPtoC=needConversion, &
                                            convertCtoP=.FALSE.)

  !!!!!----- FILL GUARDCELLS ON ALL BLOCKS, ALL LEVELS
  call Timers_start("amr_guardcell")

  !!!!! Cell-centered data first
  if ((gridDataStruct == CENTER) .OR. (gridDataStruct == CENTER_FACES)) then
    lev = 0
    ! AMReX recommended using fillpatch, which is copying *all* data, 
    ! including the GC.
    call amrex_fillpatch(unk(lev), 1.0, unk(lev), &
                                   0.0, unk(lev), &
                                   amrex_geom(lev), gr_fillPhysicalBC, &
                                   0.0, scompCC, scompCC, ncompCC)

    finest_level = amrex_get_finest_level()
    do lev=1, finest_level
       call amrex_fillpatch(unk(lev), 1.0, unk(lev-1), &
                                      0.0, unk(lev-1), &
                                      amrex_geom(lev-1), gr_fillPhysicalBC, &
                                      1.0, unk(lev  ), &
                                      0.0, unk(lev  ), &
                                      amrex_geom(lev  ), gr_fillPhysicalBC, &
                                      0.0, scompCC, scompCC, ncompCC, &
                                      amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                      lo_bc_amrex, hi_bc_amrex)
    end do

    !!!!! FINALIZE CELL-CENTERED DATA
    ! Clean data to account for possible unphysical values caused by
    ! interpolation, revert to primitive form if needed, and
    ! run EoS if needed
    call Timers_start("eos gc")

    if (present(doEos)) then
       needEos = (needEos .AND. doEos)
    else
       needEos = .FALSE.
    end if

    call Grid_getLeafIterator(itor, tiling=.FALSE.)
    if (needEos .AND. needConversion) then
       do while (itor%is_valid())
          call itor%blkMetaData(blockDesc)
          call Grid_getBlkPtr(blockDesc, solnData)

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

          call Eos_guardCells(gcEosMode, solnData, corners=.true., &
                              layers=returnLayers)

          call Grid_releaseBlkPtr(blockDesc, solnData)
          call itor%next()
       end do
    else if (needEos .AND. (.NOT. needConversion)) then
       do while (itor%is_valid())
          call itor%blkMetaData(blockDesc)
          call Grid_getBlkPtr(blockDesc, solnData)

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

          call Eos_guardCells(gcEosMode, solnData, corners=.true., &
                              layers=returnLayers)

          call Grid_releaseBlkPtr(blockDesc, solnData)
          call itor%next()
       end do
    else if (needConversion) then
       do while (itor%is_valid())
          call itor%blkMetaData(blockDesc)
          call Grid_getBlkPtr(blockDesc, solnData)

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

          call Grid_releaseBlkPtr(blockDesc, solnData)
          call itor%next()
       end do
    else
       do while (itor%is_valid())
          call itor%blkMetaData(blockDesc)
          call Grid_getBlkPtr(blockDesc, solnData)

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

          call Grid_releaseBlkPtr(blockDesc, solnData)
          call itor%next()
       end do
    end if
    call Grid_releaseLeafIterator(itor)
  end if   ! End CENTER or CENTER_FACES

#if NFACE_VARS > 0
  !!!!!----- FILL FACEVAR[XYZ] GUARDCELLS
  ! Fill FACEVARX GC if it exists and is so desired
  ! DEV: TODO Do we need C-to-P conversion here for face vars?
  if (     (gridDataStruct == CENTER_FACES) &
      .OR. (gridDataStruct == FACES) .OR. (gridDataStruct == FACEX)) then
     lev = 0
     call amrex_fillpatch(facevarx(lev), 1.0, facevarx(lev), &
                                         0.0, facevarx(lev), &
                                         amrex_geom(lev), gr_fillPhysicalFaceBC, &
                                         0.0, 1, 1, NFACE_VARS)

     do lev=1, amrex_get_finest_level()
        call amrex_fillpatch(facevarx(lev), 1.0, facevarx(lev-1), &
                                            0.0, facevarx(lev-1), &
                                            amrex_geom(lev-1), gr_fillPhysicalFaceBC, &
                                            1.0, facevarx(lev  ), &
                                            0.0, facevarx(lev  ), &
                                            amrex_geom(lev  ), gr_fillPhysicalFaceBC, &
                                            0.0, 1, 1, NFACE_VARS, &
                                            amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                            lo_bc_amrex, hi_bc_amrex) 
     end do
  end if
#if NDIM >= 2
  ! Fill FACEVARY GC if it exists and is so desired
  if (     (gridDataStruct == CENTER_FACES) &
      .OR. (gridDataStruct == FACES) .OR. (gridDataStruct == FACEY)) then
     lev = 0
     call amrex_fillpatch(facevary(lev), 1.0, facevary(lev), &
                                         0.0, facevary(lev), &
                                         amrex_geom(lev), gr_fillPhysicalFaceBC, &
                                         0.0, 1, 1, NFACE_VARS)

     do lev=1, amrex_get_finest_level()
        call amrex_fillpatch(facevary(lev), 1.0, facevary(lev-1), &
                                            0.0, facevary(lev-1), &
                                            amrex_geom(lev-1), gr_fillPhysicalFaceBC, &
                                            1.0, facevary(lev  ), &
                                            0.0, facevary(lev  ), &
                                            amrex_geom(lev  ), gr_fillPhysicalFaceBC, &
                                            0.0, 1, 1, NFACE_VARS, &
                                            amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                            lo_bc_amrex, hi_bc_amrex) 
     end do
  end if
#endif
#if NDIM == 3
  ! Fill FACEVARZ GC if it exists and is so desired
  if (     (gridDataStruct == CENTER_FACES) &
      .OR. (gridDataStruct == FACES) .OR. (gridDataStruct == FACEZ)) then
     lev = 0
     call amrex_fillpatch(facevarz(lev), 1.0, facevarz(lev), &
                                         0.0, facevarz(lev), &
                                         amrex_geom(lev), gr_fillPhysicalFaceBC, &
                                         0.0, 1, 1, NFACE_VARS)

     do lev=1, amrex_get_finest_level()
        call amrex_fillpatch(facevarz(lev), 1.0, facevarz(lev-1), &
                                            0.0, facevarz(lev-1), &
                                            amrex_geom(lev-1), gr_fillPhysicalFaceBC, &
                                            1.0, facevarz(lev  ), &
                                            0.0, facevarz(lev  ), &
                                            amrex_geom(lev  ), gr_fillPhysicalFaceBC, &
                                            0.0, 1, 1, NFACE_VARS, &
                                            amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
                                            lo_bc_amrex, hi_bc_amrex) 
     end do
  end if
#endif
#else
  if (     (gridDataStruct == FACES) .OR. (gridDataStruct == FACEX) &
      .OR. (gridDataStruct == FACEY) .OR. (gridDataStruct == FACEZ)) then
    call Driver_abortFlash("[Grid_fillGuardCells] No face data to work with")
  end if
#endif

  call Timers_stop("eos gc")
  call Timers_stop("amr_guardcell")

  gr_justExchangedGC = .TRUE.

  call Logfile_stamp(finest_level+1, &
          '[Grid_fillGuardCells] GC fill/GC EoS up to level ')

  !We now test whether we can skip the next guard cell fill.
  skipNextGcellFill = .false.
!  if(present(unitReadsMeshDataOnly)) then
!     if (unitReadsMeshDataOnly) then
!        if (gr_gcellsUpToDate) then
!           !If *all* guard cells were up to date on entry to
!           !Grid_fillGuardCells then they will continue to be up to date.
!           skipNextGcellFill = .true.
!        else
!           !Check whether we filled guardcells for all layers, all
!           !variables and all active blocks.  This ensures all
!           !guard cells are up to date for the next unit.
!           if ((gridDataStruct == CENTER_FACES .OR. &
!                (gridDataStruct == CENTER .AND. (NFACE_VARS < 1))) &
!                .and. idir == ALLDIR) then
!              skipNextGcellFill = .true.
!              if (present(minLayers)) then
!                 if (minval(layers(1:NDIM)) < guard) then
!                    skipNextGcellFill = .false.
!                 end if
!              end if
!              if (present(mask)) then
!                 if (.not.all(mask .eqv. .true.)) then
!                    skipNextGcellFill = .false.
!                 end if
!              end if
!              if (present(selectBlockType)) then
!                 if (selectBlockType /= ACTIVE_BLKS) then
!                    skipNextGcellFill = .false.
!                 end if
!              end if
!           end if
!        end if
!     end if
!  end if
  gr_gcellsUpToDate = skipNextGcellFill

  call Timers_stop("guardcell internal")

#ifdef DEBUG_GRID
  write(*,'(A,I3)') "[Grid_fillGuardcell] From level 1 to level ", &
                    finest_level+1
#endif

end subroutine Grid_fillGuardCells

!!****if* source/Grid/GridMain/AMR/Amrex/gr_restrictAllLevels
!!
!! NAME
!!  gr_restrictAllLevels
!!
!! SYNOPSIS
!!  call gr_restrictAllLevels(integer(IN) :: gridDataStruct, 
!!                            logical(IN) :: convertPtoC,
!!                            logical(IN) :: convertCtoP)
!!
!!  For each leaf block, average the leaf data associated with the given grid
!!  data structure type down to all anscestor blocks.
!!
!!  For performance reasons, this routine does not make any assumptions about
!!  whether the data is presently in primitive or conservative form.  Nor does
!!  it make assumptions about which form it should be left.
!!
!! DESCRIPTION 
!!
!! ARGUMENTS 
!!  gridDataStruct - integer constant that indicates which grid data structure 
!!                   variable's data to average.  Valid values are
!!                     CENTER             cell-centered data only
!!                     FACEX              X face-centered data only
!!                     FACEY              Y face-centered data only
!!                     FACEZ              Z face-centered data only
!!                     FACES              All face-centered data only
!!                     CENTER_FACES       cell-centered and all face-centered
!! convertPtoC - if this value is true, then all primitive form quantities on 
!!               all leaf blocks will be converted to conservative form before 
!!               averaging.
!! convertCtoP - if this value is true, then all primitive form quantities will
!!               be reverted back to primitive form after averaging.
!!  
!!***

#include "Flash.h"
#include "constants.h"

subroutine gr_restrictAllLevels(gridDataStruct, convertPtoC, convertCtoP)
  use amrex_amrcore_module,      ONLY : amrex_get_finest_level, &
                                        amrex_geom, &
                                        amrex_ref_ratio
  use amrex_multifabutil_module, ONLY : amrex_average_down

  use Grid_interface,            ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use gr_interface,              ONLY : gr_getBlkIterator, &
                                        gr_releaseBlkIterator
  use gr_amrexInterface,         ONLY : gr_primitiveToConserve, &
                                        gr_conserveToPrimitive
  use gr_iterator,               ONLY : gr_iterator_t
  use block_metadata,            ONLY : block_metadata_t
  use gr_physicalMultifabs,      ONLY : unk, &
                                        facevarx, facevary, facevarz
  use Driver_interface,          ONLY : Driver_abortFlash

  implicit none

  integer, intent(IN) :: gridDataStruct
  logical, intent(IN) :: convertPtoC
  logical, intent(IN) :: convertCtoP

  integer :: lev
  integer :: finest_level

  type(gr_iterator_t)    :: itor
  type(block_metadata_t) :: blockDesc

  real,   pointer :: solnData(:,:,:,:) => null()

  if (       (gridDataStruct /= CENTER) .AND. (gridDataStruct /= CENTER_FACES) &
       .AND. (gridDataStruct /= FACES)  .AND. (gridDataStruct /= FACEX) &
       .AND. (gridDataStruct /= FACEY)  .AND. (gridDataStruct /= FACEZ)) then
     write(*,*) "Unsupported gridDataStruct ", gridDataStruct 
     call Driver_abortFlash("[gr_restrictAllLevels]: Unsupported gridDataStruct")
  end if

  ! Work in AMReX 0-based level indexing
  finest_level = amrex_get_finest_level()

  !!!!! CELL-CENTERED DATA
  if ((gridDataStruct == CENTER) .OR. (gridDataStruct == CENTER_FACES)) then

    ! Convert primitive form to conservative form on leaves only as
    ! averaging will propagate conservative form down to ancestors
    if (convertPtoC) then
      call gr_getBlkIterator(itor, LEAF)
      do while (itor%is_valid())
        call itor%blkMetaData(blockDesc)
        call Grid_getBlkPtr(blockDesc, solnData, CENTER)
        
        call gr_primitiveToConserve(blockDesc%limitsGC(LOW,  :), &
                                    blockDesc%limitsGC(HIGH, :), &
                                    solnData, &
                                    blockDesc%limitsGC(LOW,  :), &
                                    blockDesc%limitsGC(HIGH, :), &
                                    NUNK_VARS, &
                                    UNK_VARS_BEGIN, NUNK_VARS)

        call Grid_releaseBlkPtr(blockDesc, solnData, CENTER)
        call itor%next()
      end do
      call gr_releaseBlkIterator(itor)
    end if

    ! Average from finest down to coarsest
    do lev = finest_level, 1, -1
        call amrex_average_down(unk(lev  ), &
                                unk(lev-1), &
                                amrex_geom(lev  ), &
                                amrex_geom(lev-1), &
                                UNK_VARS_BEGIN, NUNK_VARS, &
                                amrex_ref_ratio(lev-1))
    end do 

    ! Revert conservative form back to primitive form on all blocks
    if (convertCtoP) then
      call gr_getBlkIterator(itor)
      do while (itor%is_valid())
        call itor%blkMetaData(blockDesc)
        call Grid_getBlkPtr(blockDesc, solnData, CENTER)

        call gr_conserveToPrimitive(blockDesc%limitsGC(LOW,  :), &
                                    blockDesc%limitsGC(HIGH, :), &
                                    solnData, &
                                    blockDesc%limitsGC(LOW,  :), &
                                    blockDesc%limitsGC(HIGH, :), &
                                    NUNK_VARS, &
                                    UNK_VARS_BEGIN, NUNK_VARS)

        call Grid_releaseBlkPtr(blockDesc, solnData, CENTER)
        call itor%next()
      end do
      call gr_releaseBlkIterator(itor)
    end if

  end if

#if NFACE_VARS > 0
  !!!!! FACE-CENTERED DATA
  if (     (gridDataStruct == CENTER_FACES) &
      .OR. (gridDataStruct == FACES) .OR. (gridDataStruct == FACEX)) then
    do lev = finest_level, 1, -1
        call amrex_average_down(facevarx(lev  ), &
                                facevarx(lev-1), &
                                amrex_geom(lev  ), &
                                amrex_geom(lev-1), &
                                1, NFACE_VARS, &
                                amrex_ref_ratio(lev-1))
    end do 
  end if
#if NDIM >= 2
  if (     (gridDataStruct == CENTER_FACES) &
      .OR. (gridDataStruct == FACES) .OR. (gridDataStruct == FACEY)) then
    do lev = finest_level, 1, -1
        call amrex_average_down(facevary(lev  ), &
                                facevary(lev-1), &
                                amrex_geom(lev  ), &
                                amrex_geom(lev-1), &
                                1, NFACE_VARS, &
                                amrex_ref_ratio(lev-1))
    end do 
  end if
#endif
#if NDIM == 3
  if (     (gridDataStruct == CENTER_FACES) &
      .OR. (gridDataStruct == FACES) .OR. (gridDataStruct == FACEZ)) then
    do lev = finest_level, 1, -1
        call amrex_average_down(facevarz(lev  ), &
                                facevarz(lev-1), &
                                amrex_geom(lev  ), &
                                amrex_geom(lev-1), &
                                1, NFACE_VARS, &
                                amrex_ref_ratio(lev-1))
    end do 
  end if
#endif
#else
  if (     (gridDataStruct == FACES) .OR. (gridDataStruct == FACEX) &
      .OR. (gridDataStruct == FACEY) .OR. (gridDataStruct == FACEZ)) then
    call Driver_abortFlash("[gr_restrictAllLevels] No face data to work with")
  end if
#endif

end subroutine gr_restrictAllLevels


!!****if* source/Grid/GridMain/AMR/Amrex/gr_averageDownLevels
!!
!! NAME
!!  gr_averageDownLevels
!!
!! SYNOPSIS
!!  call gr_averageDownLevels(integer(IN) :: gridDataStruct)
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
!!  
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "Flash.h"

subroutine gr_averageDownLevels(gridDataStruct)
    use amrex_amrcore_module,      ONLY : amrex_get_finest_level, &
                                          amrex_geom, &
                                          amrex_ref_ratio
    use amrex_multifabutil_module, ONLY : amrex_average_down

    use gr_physicalMultifabs,      ONLY : unk, &
                                          facevarx, facevary, facevarz
    use Driver_interface,          ONLY : Driver_abortFlash

    implicit none

    integer, intent(IN) :: gridDataStruct

    integer :: lev
    integer :: finest_level

  if (       (gridDataStruct /= CENTER) .AND. (gridDataStruct /= CENTER_FACES) &
       .AND. (gridDataStruct /= FACES)  .AND. (gridDataStruct /= FACEX) &
       .AND. (gridDataStruct /= FACEY)  .AND. (gridDataStruct /= FACEZ)) then
     write(*,*) "Unsupported gridDataStruct ", gridDataStruct 
     call Driver_abortFlash("[gr_averageDownLevels]: Unsupported gridDataStruct")
  end if

    ! Work in AMReX 0-based level indexing
    finest_level = amrex_get_finest_level()

  !!!!! CELL-CENTERED DATA
  if ((gridDataStruct == CENTER) .OR. (gridDataStruct == CENTER_FACES)) then
    do lev = finest_level, 1, -1
#ifdef DEBUG_GRID
        write(*,'(A,A,I2,A,I2)') "[gr_averageDownLevels]", &
                                 "               Cell-centered from ", &
                                 lev+1, " down to ", lev
#endif
        call amrex_average_down(unk(lev  ), &
                                unk(lev-1), &
                                amrex_geom(lev  ), &
                                amrex_geom(lev-1), &
                                UNK_VARS_BEGIN, NUNK_VARS, &
                                amrex_ref_ratio(lev-1))
    end do 
  end if

#if NFACE_VARS > 0
  !!!!! FACE-CENTERED DATA
  if (     (gridDataStruct == CENTER_FACES) &
      .OR. (gridDataStruct == FACES) .OR. (gridDataStruct == FACEX) then
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
      .OR. (gridDataStruct == FACES) .OR. (gridDataStruct == FACEY) then
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
      .OR. (gridDataStruct == FACES) .OR. (gridDataStruct == FACEZ) then
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
    call Driver_abortFlash("[gr_averageDownLevels] No face data to work with")
  end if
#endif

end subroutine gr_averageDownLevels


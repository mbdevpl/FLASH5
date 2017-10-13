!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestInit/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!  This code tests that the subroutine Grid_fillGuardCells is correctly filling
!!  guardcells for all cell-centered data variables at a single level.
!!
!! NOTES
!!  This simulation *must* be configured with at least the following
!!  2D run:
!!     ./setup -auto -2d -nxb=2 -nyb=2 
!!              unitTest/Grid/Amrex/TestGcFill
!!             +noio -index-reorder
!!
!!  For the future:
!!             -unit=IO/IOMain/hdf5/serial/AM
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Driver_evolveFlash()
    use amrex_fort_module,     ONLY : amrex_spacedim
    use amrex_box_module,      ONLY : amrex_box
    use amrex_parallel_module, ONLY : amrex_parallel_myproc, &
                                      amrex_parallel_nprocs
    use amrex_multifab_module, ONLY : amrex_mfiter, &
                                      amrex_mfiter_build, &
                                      amrex_mfiter_destroy
    
    use Grid_interface,        ONLY : Grid_fillGuardCells, &
                                      Grid_getBlkPtr, Grid_releaseBlkPtr, &
                                      Grid_getSingleCellCoords
    use gr_amrexInterface,     ONLY : gr_getFinestLevel
    use gr_physicalMultifabs,  ONLY : unk
    use block_metadata,        ONLY : block_metadata_t
    use ut_testDriverMod

    implicit none

    integer, parameter :: FINEST_LEVEL_EX = 2

    integer :: finest_level
    integer :: lev
    integer :: i, j, k
    integer :: var

    type(amrex_mfiter)     :: mfi
    type(amrex_box)        :: bx
    type(block_metadata_t) :: blockDesc

    real, contiguous, pointer :: initData(:, :, :, :)

    integer :: idx(1:MDIM)
    real    :: coords(1:MDIM)

    !!!!! CONFIRM PROPER COORDINATE SYSTEM DIMENSIONALITY
    write(*,*)
    if (amrex_spacedim /= NDIM) then
        write(*,*) "Wrong dimensionality - ", amrex_spacedim, ' != ', NDIM
        write(*,*) "Recompile AMReX with correct dimensionality"
        write(*,*)
        stop
    end if

    call start_test_run

    !!!!! CONFIRM THAT GC ARE ZERO AFTER INIT 
    call gr_getFinestLevel(finest_level)
    call assertEqual(finest_level, FINEST_LEVEL_EX, "Incorrect finest level")

    do lev = 1, finest_level
        call amrex_mfiter_build(mfi, unk(lev-1), tiling=.FALSE.)
        do while (mfi%next())
            bx = mfi%tilebox()

            ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
            ! Level must be 1-based index and limits/limitsGC must be 1-based also
            blockDesc%level = lev
            blockDesc%grid_index = mfi%grid_index()
            blockDesc%limits(LOW,  :) = 1
            blockDesc%limits(HIGH, :) = 1
            blockDesc%limits(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1
            blockDesc%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1
            blockDesc%limitsGC(LOW,  :) = 1
            blockDesc%limitsGC(HIGH, :) = 1
            blockDesc%limitsGC(LOW,  1:NDIM) = blockDesc%limits(LOW,  1:NDIM) - NGUARD
            blockDesc%limitsGC(HIGH, 1:NDIM) = blockDesc%limits(HIGH, 1:NDIM) + NGUARD

            call Grid_getBlkPtr(blockDesc, initData)
            
            associate(lo   => blockDesc%limits(LOW,  :), &
                      hi   => blockDesc%limits(HIGH, :), &
                      loGC => blockDesc%limitsGC(LOW, :), &
                      hiGC => blockDesc%limitsGC(HIGH, :))
                do     k = loGC(KAXIS), hiGC(KAXIS)
                  do   j = loGC(JAXIS), hiGC(JAXIS)
                    do i = loGC(IAXIS), hiGC(IAXIS)
                      do var=UNK_VARS_BEGIN, UNK_VARS_END
                        if (      (lo(IAXIS) <= i) .AND. (i <= hi(IAXIS)) &
                            .AND. (lo(JAXIS) <= j) .AND. (j <= hi(JAXIS)) &
                            .AND. (lo(KAXIS) <= k) .AND. (k <= hi(KAXIS))) then
                            ! Interior has data
                            idx = [i - lo(IAXIS) + 1, &
                                   j - lo(JAXIS) + 1, &
                                   k - lo(KAXIS) + 1]
                            call Grid_getSingleCellCoords(idx, blockDesc, &
                                                  CENTER, INTERIOR, coords)
                            call assertEqual(initData(i, j, k, var), &
                                             DBLE((coords(IAXIS)+coords(JAXIS)) * var), &
                                             "Bad data")
                        else
                            ! Guardcells not populated yet 
                            call assertEqual(initData(i, j, k, var), 0.0d0, &
                                             "GC not zero")
                        end if
                      end do
                    end do
                  end do
                end do
            end associate

            call Grid_releaseBlkPtr(blockDesc, initData)
        end do
    end do

    call Grid_fillGuardCells(CENTER, ALLDIR)
 
    !!!!! CONFIRM THAT GC ARE NOW FILLED APPROPRIATELY
    call gr_getFinestLevel(finest_level)
    call assertEqual(finest_level, FINEST_LEVEL_EX, "Incorrect finest level")

    do lev = 1, finest_level
        call amrex_mfiter_build(mfi, unk(lev-1), tiling=.FALSE.)
        do while (mfi%next())
            bx = mfi%tilebox()

            ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
            ! Level must be 1-based index and limits/limitsGC must be 1-based also
            blockDesc%level = lev
            blockDesc%grid_index = mfi%grid_index()
            blockDesc%limits(LOW,  :) = 1
            blockDesc%limits(HIGH, :) = 1
            blockDesc%limits(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1
            blockDesc%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1
            blockDesc%limitsGC(LOW,  :) = 1
            blockDesc%limitsGC(HIGH, :) = 1
            blockDesc%limitsGC(LOW,  1:NDIM) = blockDesc%limits(LOW,  1:NDIM) - NGUARD
            blockDesc%limitsGC(HIGH, 1:NDIM) = blockDesc%limits(HIGH, 1:NDIM) + NGUARD

            associate(lo   => blockDesc%limits(LOW,  :), &
                      hi   => blockDesc%limits(HIGH, :), &
                      loGC => blockDesc%limitsGC(LOW, :), &
                      hiGC => blockDesc%limitsGC(HIGH, :))
                
                call Grid_getBlkPtr(blockDesc, initData)

                do     k = loGC(KAXIS), hiGC(KAXIS)
                  do   j = loGC(JAXIS), hiGC(JAXIS)
                    do i = loGC(IAXIS), hiGC(IAXIS)
                      idx = [i - lo(IAXIS) + 1, &
                             j - lo(JAXIS) + 1, &
                             k - lo(KAXIS) + 1]
                      call Grid_getSingleCellCoords(idx, blockDesc, &
                                                    CENTER, INTERIOR, coords)

                      do var=UNK_VARS_BEGIN, UNK_VARS_END
                        if (      (lo(IAXIS) <= i) .AND. (i <= hi(IAXIS)) &
                            .AND. (lo(JAXIS) <= j) .AND. (j <= hi(JAXIS)) &
                            .AND. (lo(KAXIS) <= k) .AND. (k <= hi(KAXIS))) then
                            ! Interior has data
                            call assertEqual(initData(i, j, k, var), &
                                 DBLE((coords(IAXIS) + coords(JAXIS)) * var), &
                                 "Bad data")
                        else
                          ! Spot check a few values for now
                          if (lev == 1) then
                            if (     ((i == 6) .AND. (j == 6)) &
                                .OR. ((i == 6) .AND. (j == 2))) then
                              call assertEqual(initData(i, j, k, var), &
                                               1.5d0*var, &
                                               "Incorrect coarse GC value")
                            else if (     ((i == 5) .AND. (j == 1)) &
                                     .OR. ((i == 5) .AND. (j == 5))) then
                              call assertEqual(initData(i, j, k, var), &
                                               0.5d0*var, &
                                               "Incorrect coarse GC value")
                            else if ((i == 0) .AND. (j == 5)) then
                              call assertEqual(initData(i, j, k, var), &
                                               2.0d0*var, &
                                               "Incorrect coarse GC value")
                            end if
                          else if (lev == 2) then
                            if((i == -1) .AND. (j == -1)) then
                              call assertEqual(initData(i, j, k, var), &
                                               3.25d0*var, &
                                               "Incorrect fine GC value")
                            else if ((i == 0) .AND. (j == 5)) then
                              call assertEqual(initData(i, j, k, var), &
                                               3.0d0*var, &
                                               "Incorrect fine GC value")
                            else if ((i == 0) .AND. (j == 6)) then
                              call assertEqual(initData(i, j, k, var), &
                                               3.25d0*var, &
                                               "Incorrect fine GC value")
                            else if ((i == 3) .AND. (j == 9)) then
                              call assertEqual(initData(i, j, k, var), &
                                               0.75d0*var, &
                                               "Incorrect fine GC value")
                            else if ((i == 3) .AND. (j == 10)) then
                              call assertEqual(initData(i, j, k, var), &
                                               1.0d0*var, &
                                               "Incorrect fine GC value")
                            end if
                          end if
                        end if
                      end do
                    end do
                  end do
                end do
            end associate
        end do
    end do

    call finish_test_run

end subroutine Driver_evolveFlash


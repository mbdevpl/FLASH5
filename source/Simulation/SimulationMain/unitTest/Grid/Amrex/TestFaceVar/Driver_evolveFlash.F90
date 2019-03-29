!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestFaceVar/Driver_evolveFlash
!!
!! NAME
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!  This test is meant to confirm proper guardcell filling of face variables as
!!  performed by the AMReX-based implementation of the Grid unit.  It confirms
!!  that all face variables have correct initial condition data, including
!!  face-centered data in the guardcells.  It then requests that all guardcells
!!  be filled for both X and Y face-centered data and confirms correct results.
!!
!!  Note that 
!!  - this simulation is run with only one refinement level so that we
!!  are not confirming full functionality of the the guardcell fill, but rather
!!  that the gr_fillPhysicalBC/single-level AMReX fillpatch routines are working.
!!
!! NOTES
!!  This simulation *must* be configured with at least the following
!!     ./setup -auto -2d -nxb=8 -nyb=4 
!!              unitTest/Grid/Amrex/TestFaceVar
!!             +noio -index-reorder
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Driver_evolveFlash()
    use amrex_fort_module,     ONLY : amrex_spacedim

    use Grid_interface,        ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                                      Grid_getLeafIterator, &
                                      Grid_releaseLeafIterator, &
                                      Grid_fillGuardCells, &
                                      Grid_getBlkBC
    use leaf_iterator,         ONLY : leaf_iterator_t
    use block_metadata,        ONLY : block_metadata_t
    use ut_testDriverMod

    implicit none

    integer, parameter :: N_CELLS_X = 16
    integer, parameter :: N_CELLS_Y = 8

    type(leaf_iterator_t)  :: itor
    type(block_metadata_t) :: block
    real, pointer          :: solnData(:, :, :, :) => null()
 
    integer :: i, j, var

    integer :: xBlkMin
    integer :: xBlkMax
    integer :: yBlkMin
    integer :: yBlkMax
    
    integer :: n_blocks
    integer :: boundaryType(LOW:HIGH, 1:MDIM)

    !!!!! CONFIRM PROPER COORDINATE SYSTEM
    ! Dimensionality
    write(*,*)
    if (amrex_spacedim /= NDIM) then
        write(*,*) "Wrong dimensionality - ", amrex_spacedim, ' != ', NDIM
        write(*,*) "Recompile AMReX with correct dimensionality"
        write(*,*)
        stop
    end if

    call start_test_run

    !!!!! CONFIRM PROPER GENERAL SETUP AS REQUIRED BY TEST
    call assertEqual(NDIM, 2, "Domain not 2D")
    call assertEqual(NFACE_VARS, 2, "Invalid number of facevars")

    !!!!! CONFIRM PROPER INITIAL CONDITIONS
    ! Facevar GCs not set during initialization
    n_blocks = 0
    call Grid_getLeafIterator(itor)
    call itor%blkMetaData(block)
    xBlkMin = block%limits(LOW,  IAXIS)
    xBlkMax = block%limits(HIGH, IAXIS)
    yBlkMin = block%limits(LOW,  JAXIS)
    yBlkMax = block%limits(HIGH, JAXIS)
    do while (itor%is_valid())
      call itor%blkMetaData(block)
      xBlkMin = MIN(xBlkMin, block%limits(LOW,  IAXIS))
      yBlkMin = MIN(yBlkMin, block%limits(LOW,  JAXIS))
      xBlkMax = MAX(xBlkMax, block%limits(HIGH, IAXIS))
      yBlkMax = MAX(yBlkMax, block%limits(HIGH, JAXIS))
      n_blocks = n_blocks + 1

      associate(lo => block%limitsGC(LOW, :), &
                hi => block%limitsGC(HIGH, :))
      call Grid_getBlkPtr(block, solnData, FACEX)
        do   j = lo(JAXIS), hi(JAXIS)
          do i = lo(IAXIS), hi(IAXIS)+1
            do var = 1, NFACE_VARS 
              if (      (lo(IAXIS)+NGUARD <= i) .AND. (i <= hi(IAXIS)-NGUARD+1) &
                  .AND. (lo(JAXIS)+NGUARD <= j) .AND. (j <= hi(JAXIS)-NGUARD)) then
                call assertEqual(solnData(i, j, 1, var), &
                                 1.3*i*i, &
                                 "Incorrect initial condition")
              else
                call assertEqual(solnData(i, j, 1, var), 0.0, &
                                 "Guardcell data not zero")
              end if
            end do
          end do
        end do
      call Grid_releaseBlkPtr(block, solnData, FACEX)

      call Grid_getBlkPtr(block, solnData, FACEY)
        do   j = lo(JAXIS), hi(JAXIS)+1
          do i = lo(IAXIS), hi(IAXIS)
            do var = 1, NFACE_VARS
              if (      (lo(IAXIS)+NGUARD <= i) .AND. (i <= hi(IAXIS)-NGUARD) &
                  .AND. (lo(JAXIS)+NGUARD <= j) .AND. (j <= hi(JAXIS)-NGUARD+1)) then
                call assertEqual(solnData(i, j, 1, var), &
                                 i*i + 1.2*j*j, &
                                 "Incorrect initial condition")
              else
                call assertEqual(solnData(i, j, 1, var), 0.0, &
                                 "Guardcell data not zero")
              end if
            end do
          end do
        end do
      call Grid_releaseBlkPtr(block,solnData,FACEY)
 
      end associate
      call itor%next()
    end do
    call Grid_releaseLeafIterator(itor)

    !!!!! CONFIRM PROPER BLOCK/CELL SETUP AS REQUIRED BY TEST
    call assertEqual(1, xBlkMin, "Incorrect starting x cell index")
    call assertEqual(N_CELLS_X, xBlkMax, "Invalid number of X cells")
    call assertEqual(1, yBlkMin, "Incorrect starting y cell index")
    call assertEqual(N_CELLS_Y, yBlkMax, "Invalid number of Y cells")

    ! Confirm that domain is decomposed on more than one block so that
    ! the following tests confirm correct tranfer of interior data to blocks
    ! with a face that is not on the boundary
    call assertTrue(n_blocks == 4, "Domain should have 4 blocks")

    !!!!! TRANSFER GUARDCELL DATA & FILL IN DATA OUTSIDE DOMAIN
    ! Set data on faces outside interior and not on boundaries
    ! using the specialized BC apply routine
    ! 
    ! NOTE: The data on these faces follow a different pattern from
    ! the interior/boundaries.  Therefore, we can confirm that we are
    ! not overwriting interior/boundary data.
    call Grid_fillGuardCells(FACEX, ALLDIR)
    call Grid_fillGuardCells(FACEY, ALLDIR)

    !!!!! CONFIRM PROPER FILLING OF FACE-CENTERED DATA
    call Grid_getLeafIterator(itor)
    do while (itor%is_valid())
      call itor%blkMetaData(block)

      call Grid_getBlkBC(block, boundaryType)

      associate(lo => block%limitsGC(LOW, :), &
                hi => block%limitsGC(HIGH, :))
        call Grid_getBlkPtr(block, solnData, FACEX)
        do   var = 1, NFACE_VARS 
          do   j = lo(JAXIS), hi(JAXIS)
            do i = lo(IAXIS), hi(IAXIS)+1
              if      (      (lo(IAXIS)+NGUARD <  i) .AND. (i <  hi(IAXIS)-NGUARD+1) &
                       .AND. (lo(JAXIS)+NGUARD <= j) .AND. (j <= hi(JAXIS)-NGUARD)) then
                call assertEqual(solnData(i, j, 1, var), 1.3*i*i, &
                                 "Incorrect data in interior")
              else if (      (boundaryType(LOW,  IAXIS) /= NOT_BOUNDARY) &
                       .AND. ( lo(IAXIS)+NGUARD == i) &
                       .AND. (1 <= j) .AND. (j <= N_CELLS_Y)) then
                call assertEqual(solnData(i, j, 1, var), 0.0, &
                                 "Incorrect data in low boundary")
              else if (      (boundaryType(HIGH, IAXIS) /= NOT_BOUNDARY) &
                       .AND. ( hi(IAXIS)-NGUARD+1 == i) &
                       .AND. (1 <= j) .AND. (j <= N_CELLS_Y)) then
                call assertEqual(solnData(i, j, 1, var), 0.0, &
                                 "Incorrect data in high boundary")
              else if (      (boundaryType(LOW,  IAXIS) /= NOT_BOUNDARY) &
                       .AND. (i < lo(IAXIS)+NGUARD)) then
                call assertEqual(solnData(i, j, 1, var), 1.1*var, &
                                 "Incorrect data for GC outside domain")
              else if (      (boundaryType(HIGH, IAXIS) /= NOT_BOUNDARY) &
                       .AND. (i > hi(IAXIS)-NGUARD+1)) then
                call assertEqual(solnData(i, j, 1, var), 1.1*var, &
                                 "Incorrect data for GC outside domain")
              else if (      (boundaryType(LOW,  JAXIS) /= NOT_BOUNDARY) &
                       .AND. (j < lo(JAXIS)+NGUARD)) then
                call assertEqual(solnData(i, j, 1, var), 1.1*var, &
                                 "Incorrect data for GC outside domain")
              else if (      (boundaryType(HIGH,  JAXIS) /= NOT_BOUNDARY) &
                       .AND. (j > hi(JAXIS)-NGUARD)) then
                call assertEqual(solnData(i, j, 1, var), 1.1*var, &
                                 "Incorrect data for GC outside domain")
              else
                call assertEqual(solnData(i, j, 1, var), 1.3*i*i, &
                                 "Incorrect data for GC that is in interior")
              end if
            end do
          end do
        end do
        call Grid_releaseBlkPtr(block, solnData, FACEX)

        call Grid_getBlkPtr(block, solnData, FACEY)
        do   var = 1, NFACE_VARS
          do   j = lo(JAXIS), hi(JAXIS)+1
            do i = lo(IAXIS), hi(IAXIS)
              if      (      (lo(IAXIS)+NGUARD <= i) .AND. (i <= hi(IAXIS)-NGUARD) &
                       .AND. (lo(JAXIS)+NGUARD <  j) .AND. (j <  hi(JAXIS)-NGUARD+1)) then
                call assertEqual(solnData(i, j, 1, var), i*i + 1.2*j*j, &
                                 "Incorrect initial condition")
              else if (      (boundaryType(LOW,  JAXIS) /= NOT_BOUNDARY) &
                       .AND. ( lo(JAXIS)+NGUARD == j) &
                       .AND. (1 <= i) .AND. (i <= N_CELLS_X)) then
                call assertEqual(solnData(i, j, 1, var), 0.0, &
                                 "Incorrect data in low boundary")
              else if (      (boundaryType(HIGH, JAXIS) /= NOT_BOUNDARY) &
                       .AND. ( hi(JAXIS)-NGUARD+1 == j) &
                       .AND. (1 <= i) .AND. (i <= N_CELLS_X)) then
                call assertEqual(solnData(i, j, 1, var), 0.0, &
                                 "Incorrect data in high boundary")
              else if (      (boundaryType(LOW,  IAXIS) /= NOT_BOUNDARY) &
                       .AND. (i < lo(IAXIS)+NGUARD)) then
                call assertEqual(solnData(i, j, 1, var), 1.1*var, &
                                 "Incorrect data for GC outside domain")
              else if (      (boundaryType(HIGH, IAXIS) /= NOT_BOUNDARY) &
                       .AND. (i > hi(IAXIS)-NGUARD)) then
                call assertEqual(solnData(i, j, 1, var), 1.1*var, &
                                 "Incorrect data for GC outside domain")
              else if (      (boundaryType(LOW,  JAXIS) /= NOT_BOUNDARY) &
                       .AND. (j < lo(JAXIS)+NGUARD)) then
                call assertEqual(solnData(i, j, 1, var), 1.1*var, &
                                 "Incorrect data for GC outside domain")
              else if (      (boundaryType(HIGH,  JAXIS) /= NOT_BOUNDARY) &
                       .AND. (j > hi(JAXIS)-NGUARD+1)) then
                call assertEqual(solnData(i, j, 1, var), 1.1*var, &
                                 "Incorrect data for GC outside domain")
              else
                call assertEqual(solnData(i, j, 1, var), i*i + 1.2*j*j, &
                                 "Incorrect initial condition in guardcell")
              end if
            end do
          end do
        end do
        call Grid_releaseBlkPtr(block, solnData, FACEY)

      end associate
      call itor%next()
    end do
    call Grid_releaseLeafIterator(itor)

    call finish_test_run

end subroutine Driver_evolveFlash


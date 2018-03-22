!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amr/TestCyl2/Driver_evolveFlash
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
!!  A subset of simulation configuration data is loaded into AMReX at
!!  initialization and is therefore owned by AMReX.  As a result, AMReX is used
!!  to provide these data values to client code through the associated public
!!  Grid_* and local gr_* interface accessor routines.
!!
!!  This code tests that AMReX is properly initialized for a 2D cylindrical domain by
!!  verifying correct results as obtained through the accessor routines.
!!
!! NOTES
!!  This simulation *must* be configured with at least the following
!!  2D run:
!!     ./setup -auto -2d -nxb=8 -nyb=4 
!!              unitTest/Grid/Amr/TestCyl2d
!!             +noio -index-reorder
!!
!!  For the future:
!!             -unit=IO/IOMain/hdf5/serial/AM
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Driver_evolveFlash()
    use Grid_interface,        ONLY : Grid_getDomainBoundBox, &
                                      Grid_getBlkBoundBox, &
                                      Grid_getSingleCellCoords, &
                                      Grid_getGeometry, &
                                      Grid_getDeltas, &
                                      Grid_getMaxRefinement, &
                                      Grid_getBlkPtr, Grid_releaseBlkPtr, &
                                      Grid_updateRefinement, &
                                      Grid_getLeafIterator, &
                                      Grid_releaseLeafIterator, &
                                      Grid_getSingleCellVol
    use Grid_data,             ONLY : gr_iguard, gr_jguard, gr_kguard, &
                                      gr_meshMe, &
                                      gr_numRefineVarsMax, gr_numRefineVars, &
                                      gr_refine_var, &
                                      gr_refine_cutoff, gr_derefine_cutoff, &
                                      gr_refine_filter, &
                                      gr_enforceMaxRefinement, &
                                      gr_eosMode, &
                                      gr_eosModeInit
    use leaf_iterator,         ONLY : leaf_iterator_t
    use block_metadata,        ONLY : block_metadata_t
    use ut_testDriverMod

    implicit none

    !!!!! EXPECTED RESULTS BASED ON flash.par AND SETUP VALUES GIVEN ABOVE
    integer,  parameter :: NXCELL_EX   = 64
    integer,  parameter :: NYCELL_EX   = 64
    integer,  parameter :: NZCELL_EX   =  4
    integer,  parameter :: NXBLK_EX    =  8
    integer,  parameter :: NYBLK_EX    = 16
    integer,  parameter :: NZBLK_EX    =  2
    real,     parameter :: XMIN_EX     =  0.00d0
    real,     parameter :: XMAX_EX     =  2.00d0
    real,     parameter :: YMIN_EX     = -1.50d0
    real,     parameter :: YMAX_EX     =  4.50d0
    real,     parameter :: XDELTA_EX   = (XMAX_EX-XMIN_EX)/NXCELL_EX
    real,     parameter :: YDELTA_EX   = (YMAX_EX-YMIN_EX)/NYCELL_EX
    integer,  parameter :: MAXLEVEL_EX =  1
    integer,  parameter :: XL_BC_EX    = REFLECTING
    integer,  parameter :: XH_BC_EX    = OUTFLOW
    integer,  parameter :: YL_BC_EX    = DIRICHLET
    integer,  parameter :: YH_BC_EX    = DIODE

    integer :: geometry = -100
    real    :: domain(LOW:HIGH, MDIM) = 0.0d0
    integer :: domainBC(LOW:HIGH, MDIM) = PERIODIC 
    real    :: deltas(1:MDIM) = 0.0d0
    real    :: x_expected = 0.0d0
    real    :: y_expected = 0.0d0
    integer :: max_level = -1

    type(leaf_iterator_t)  :: itor
    type(block_metadata_t) :: block
    real, pointer          :: solnData(:, :, :, :) => null()
    integer                :: n_blocks = 0
    integer                :: blkLimits(LOW:HIGH, 1:MDIM) = 0
    integer                :: blkLimitsGC(LOW:HIGH, 1:MDIM) = 0
    integer                :: blkGC(LOW:HIGH, 1:MDIM) = 0
    integer                :: blkSize(1:MDIM) = 0
    integer                :: xBlkMin = 0
    integer                :: xBlkMax = 0
    integer                :: yBlkMin = 0
    integer                :: yBlkMax = 0
    integer                :: zBlkMin = 0
    integer                :: zBlkMax = 0
    real                   :: xMin = 0.0d0
    real                   :: xMax = 0.0d0
    real                   :: yMin = 0.0d0
    real                   :: yMax = 0.0d0
    real                   :: zMin = 0.0d0
    real                   :: zMax = 0.0d0
    real                   :: boundBox(LOW:HIGH, 1:MDIM) = 0.0d0
    real                   :: c_lo(1:MDIM) = 0.0d0
    real                   :: c_hi(1:MDIM) = 0.0d0
    real                   :: c_gc(1:MDIM) = 0.0d0 
    real                   :: c_itr(1:MDIM) = 0.0d0
    real                   :: x_coords(4) = 0.0d0
    real                   :: y_coords(4) = 0.0d0
    real                   :: x_coords_gc(8) = 0.0d0
    real                   :: y_coords_gc(8) = 0.0d0
    real                   :: volume = 0.0d0
    real                   :: r = 0.0d0

    integer :: rank = -1
    integer :: ilev = 0
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    integer :: var = 0

    call start_test_run

    ! Physical domain
    call Grid_getGeometry(geometry)
    call assertEqual(geometry, CYLINDRICAL, "Incorrect coordinate system type")

    call Grid_getDomainBoundBox(domain)
    call assertEqual(domain(LOW,  IAXIS), XMIN_EX,"Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, IAXIS), XMAX_EX,"Incorrect high X-coordinate")
    call assertEqual(domain(LOW,  JAXIS), YMIN_EX,"Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, JAXIS), YMAX_EX,"Incorrect high Y-coordinate")
    call assertEqual(domain(LOW,  KAXIS), 0.0d0,  "Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, KAXIS), 0.0d0,  "Incorrect high Z-coordinate")

    !!!!! CONFIRM PROPER REFINEMENT
    call Grid_getMaxRefinement(max_level, mode=1)
    call assertEqual(max_level, MAXLEVEL_EX, "Incorrect maximum refine level")

    do ilev = 1, max_level
        call Grid_getDeltas(ilev, deltas)

        x_expected = XDELTA_EX / 2.0d0**(ilev - 1)
        y_expected = YDELTA_EX / 2.0d0**(ilev - 1)
        call assertEqual(deltas(IAXIS),x_expected,"Incorrect high X-coordinate")
        call assertEqual(deltas(JAXIS),y_expected,"Incorrect high Y-coordinate")
        call assertEqual(deltas(KAXIS),0.0d0,     "Incorrect high Z-coordinate")
    end do

    !!!!! CONFIRM PROPER BLOCK/CELL STRUCTURE
    ! Walk across all blocks to test and collect info
    n_blocks = 0
    
    call Grid_getLeafIterator(itor)

    call itor%blkMetaData(block)
    xBlkMin = block%limits(LOW,  IAXIS)
    xBlkMax = block%limits(HIGH, IAXIS)
    yBlkMin = block%limits(LOW,  JAXIS)
    yBlkMax = block%limits(HIGH, JAXIS)
    zBlkMin = block%limits(LOW,  KAXIS)
    zBlkMax = block%limits(HIGH, KAXIS)
    ! DEV: TODO Do better than this
    xMin = 1.0d10
    xMax = -xMin
    yMin = 1.0d10
    yMax = -yMin
    zMin = 1.0d10
    zMax = -zMin
    do while (itor%is_valid())
        n_blocks = n_blocks + 1
        call itor%blkMetaData(block)

        call Grid_getBlkBoundBox(block, boundBox)
        xMin = MIN(xMin, boundBox(LOW,  IAXIS))
        xMax = MAX(xMax, boundBox(HIGH, IAXIS))
        yMin = MIN(yMin, boundBox(LOW,  JAXIS))
        yMax = MAX(yMax, boundBox(HIGH, JAXIS))
        zMin = MIN(zMin, boundBox(LOW,  KAXIS))
        zMax = MAX(zMax, boundBox(HIGH, KAXIS))

        call assertEqual(block%level, 1, "Incorrect block level")

        ! Check guard cells along all directions
        blkLimits   = block%limits
        blkLimitsGC = block%limitsGC
        blkGC(LOW, :) = blkLimits(LOW, :) - blkLimitsGC(LOW, :)
        blkGC(HIGH, :) = blkLimitsGC(HIGH, :) - blkLimits(HIGH, :)
        call assertEqual(blkGC(LOW,  IAXIS), gr_iguard, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(HIGH, IAXIS), gr_iguard, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(LOW,  JAXIS), gr_jguard, &
                         "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(HIGH, JAXIS), gr_jguard, &
                         "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(LOW,  KAXIS), 0, "Incorrect guard cell along Z-axis")
        call assertEqual(blkGC(HIGH, KAXIS), 0, "Incorrect guard cell along Z-axis")

        ! Correct cells per block along each direction
        blkSize = blkLimits(HIGH, :) - blkLimits(LOW, :) + 1
        call assertEqual(blkSize(IAXIS), NXCELL_EX / NXBLK_EX, &
                         "Incorrect cells per block along X-axis")
        call assertEqual(blkSize(JAXIS), NYCELL_EX / NYBLK_EX, &
                         "Incorrect cells per block along Y-axis")
        call assertEqual(blkSize(KAXIS), 1, "Incorrect cells per block along Z-axis")

        xBlkMin = MIN(xBlkMin, blkLimits(LOW,  IAXIS))
        yBlkMin = MIN(yBlkMin, blkLimits(LOW,  JAXIS))
        zBlkMin = MIN(zBlkMin, blkLimits(LOW,  KAXIS))
        xBlkMax = MAX(xBlkMax, blkLimits(HIGH, IAXIS))
        yBlkMax = MAX(yBlkMax, blkLimits(HIGH, JAXIS))
        zBlkMax = MAX(zBlkMax, blkLimits(HIGH, KAXIS))

        call itor%next()
    end do

    call Grid_releaseLeafIterator(itor)
    
    ! Confirm proper number of blocks and cells
    call assertEqual(xBlkMin, 1, "Incorrect origin X-coordinate")
    call assertEqual(yBlkMin, 1, "Incorrect origin Y-coordinate")
    call assertEqual(zBlkMin, 1, "Incorrect origin Z-coordinate")

    ! FIXME: This only works for one processor.  Need a reduction here.
    call assertEqual(n_blocks, NXBLK_EX*NYBLK_EX, &
                     "Incorrect total number of blocks")

    call assertEqual(xBlkMax, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    call assertEqual(yBlkMax, NYCELL_EX, &
                     "Incorrect total number of cells along Y-axis")
    call assertEqual(zBlkMax, 1, "More than one cell along Z-axis")
    
    call assertEqual(xMin, XMIN_EX, "Incorrect minimum X-coordinate found")
    call assertEqual(xMax, XMAX_EX, "Incorrect maximum X-coordinate found")
    call assertEqual(yMin, YMIN_EX, "Incorrect minimum Y-coordinate found")
    call assertEqual(yMax, YMAX_EX, "Incorrect maximum Y-coordinate found")
    call assertEqual(zMin, 1.0d0,   "Incorrect minimum Z-coordinate found")
    call assertEqual(zMax, 1.0d0,   "Incorrect maximum Z-coordinate found")

    !!!!! CONFIRM PROPER BC
    call Grid_getDomainBC(domainBC)
    call assertEqual(domainBC(LOW,  IAXIS), XL_BC_EX, "Incorrect X-left BC")
    call assertEqual(domainBC(HIGH, IAXIS), XH_BC_EX, "Incorrect X-right BC")
    call assertEqual(domainBC(LOW,  JAXIS), YL_BC_EX, "Incorrect Y-left BC")
    call assertEqual(domainBC(HIGH, JAXIS), YH_BC_EX, "Incorrect Y-right BC")

    !!!!! CONFIRM PROPER INITIAL CONDITIONS
    call Grid_getLeafIterator(itor)
    do while (itor%is_valid())
        call itor%blkMetaData(block)
        call Grid_getBlkPtr(block, solnData)

        associate(lo => block%limits(LOW, :), &
                  hi => block%limits(HIGH, :))
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        do var = UNK_VARS_BEGIN, UNK_VARS_END 
                            call assertEqual(solnData(i, j, k, var), &
                                             1.1d0 * var, &
                                             "Incorrect initial condition")
                        end do
                    end do
                end do
            end do
        end associate

        call Grid_releaseBlkPtr(block, solnData)

        call itor%next()
    end do
    call Grid_releaseLeafIterator(itor)

    !!!!! CONFIRM CELL COORDINATE ACCESSORS
    ! Find coordinates of lo/hi
    block%level = 1
    block%limits(LOW,  :) = [1, 1, 1]
    call Grid_getSingleCellCoords([1, 1, 1], block, LEFT_EDGE, INTERIOR, c_lo)
    block%limits(LOW,  :) = [57, 61, 2]
    call Grid_getSingleCellCoords([8, 4, 1], block, RIGHT_EDGE, INTERIOR, c_hi)
    ! Find coordinates of cell as a guard cell of one block ...
    block%limits(LOW,  :) = [57, 61, 2]
    call Grid_getSingleCellCoords([1, 1, 1], block, CENTER, EXTERIOR, c_gc)
    ! and as an interior cell of its neighboring block
    block%limits(LOW,  :) = [49, 57, 2]
    call Grid_getSingleCellCoords([7, 3, 1], block, CENTER, INTERIOR, c_itr)

    call assertEqual(c_gc(IAXIS), c_itr(IAXIS), "Invalid cell X-coordinate")
    call assertEqual(c_gc(JAXIS), c_itr(JAXIS), "Invalid cell Y-coordinate")
    call assertEqual(c_gc(KAXIS), c_itr(KAXIS), "Invalid cell Z-coordinate")

    call assertEqual(c_lo(IAXIS), XMIN_EX, "Invalid cell X-coordinate")
    call assertEqual(c_lo(JAXIS), YMIN_EX, "Invalid cell Y-coordinate")
    call assertEqual(c_lo(KAXIS), 0.0d0,   "Invalid cell Z-coordinate")
    
    call assertEqual(c_hi(IAXIS), XMAX_EX, "Invalid cell X-coordinate")
    call assertEqual(c_hi(JAXIS), YMAX_EX, "Invalid cell Y-coordinate")
    call assertEqual(c_hi(KAXIS), 0.0d0,   "Invalid cell Z-coordinate")

    call assertEqual(c_gc(IAXIS), XMAX_EX - 9.5*XDELTA_EX, "Invalid cell X-coordinate")
    call assertEqual(c_gc(JAXIS), YMAX_EX - 5.5*YDELTA_EX, "Invalid cell Y-coordinate")
    call assertEqual(c_gc(KAXIS), 0.0d0,      "Invalid cell Z-coordinate")

    ! Check cell coordinates by block starting at lower corner
    block%limits(LOW, :) = [1, 1, 1]
    call Grid_getCellCoords(IAXIS, block, LEFT_EDGE, .FALSE., &
                            x_coords, SIZE(x_coords))
    call Grid_getCellCoords(JAXIS, block, LEFT_EDGE, .FALSE., &
                            y_coords, SIZE(y_coords))
    do j = 1, SIZE(x_coords)
        call assertEqual(x_coords(j), XMIN_EX + (j-1)*XDELTA_EX, "Bad X-coordinate")
        call assertEqual(y_coords(j), YMIN_EX + (j-1)*YDELTA_EX, "Bad Y-coordinate")
    end do

    ! Check cell coordinates by block starting at upper corner
    block%limits(LOW, :) = [61, 61, 1]
    call Grid_getCellCoords(IAXIS, block, RIGHT_EDGE, .FALSE., &
                            x_coords, SIZE(x_coords))
    call Grid_getCellCoords(JAXIS, block, RIGHT_EDGE, .FALSE., &
                            y_coords, SIZE(y_coords))
    do j = 1, SIZE(x_coords)
        call assertEqual(x_coords(j), XMAX_EX - (4-j)*XDELTA_EX, "Bad X-coordinate")
        call assertEqual(y_coords(j), YMAX_EX - (4-j)*YDELTA_EX, "Bad Y-coordinate")
    end do

    block%limits(LOW, :) = [61, 61, 1]
    call Grid_getCellCoords(IAXIS, block, RIGHT_EDGE, .TRUE., &
                            x_coords_gc, SIZE(x_coords_gc))
    call Grid_getCellCoords(JAXIS, block, RIGHT_EDGE, .TRUE., &
                            y_coords_gc, SIZE(y_coords_gc))
    do j = 1, SIZE(x_coords_gc)
        call assertEqual(x_coords_gc(j), XMAX_EX - (4-j+NGUARD)*XDELTA_EX, "Bad X-coordinate")
        call assertEqual(y_coords_gc(j), YMAX_EX - (4-j+NGUARD)*YDELTA_EX, "Bad Y-coordinate")
    end do
 
    !!!!! CELL VOLUMES
    associate(dr => XDELTA_EX, &
              dz => YDELTA_EX)
        ! Cell volume should depend on r
        block%level = 1
        block%limits(LOW,  :) = [  1,   1, 1]
        block%limits(HIGH, :) = [NXB, NYB, 1]
        call Grid_getSingleCellVol(block, [1, 1, 1], volume, INTERIOR)
        r = 0.5d0*dr
        call assertEqual(volume, 2.0d0*PI*r*dr*dz, "Invalid cell volume")
        call Grid_getSingleCellVol(block, [2, 1, 1], volume, INTERIOR)
        r = 1.5d0*dr
        call assertEqual(volume, 2.0d0*PI*r*dr*dz, "Invalid cell volume")
        call Grid_getSingleCellVol(block, [3, 1, 1], volume, INTERIOR)
        r = 2.5d0*dr
        call assertEqual(volume, 2.0d0*PI*r*dr*dz, "Invalid cell volume")
        call Grid_getSingleCellVol(block, [4, 1, 1], volume, INTERIOR)
        r = 3.5d0*dr
        call assertEqual(volume, 2.0d0*PI*r*dr*dz, "Invalid cell volume")

        ! Moving through cells along z should not change volume
        call Grid_getSingleCellVol(block, [1, 2, 1], volume, INTERIOR)
        r = 0.5d0*dr
        call assertEqual(volume, 2.0d0*PI*r*dr*dz, "Invalid cell volume")
        call Grid_getSingleCellVol(block, [4, 2, 1], volume, INTERIOR)
        r = 3.5d0*dr
        call assertEqual(volume, 2.0d0*PI*r*dr*dz, "Invalid cell volume")
    
        ! Look at guardcells at negative radius
        block%level = 1
        block%limits(LOW,  :) = [  1,   1, 1]
        block%limits(HIGH, :) = [NXB, NYB, 1]
        call Grid_getSingleCellVol(block, [2, 1, 1], volume, EXTERIOR)
        r = 0.5d0*dr
        call assertEqual(volume, 2.0d0*PI*r*dr*dz, "Invalid cell volume")
        call Grid_getSingleCellVol(block, [1, 1, 1], volume, EXTERIOR)
        r = 1.5d0*dr
        call assertEqual(volume, 2.0d0*PI*r*dr*dz, "Invalid cell volume")
       
        ! No change along z
        call Grid_getSingleCellVol(block, [2, 2, 1], volume, EXTERIOR)
        r = 0.5d0*dr
        call assertEqual(volume, 2.0d0*PI*r*dr*dz, "Invalid cell volume")
        call Grid_getSingleCellVol(block, [1, 2, 1], volume, EXTERIOR)
        r = 1.5d0*dr
        call assertEqual(volume, 2.0d0*PI*r*dr*dz, "Invalid cell volume")
    end associate

    !!!!! TODO: Check face areas?

    call finish_test_run

end subroutine Driver_evolveFlash


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
                                      Grid_getSingleCellCoords, &
                                      Grid_getCellCoords, &
                                      Grid_getCellVolumes, &
                                      Grid_getGeometry, &
                                      Grid_getDeltas, &
                                      Grid_getMaxRefinement, &
                                      Grid_getTileIterator, &
                                      Grid_releaseTileIterator, &
                                      Grid_getSingleCellVol
    use Grid_iterator,         ONLY : Grid_iterator_t
    use Grid_tile,             ONLY : Grid_tile_t
    use ut_testDriverMod

    implicit none

    !!!!! EXPECTED RESULTS BASED ON flash.par AND SETUP VALUES GIVEN ABOVE
    integer,  parameter :: NXCELL_EX   = 64
    integer,  parameter :: NYCELL_EX   = 64
    integer,  parameter :: NZCELL_EX   =  4
    integer,  parameter :: NXBLK_EX    =  8
    integer,  parameter :: NYBLK_EX    = 16
    integer,  parameter :: NZBLK_EX    =  2
    real,     parameter :: XMIN_EX     =  0.00
    real,     parameter :: XMAX_EX     =  2.00
    real,     parameter :: YMIN_EX     = -1.50
    real,     parameter :: YMAX_EX     =  4.50
    real,     parameter :: XDELTA_EX   = (XMAX_EX-XMIN_EX)/NXCELL_EX
    real,     parameter :: YDELTA_EX   = (YMAX_EX-YMIN_EX)/NYCELL_EX
    integer,  parameter :: MAXLEVEL_EX =  1
    integer,  parameter :: XL_BC_EX    = REFLECTING
    integer,  parameter :: XH_BC_EX    = OUTFLOW
    integer,  parameter :: YL_BC_EX    = DIRICHLET
    integer,  parameter :: YH_BC_EX    = DIODE

    integer :: geometry
    real    :: domain(LOW:HIGH, MDIM)
    integer :: domainBC(LOW:HIGH, MDIM)
    real    :: deltas(1:MDIM)
    real    :: x_expected
    real    :: y_expected
    integer :: max_level

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t)     :: tileDesc
    real, pointer         :: solnData(:, :, :, :)
    integer               :: n_blocks
    integer               :: blkLimits(LOW:HIGH, 1:MDIM)
    integer               :: blkLimitsGC(LOW:HIGH, 1:MDIM)
    integer               :: blkGC(LOW:HIGH, 1:MDIM)
    integer               :: blkSize(1:MDIM)
    integer               :: xBlkMin
    integer               :: xBlkMax
    integer               :: yBlkMin
    integer               :: yBlkMax
    integer               :: zBlkMin
    integer               :: zBlkMax
    real                  :: xMin
    real                  :: xMax
    real                  :: yMin
    real                  :: yMax
    real                  :: zMin
    real                  :: zMax
    real                  :: boundBox(LOW:HIGH, 1:MDIM)
    real                  :: c_lo(1:MDIM)
    real                  :: c_hi(1:MDIM)
    real                  :: volume
    real                  :: r
 
    real, allocatable :: x_coords(:)
    real, allocatable :: y_coords(:)
    real, allocatable :: z_coords(:)
    real, allocatable :: volumes(:, :, :)
    
    integer :: offset(1:MDIM)
    integer :: lo(1:MDIM)
    integer :: hi(1:MDIM)

    integer :: rank
    integer :: ilev
    integer :: i, j, k, var

    nullify(solnData)

    call start_test_run

    ! Physical domain
    call Grid_getGeometry(geometry)
    call assertEqual(geometry, CYLINDRICAL, "Incorrect coordinate system type")

    call Grid_getDomainBoundBox(domain)
    call assertEqual(domain(LOW,  IAXIS), XMIN_EX,"Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, IAXIS), XMAX_EX,"Incorrect high X-coordinate")
    call assertEqual(domain(LOW,  JAXIS), YMIN_EX,"Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, JAXIS), YMAX_EX,"Incorrect high Y-coordinate")
    call assertEqual(domain(LOW,  KAXIS), 0.0,    "Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, KAXIS), 0.0,    "Incorrect high Z-coordinate")

    !!!!! CONFIRM PROPER REFINEMENT
    call Grid_getMaxRefinement(max_level, mode=1)
    call assertEqual(max_level, MAXLEVEL_EX, "Incorrect maximum refine level")

    do ilev = 1, max_level
        call Grid_getDeltas(ilev, deltas)

        x_expected = XDELTA_EX / 2.0**(ilev - 1)
        y_expected = YDELTA_EX / 2.0**(ilev - 1)
        call assertEqual(deltas(IAXIS),x_expected,"Incorrect high X-coordinate")
        call assertEqual(deltas(JAXIS),y_expected,"Incorrect high Y-coordinate")
        call assertEqual(deltas(KAXIS),0.0,       "Incorrect high Z-coordinate")
    end do

    !!!!! CONFIRM PROPER BLOCK/CELL STRUCTURE
    ! Walk across all blocks to test and collect info
    n_blocks = 0

    ! The tests that use this iterator are testing for block info
    ! Do not use tiling here
    call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)

    call itor%currentTile(tileDesc)
    xBlkMin = tileDesc%limits(LOW,  IAXIS)
    xBlkMax = tileDesc%limits(HIGH, IAXIS)
    yBlkMin = tileDesc%limits(LOW,  JAXIS)
    yBlkMax = tileDesc%limits(HIGH, JAXIS)
    zBlkMin = tileDesc%limits(LOW,  KAXIS)
    zBlkMax = tileDesc%limits(HIGH, KAXIS)
    ! DEV: TODO Do better than this
    xMin = 1.0e10
    xMax = -xMin
    yMin = 1.0e10
    yMax = -yMin
    zMin = 1.0e10
    zMax = -zMin
    do while (itor%isValid())
        n_blocks = n_blocks + 1
        call itor%currentTile(tileDesc)

        call tileDesc%boundBox(boundBox)
        xMin = MIN(xMin, boundBox(LOW,  IAXIS))
        xMax = MAX(xMax, boundBox(HIGH, IAXIS))
        yMin = MIN(yMin, boundBox(LOW,  JAXIS))
        yMax = MAX(yMax, boundBox(HIGH, JAXIS))
        zMin = MIN(zMin, boundBox(LOW,  KAXIS))
        zMax = MAX(zMax, boundBox(HIGH, KAXIS))

        call assertEqual(tileDesc%level, 1, "Incorrect block level")

        ! Check guard cells along all directions
        blkLimits   = tileDesc%limits
        blkLimitsGC = tileDesc%blkLimitsGC
        blkGC(LOW, :) = blkLimits(LOW, :) - blkLimitsGC(LOW, :)
        blkGC(HIGH, :) = blkLimitsGC(HIGH, :) - blkLimits(HIGH, :)
        call assertEqual(blkGC(LOW,  IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(HIGH, IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(LOW,  JAXIS), NGUARD, &
                         "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(HIGH, JAXIS), NGUARD, &
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
    call assertEqual(zMin, 1.0,     "Incorrect minimum Z-coordinate found")
    call assertEqual(zMax, 1.0,     "Incorrect maximum Z-coordinate found")

    !!!!! CONFIRM PROPER BC
    call Grid_getDomainBC(domainBC)
    call assertEqual(domainBC(LOW,  IAXIS), XL_BC_EX, "Incorrect X-left BC")
    call assertEqual(domainBC(HIGH, IAXIS), XH_BC_EX, "Incorrect X-right BC")
    call assertEqual(domainBC(LOW,  JAXIS), YL_BC_EX, "Incorrect Y-left BC")
    call assertEqual(domainBC(HIGH, JAXIS), YH_BC_EX, "Incorrect Y-right BC")

    !!!!! CONFIRM PROPER INITIAL CONDITIONS
    call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
    do while (itor%isValid())
        call itor%currentTile(tileDesc)
        call tileDesc%getDataPtr(solnData, CENTER)

        associate(lo => tileDesc%limits(LOW, :), &
                  hi => tileDesc%limits(HIGH, :))
            do           var = UNK_VARS_BEGIN, UNK_VARS_END 
                do         k = lo(KAXIS), hi(KAXIS)
                    do     j = lo(JAXIS), hi(JAXIS)
                        do i = lo(IAXIS), hi(IAXIS)
                            call assertEqual(solnData(i, j, k, var), &
                                             1.1 * var, &
                                             "Incorrect initial condition")
                        end do
                    end do
                end do
            end do
        end associate

        call tileDesc%releaseDataPtr(solnData, CENTER)

        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

    !!!!! CONFIRM CELL COORDINATE ACCESSORS
    ! Find coordinates of lo/hi
    call Grid_getSingleCellCoords([1, 1, 1], 1, LEFT_EDGE, c_lo)
    call Grid_getSingleCellCoords([NXCELL_EX, NYCELL_EX, NZCELL_EX], 1, RIGHT_EDGE, c_hi)

    call assertEqual(c_lo(IAXIS), XMIN_EX, "Invalid cell X-coordinate")
    call assertEqual(c_lo(JAXIS), YMIN_EX, "Invalid cell Y-coordinate")
    call assertEqual(c_lo(KAXIS), 0.0,   "Invalid cell Z-coordinate")
    
    call assertEqual(c_hi(IAXIS), XMAX_EX, "Invalid cell X-coordinate")
    call assertEqual(c_hi(JAXIS), YMAX_EX, "Invalid cell Y-coordinate")
    call assertEqual(c_hi(KAXIS), 0.0,   "Invalid cell Z-coordinate")
    
    call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
    do while (itor%isValid())
       call itor%currentTile(tileDesc)

       lo = tileDesc%limits(LOW,  :)
       hi = tileDesc%limits(HIGH, :)
       allocate(x_coords(lo(IAXIS):hi(IAXIS)))
       allocate(y_coords(lo(JAXIS):hi(JAXIS)))
       allocate(z_coords(lo(KAXIS):hi(KAXIS)))

       call Grid_getCellCoords(IAXIS, LEFT_EDGE, tileDesc%level, &
                               lo, hi, x_coords)
       call Grid_getCellCoords(JAXIS, LEFT_EDGE, tileDesc%level, &
                               lo, hi, y_coords)
       call Grid_getCellCoords(KAXIS, LEFT_EDGE, tileDesc%level, &
                               lo, hi, z_coords)
       do       k = lo(KAXIS), hi(KAXIS)
          do    j = lo(JAXIS), hi(JAXIS)
             do i = lo(IAXIS), hi(IAXIS)
                  call assertEqual(x_coords(i), XMIN_EX + (i-1)*XDELTA_EX, "Bad X-coordinate")
                  call assertEqual(y_coords(j), YMIN_EX + (j-1)*YDELTA_EX, "Bad Y-coordinate")
                  call assertEqual(z_coords(k), 0.0,                       "Bad Z-coordinate")
                  
                  call Grid_getSingleCellCoords([i, j, k], tileDesc%level, LEFT_EDGE, c_lo)
                  call assertEqual(x_coords(i), c_lo(IAXIS), "X-coordinate doesn't match")
                  call assertEqual(y_coords(j), c_lo(JAXIS), "Y-coordinate doesn't match")
                  call assertEqual(z_coords(k), c_lo(KAXIS), "Z-coordinate doesn't match")
             end do
          end do
       end do

       call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                               lo, hi, x_coords)
       call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                               lo, hi, y_coords)
       call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                               lo, hi, z_coords)
       do       k = lo(KAXIS), hi(KAXIS)
          do    j = lo(JAXIS), hi(JAXIS)
             do i = lo(IAXIS), hi(IAXIS)
                  call assertEqual(x_coords(i), XMIN_EX + (i-0.5)*XDELTA_EX, "Bad X-coordinate")
                  call assertEqual(y_coords(j), YMIN_EX + (j-0.5)*YDELTA_EX, "Bad Y-coordinate")
                  call assertEqual(z_coords(k), 0.0,                         "Bad Z-coordinate")
                  
                  call Grid_getSingleCellCoords([i, j, k], tileDesc%level, CENTER, c_lo)
                  call assertEqual(x_coords(i), c_lo(IAXIS), "X-coordinate doesn't match")
                  call assertEqual(y_coords(j), c_lo(JAXIS), "Y-coordinate doesn't match")
                  call assertEqual(z_coords(k), c_lo(KAXIS), "Z-coordinate doesn't match")
             end do
          end do
       end do

       call Grid_getCellCoords(IAXIS, RIGHT_EDGE, tileDesc%level, &
                               lo, hi, x_coords)
       call Grid_getCellCoords(JAXIS, RIGHT_EDGE, tileDesc%level, &
                               lo, hi, y_coords)
       call Grid_getCellCoords(KAXIS, RIGHT_EDGE, tileDesc%level, &
                               lo, hi, z_coords)
       do       k = lo(KAXIS), hi(KAXIS)
          do    j = lo(JAXIS), hi(JAXIS)
             do i = lo(IAXIS), hi(IAXIS)
                  call assertEqual(x_coords(i), XMIN_EX + i*XDELTA_EX, "Bad X-coordinate")
                  call assertEqual(y_coords(j), YMIN_EX + j*YDELTA_EX, "Bad Y-coordinate")
                  call assertEqual(z_coords(k), 0.0,                   "Bad Z-coordinate")
                  
                  call Grid_getSingleCellCoords([i, j, k], tileDesc%level, RIGHT_EDGE, c_lo)
                  call assertEqual(x_coords(i), c_lo(IAXIS), "X-coordinate doesn't match")
                  call assertEqual(y_coords(j), c_lo(JAXIS), "Y-coordinate doesn't match")
                  call assertEqual(z_coords(k), c_lo(KAXIS), "Z-coordinate doesn't match")
             end do
          end do
       end do

       deallocate(x_coords)
       deallocate(y_coords)
       deallocate(z_coords)

       call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

    !!!!! CELL VOLUMES
    associate(dr => XDELTA_EX, &
              dz => YDELTA_EX)
        ! Cell volume should depend on r
        call Grid_getSingleCellVol([1, 1, 1], 1, volume)
        r = 0.5*dr
        call assertEqual(volume, 2.0*PI*r*dr*dz, "Invalid cell volume")
        call Grid_getSingleCellVol([2, 1, 1], 1, volume)
        r = 1.5*dr
        call assertEqual(volume, 2.0*PI*r*dr*dz, "Invalid cell volume")
        call Grid_getSingleCellVol([3, 1, 1], 1, volume)
        r = 2.5*dr
        call assertEqual(volume, 2.0*PI*r*dr*dz, "Invalid cell volume")
        call Grid_getSingleCellVol([4, 1, 1], 1, volume)
        r = 3.5*dr
        call assertEqual(volume, 2.0*PI*r*dr*dz, "Invalid cell volume")

        ! Moving through cells along z should not change volume
        call Grid_getSingleCellVol([1, 2, 1], 1, volume)
        r = 0.5*dr
        call assertEqual(volume, 2.0*PI*r*dr*dz, "Invalid cell volume")
        call Grid_getSingleCellVol([4, 2, 1], 1, volume)
        r = 3.5*dr
        call assertEqual(volume, 2.0*PI*r*dr*dz, "Invalid cell volume")
    end associate

    call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
    do while (itor%isValid())
       call itor%currentTile(tileDesc)

       lo = tileDesc%limits(LOW,  :)
       hi = tileDesc%limits(HIGH, :)
       allocate(volumes(lo(IAXIS):hi(IAXIS), &
                        lo(JAXIS):hi(JAXIS), &
                        lo(KAXIS):hi(KAXIS)))

       call Grid_getCellVolumes(tileDesc%level, &
                                lbound(volumes), ubound(volumes), &
                                volumes)
       do       k = lo(KAXIS), hi(KAXIS)
          do    j = lo(JAXIS), hi(JAXIS)
             do i = lo(IAXIS), hi(IAXIS)
                  call Grid_getSingleCellVol([i, j, k], tileDesc%level, volume)
                  call assertEqual(volumes(i, j, k), volume, "Bad volume")
             end do
          end do
       end do

       deallocate(volumes)

       call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

    !!!!! TODO: Check face areas?

    call finish_test_run

end subroutine Driver_evolveFlash


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
!!  A subset of simulation configuration data is loaded into AMReX at
!!  initialization and is therefore owned by AMReX.  As a result, AMReX is used
!!  to provide these data values to client code through the associated public
!!  Grid_* and local gr_* interface accessor routines.
!!
!!  This code tests that AMReX is properly initialized for a Cartesian domain by
!!  verifying correct results as obtained through the accessor routines.
!!
!! NOTES
!!  This simulation *must* be configured with at least the following
!!  2D run:
!!     ./setup -auto -2d -nxb=8 -nyb=4 
!!              unitTest/Grid/Amrex/TestInit 
!!             +noio -index-reorder
!!  3D run:
!!     ./setup -auto -3d -nxb=8 -nyb=4 -nzb=2
!!              unitTest/Grid/Amrex/TestInit 
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
    use amrex_parallel_module, ONLY : amrex_parallel_myproc, &
                                      amrex_parallel_nprocs

    use Grid_interface,        ONLY : Grid_getDomainBoundBox, &
                                      Grid_getBlkBoundBox, &
                                      Grid_getSingleCellCoords, &
                                      Grid_getGeometry, &
                                      Grid_getDeltas, &
                                      Grid_getMaxRefinement, &
                                      Grid_getBlkPtr, Grid_releaseBlkPtr, &
                                      Grid_updateRefinement, &
                                      Grid_getLeafIterator, &
                                      Grid_releaseLeafIterator
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
    use block_metadata,        ONLY : block_metadata_t, bmd_print
    use ut_testDriverMod

    implicit none

    !!!!! EXPECTED RESULTS BASED ON flash.par AND SETUP VALUES GIVEN ABOVE
    integer,  parameter :: NXCELL_EX   = 64
    integer,  parameter :: NYCELL_EX   = 64
    integer,  parameter :: NZCELL_EX   =  4
    integer,  parameter :: NXBLK_EX    =  8
    integer,  parameter :: NYBLK_EX    = 16
    integer,  parameter :: NZBLK_EX    =  2
    real,     parameter :: XMIN_EX     = -1.00d0
    real,     parameter :: XMAX_EX     =  2.00d0
    real,     parameter :: YMIN_EX     = -1.50d0
    real,     parameter :: YMAX_EX     =  4.50d0
    real,     parameter :: ZMIN_EX     =  0.50d0
    real,     parameter :: ZMAX_EX     =  0.75d0
    real,     parameter :: XDELTA_EX   = (XMAX_EX-XMIN_EX)/NXCELL_EX
    real,     parameter :: YDELTA_EX   = (YMAX_EX-YMIN_EX)/NYCELL_EX
    real,     parameter :: ZDELTA_EX   = (ZMAX_EX-ZMIN_EX)/NZCELL_EX
    integer,  parameter :: MAXLEVEL_EX =  4
    integer,  parameter :: XL_BC_EX    = OUTFLOW
    integer,  parameter :: XH_BC_EX    = REFLECTING
    integer,  parameter :: YL_BC_EX    = DIRICHLET
    integer,  parameter :: YH_BC_EX    = DIODE
    integer,  parameter :: ZL_BC_EX    = PERIODIC
    integer,  parameter :: ZH_BC_EX    = DIRICHLET

    integer :: geometry = -100
    real    :: domain(LOW:HIGH, MDIM) = 0.0d0
    integer :: domainBC(LOW:HIGH, MDIM) = PERIODIC 
    real    :: deltas(1:MDIM) = 0.0d0
    real    :: x_expected = 0.0d0
    real    :: y_expected = 0.0d0
    real    :: z_expected = 0.0d0
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

    integer :: rank = -1
    integer :: ilev = 0
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    integer :: var = 0

    rank = amrex_parallel_myproc()

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

    !!!!! CONFIRM MPI SETUP
    call assertEqual(rank, gr_meshMe, "AMReX/FLASH ranks are different")
 
    ! Physical domain
    call Grid_getGeometry(geometry)
    call assertEqual(geometry, CARTESIAN, "Incorrect coordinate system type")

    call Grid_getDomainBoundBox(domain)
#if NDIM == 1
    call assertEqual(domain(LOW,  IAXIS), XMIN_EX,"Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, IAXIS), XMAX_EX,"Incorrect high X-coordinate")
    call assertEqual(domain(LOW,  JAXIS), 0.0d0,  "Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, JAXIS), 0.0d0,  "Incorrect high Y-coordinate")
    call assertEqual(domain(LOW,  KAXIS), 0.0d0,  "Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, KAXIS), 0.0d0,  "Incorrect high Z-coordinate")
#elif NDIM == 2
    call assertEqual(domain(LOW,  IAXIS), XMIN_EX,"Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, IAXIS), XMAX_EX,"Incorrect high X-coordinate")
    call assertEqual(domain(LOW,  JAXIS), YMIN_EX,"Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, JAXIS), YMAX_EX,"Incorrect high Y-coordinate")
    call assertEqual(domain(LOW,  KAXIS), 0.0d0,  "Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, KAXIS), 0.0d0,  "Incorrect high Z-coordinate")
#elif NDIM == 3 
    call assertEqual(domain(LOW,  IAXIS), XMIN_EX,"Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, IAXIS), XMAX_EX,"Incorrect high X-coordinate")
    call assertEqual(domain(LOW,  JAXIS), YMIN_EX,"Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, JAXIS), YMAX_EX,"Incorrect high Y-coordinate")
    call assertEqual(domain(LOW,  KAXIS), ZMIN_EX,"Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, KAXIS), ZMAX_EX,"Incorrect high Z-coordinate")
#endif

    !!!!! CONFIRM PROPER REFINEMENT
    call Grid_getMaxRefinement(max_level, mode=1)
    call assertEqual(max_level, MAXLEVEL_EX, "Incorrect maximum refine level")

    do ilev = 1, max_level
        call Grid_getDeltas(ilev, deltas)

        x_expected = XDELTA_EX / 2.0d0**(ilev - 1)
        y_expected = YDELTA_EX / 2.0d0**(ilev - 1)
        z_expected = ZDELTA_EX / 2.0d0**(ilev - 1)
#if NDIM == 1
        call assertEqual(deltas(IAXIS),x_expected,"Incorrect high X-coordinate")
        call assertEqual(deltas(JAXIS),0.0d0,     "Incorrect high Y-coordinate")
        call assertEqual(deltas(KAXIS),0.0d0,     "Incorrect high Z-coordinate")
#elif NDIM == 2
        call assertEqual(deltas(IAXIS),x_expected,"Incorrect high X-coordinate")
        call assertEqual(deltas(JAXIS),y_expected,"Incorrect high Y-coordinate")
        call assertEqual(deltas(KAXIS),0.0d0,     "Incorrect high Z-coordinate")
#elif NDIM == 3 
        call assertEqual(deltas(IAXIS),x_expected,"Incorrect high X-coordinate")
        call assertEqual(deltas(JAXIS),y_expected,"Incorrect high Y-coordinate")
        call assertEqual(deltas(KAXIS),z_expected,"Incorrect high Z-coordinate")
#endif
    end do

    ! DEV: TODO Once unittest is refining mesh, check Grid_getMaxRefinement
    ! with mode that checks actual number of levels in use
    
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
!        write(*,*) '-----------------------------------------'
!        write(*,*) 'Level = ', block%level
!        write(*,*) 'GID   = ', block%grid_index
!        write(*,*) 'Lower:  ', block%limits(LOW, :)
!        write(*,*) '        ', boundBox(LOW, :)
!        write(*,*) 'Upper:  ', block%limits(HIGH, :)
!        write(*,*) '        ', boundBox(HIGH, :)
        xMin = MIN(xMin, boundBox(LOW,  IAXIS))
        xMax = MAX(xMax, boundBox(HIGH, IAXIS))
        yMin = MIN(yMin, boundBox(LOW,  JAXIS))
        yMax = MAX(yMax, boundBox(HIGH, JAXIS))
        zMin = MIN(zMin, boundBox(LOW,  KAXIS))
        zMax = MAX(zMax, boundBox(HIGH, KAXIS))

        ! DEVNOTE: Should we leave this unittest with simple data
        ! that does not refine so that testing the block structure is easy?
        ! All blocks on coarsest level since no refining
!        call assertEqual(block%level, 1, "Incorrect block level")

        ! Check guard cells along all directions
        blkLimits   = block%limits
        blkLimitsGC = block%limitsGC
        blkGC(LOW, :) = blkLimits(LOW, :) - blkLimitsGC(LOW, :)
        blkGC(HIGH, :) = blkLimitsGC(HIGH, :) - blkLimits(HIGH, :)
#if NDIM == 1
        call assertEqual(blkGC(LOW,  IAXIS), gr_iguard, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(HIGH, IAXIS), gr_iguard, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(LOW,  JAXIS), 0, "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(HIGH, JAXIS), 0, "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(LOW,  KAXIS), 0, "Incorrect guard cell along Z-axis")
        call assertEqual(blkGC(HIGH, KAXIS), 0, "Incorrect guard cell along Z-axis")
#elif NDIM == 2
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
#elif NDIM == 3
        call assertEqual(blkGC(LOW,  IAXIS), gr_iguard, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(HIGH, IAXIS), gr_iguard, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(LOW,  JAXIS), gr_jguard, &
                         "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(HIGH, JAXIS), gr_jguard, &
                         "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(LOW,  KAXIS), gr_kguard, &
                         "Incorrect guard cell along Z-axis")
        call assertEqual(blkGC(HIGH, KAXIS), gr_kguard, &
                         "Incorrect guard cell along Z-axis")
#endif

        ! Correct cells per block along each direction
        blkSize = blkLimits(HIGH, :) - blkLimits(LOW, :) + 1
#if NDIM == 1
        call assertEqual(blkSize(IAXIS), NXCELL_EX / NXBLK_EX, &
                         "Incorrect cells per block along X-axis")
        call assertEqual(blkSize(JAXIS), 1, "Incorrect cells per block along Y-axis")
        call assertEqual(blkSize(KAXIS), 1, "Incorrect cells per block along Z-axis")
#elif NDIM == 2
        call assertEqual(blkSize(IAXIS), NXCELL_EX / NXBLK_EX, &
                         "Incorrect cells per block along X-axis")
        call assertEqual(blkSize(JAXIS), NYCELL_EX / NYBLK_EX, &
                         "Incorrect cells per block along Y-axis")
        call assertEqual(blkSize(KAXIS), 1, "Incorrect cells per block along Z-axis")
#elif NDIM == 3
        call assertEqual(blkSize(IAXIS), NXCELL_EX / NXBLK_EX, &
                         "Incorrect cells per block along X-axis")
        call assertEqual(blkSize(JAXIS), NYCELL_EX / NYBLK_EX, &
                         "Incorrect cells per block along Y-axis")
        call assertEqual(blkSize(KAXIS), NZCELL_EX / NZBLK_EX, &
                         "Incorrect cells per block along Z-axis")
#endif

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
#if NDIM == 1
    call assertEqual(n_blocks, NXBLK_EX, &
                     "Incorrect total number of blocks")

    call assertEqual(xBlkMax, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    call assertEqual(yBlkMax, 1, "More than one cell along Y-axis")
    call assertEqual(zBlkMax, 1, "More than one cell along Z-axis")
    
    call assertEqual(xMin, XMIN_EX, "Incorrect minimum X-coordinate found")
    call assertEqual(xMax, XMAX_EX, "Incorrect maximum X-coordinate found")
    call assertEqual(yMin, 1.0d0,   "Incorrect minimum Y-coordinate found")
    call assertEqual(yMax, 1.0d0,   "Incorrect maximum Y-coordinate found")
    call assertEqual(zMin, 1.0d0,   "Incorrect minimum Z-coordinate found")
    call assertEqual(zMax, 1.0d0,   "Incorrect maximum Z-coordinate found")
#elif NDIM == 2
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
#elif NDIM == 3
    call assertEqual(n_blocks, NXBLK_EX*NYBLK_EX*NZBLK_EX, &
                     "Incorrect total number of blocks")
    
    call assertEqual(xBlkMax, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    call assertEqual(yBlkMax, NYCELL_EX, &
                     "Incorrect total number of cells along Y-axis")
    call assertEqual(zBlkMax, NZCELL_EX, &
                     "Incorrect total number of cells along Z-axis")
    
    call assertEqual(xMin, XMIN_EX, "Incorrect minimum X-coordinate found")
    call assertEqual(xMax, XMAX_EX, "Incorrect maximum X-coordinate found")
    call assertEqual(yMin, YMIN_EX, "Incorrect minimum Y-coordinate found")
    call assertEqual(yMax, YMAX_EX, "Incorrect maximum Y-coordinate found")
    call assertEqual(zMin, ZMIN_EX, "Incorrect minimum Z-coordinate found")
    call assertEqual(zMax, ZMAX_EX, "Incorrect maximum Z-coordinate found")
#endif

    !!!!! CONFIRM PROPER BC
    call Grid_getDomainBC(domainBC)
    call assertEqual(domainBC(LOW,  IAXIS), XL_BC_EX, "Incorrect X-left BC")
    call assertEqual(domainBC(HIGH, IAXIS), XH_BC_EX, "Incorrect X-right BC")
    call assertEqual(domainBC(LOW,  JAXIS), YL_BC_EX, "Incorrect Y-left BC")
    call assertEqual(domainBC(HIGH, JAXIS), YH_BC_EX, "Incorrect Y-right BC")
    call assertEqual(domainBC(LOW,  KAXIS), ZL_BC_EX, "Incorrect Z-left BC")
    call assertEqual(domainBC(HIGH, KAXIS), ZH_BC_EX, "Incorrect Z-right BC")

    !!!!! CONFIRM REFINEMENT SETUP
    ! TODO: Get nrefs from AMReX

    !!!!! CONFIRM PROPER INITIAL CONDITIONS
    call Grid_getLeafIterator(itor)
    do while (itor%is_valid())
        call itor%blkMetaData(block)

        associate(lo => block%limits(LOW, :), &
                  hi => block%limits(HIGH, :))
        call Grid_getBlkPtr(block, solnData)
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        do var = UNK_VARS_BEGIN, UNK_VARS_END 
                            call assertEqual(solnData(i, j, k, var), &
                                             1.1d0 * var, &
                                             "Incorrect initial condition in unk")
                        end do
                    end do
                end do
            end do
        call Grid_releaseBlkPtr(block, solnData)
    !!!!! CONFIRM PROPER INITIAL CONDITIONS for FACEVARS
#if NFACE_VARS>0
        call Grid_getBlkPtr(block,solnData,FACEX)
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)+1
                        do var = 1, NFACE_VARS 
                            call assertEqual(solnData(i, j, k, var), &
                                             1.3*i*i, &
                                             "Incorrect initial condition")
                        end do
                    end do
                end do
            end do
        call Grid_releaseBlkPtr(block,solnData,FACEX)
#if NDIM>1
        call Grid_getBlkPtr(block,solnData,FACEY)
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)+1
                    do i = lo(IAXIS), hi(IAXIS)
                        do var = 1, NFACE_VARS
                            call assertEqual(solnData(i, j, k, var), &
                                             i*i+1.2**j*j, &
                                             "Incorrect initial condition")
                        end do
                    end do
                end do
            end do
        call Grid_releaseBlkPtr(block,solnData,FACEY)
#endif
#if NDIM>2
        call Grid_getBlkPtr(block,solnData,FACEZ)
            do         k = lo(KAXIS), hi(KAXIS)+1
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        do var = 1, NFACE_VARS
                            call assertEqual(solnData(i, j, k, var), &
                                             i*i+j*j+1.1*k*k, &
                                             "Incorrect initial condition")
                        end do
                    end do
                end do
            end do
        call Grid_releaseBlkPtr(block,solnData,FACEZ)
#endif
#endif
        end associate
        call itor%next()
    end do
    call Grid_releaseLeafIterator(itor)

    !!!!! CONFIRM CELL COORDINATE ACCESSORS
    ! DEV: TODO Apply proper Z-info and test in 3D
    ! Find coordinates of lo/hi
    block%level = 1
    block%limits(LOW,  :) = [1, 1, 1]
    call Grid_getSingleCellCoords([1, 1, 1], block, LEFT_EDGE, INTERIOR, c_lo)
    block%limits(LOW,  :) = [61, 61, 2]
    call Grid_getSingleCellCoords([4, 4, 1], block, RIGHT_EDGE, INTERIOR, c_hi)
    ! Find coordinates of cell as a guard cell of one block ...
    block%limits(LOW,  :) = [61, 61, 2]
    call Grid_getSingleCellCoords([1, 1, 1], block, CENTER, EXTERIOR, c_gc)
    ! and as an interior cell of its neighboring block
    block%limits(LOW,  :) = [57, 57, 2]
    call Grid_getSingleCellCoords([3, 3, 1], block, CENTER, INTERIOR, c_itr)

    call assertEqual(c_gc(IAXIS), c_itr(IAXIS), "Invalid cell X-coordinate")
    call assertEqual(c_gc(JAXIS), c_itr(JAXIS), "Invalid cell Y-coordinate")
    call assertEqual(c_gc(KAXIS), c_itr(KAXIS), "Invalid cell Z-coordinate")

#if NDIM == 1
    call assertEqual(c_lo(IAXIS), XMIN_EX, "Invalid cell X-coordinate")
    call assertEqual(c_lo(JAXIS), 0.0d0,   "Invalid cell Y-coordinate")
    call assertEqual(c_lo(KAXIS), 0.0d0,   "Invalid cell Z-coordinate")
    
    call assertEqual(c_hi(IAXIS), XMAX_EX, "Invalid cell X-coordinate")
    call assertEqual(c_hi(JAXIS), 0.0d0,   "Invalid cell Y-coordinate")
    call assertEqual(c_hi(KAXIS), 0.0d0,   "Invalid cell Z-coordinate")
#elif NDIM == 2
    call assertEqual(c_lo(IAXIS), XMIN_EX, "Invalid cell X-coordinate")
    call assertEqual(c_lo(JAXIS), YMIN_EX, "Invalid cell Y-coordinate")
    call assertEqual(c_lo(KAXIS), 0.0d0,   "Invalid cell Z-coordinate")
    
    call assertEqual(c_hi(IAXIS), XMAX_EX, "Invalid cell X-coordinate")
    call assertEqual(c_hi(JAXIS), YMAX_EX, "Invalid cell Y-coordinate")
    call assertEqual(c_hi(KAXIS), 0.0d0,   "Invalid cell Z-coordinate")

    call assertEqual(c_gc(IAXIS), 1.7421875, "Invalid cell X-coordinate")
    call assertEqual(c_gc(JAXIS), 3.9843750, "Invalid cell Y-coordinate")
    call assertEqual(c_gc(KAXIS), 0.0d0,      "Invalid cell Z-coordinate")
#elif NDIM == 3
    call assertEqual(c_lo(IAXIS), XMIN_EX, "Invalid cell X-coordinate")
    call assertEqual(c_lo(JAXIS), YMIN_EX, "Invalid cell Y-coordinate")
    call assertEqual(c_lo(KAXIS), ZMIN_EX, "Invalid cell Z-coordinate")
    
    call assertEqual(c_hi(IAXIS), XMAX_EX, "Invalid cell X-coordinate")
    call assertEqual(c_hi(JAXIS), YMAX_EX, "Invalid cell Y-coordinate")
    call assertEqual(c_hi(KAXIS), ZMAX_EX, "Invalid cell Z-coordinate")
#endif

    ! DEV: TODO Implement for different NDIM and with proper Z-info
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

    !!!!! CONFIRM REFINEMENT SETUP
    ! uses default value
    call assertEqual(gr_numRefineVarsMax, 4, "Incorrect max refinement variables")
    call assertEqual(gr_numRefineVars, 3, "Incorrect max refinement variables")

    call assertEqual(gr_refine_var(1), DENS_VAR, "First refine var not DENS")
    call assertEqual(gr_refine_var(2), TEMP_VAR, "Second refine var not TEMP")
    call assertEqual(gr_refine_var(3), ENER_VAR, "Third refine var not ENER")

    call assertEqual(gr_refine_cutoff(1), 0.8d0, "Incorrect DENS refine cutoff")
    call assertEqual(gr_refine_cutoff(2), 0.5d0, "Incorrect TEMP refine cutoff")
    call assertEqual(gr_refine_cutoff(3), 0.6d0, "Incorrect ENER refine cutoff")

    call assertEqual(gr_derefine_cutoff(1), 0.45d0, &
                     "Incorrect DENS derefine cutoff")
    call assertEqual(gr_derefine_cutoff(2), 0.325d0, &
                     "Incorrect TEMP derefine cutoff")
    call assertEqual(gr_derefine_cutoff(3), 0.35d0, &
                     "Incorrect ENER derefine cutoff")

    call assertEqual(gr_refine_filter(1), 0.05d0, &
                     "Incorrect DENS derefine cutoff")
    call assertEqual(gr_refine_filter(2), 0.025d0, &
                     "Incorrect TEMP derefine cutoff")
    call assertEqual(gr_refine_filter(3), 0.035d0, &
                     "Incorrect ENER derefine cutoff")

    call assertFalse(gr_enforceMaxRefinement, "gr_enforceMaxRefinement True")

    !!!!! CONFIRM EoS SETUP
    call assertEqual(gr_eosMode, MODE_DENS_EI, &
                     "Incorrect eosMode")
    call assertEqual(gr_eosModeInit, MODE_DENS_TEMP, &
                     "Incorrect eosModeInit")

    call finish_test_run

end subroutine Driver_evolveFlash


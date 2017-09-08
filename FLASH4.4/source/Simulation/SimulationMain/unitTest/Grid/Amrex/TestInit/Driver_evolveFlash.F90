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
                                      Grid_getGeometry, &
                                      Grid_getDeltas, &
                                      Grid_getMaxRefinement, &
                                      Grid_getBlkPtr, Grid_releaseBlkPtr
    use Grid_data,             ONLY : gr_iguard, gr_jguard, gr_kguard, &
                                      gr_meshMe
    use block_iterator,        ONLY : block_iterator_t
    use block_metadata,        ONLY : block_metadata_t, bmd_print

    implicit none

    interface assertEqual
        procedure :: assertEqualInt
        procedure :: assertEqualReal
    end interface assertEqual

    ! These values should match those values in flash.par
    integer,  parameter :: NXCELL_EX   = 64
    integer,  parameter :: NYCELL_EX   = 64
    integer,  parameter :: NZCELL_EX   =  4
    ! DEVNOTE: FIXME Not able to configure different block numbers with octree
!    integer,  parameter :: NXBLK_EX    =  8
!    integer,  parameter :: NYBLK_EX    = 16
!    integer,  parameter :: NZBLK_EX    =  2
    integer,  parameter :: NXBLK_EX    = 16
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

    integer :: n_tests = 0
    integer :: n_failed = 0
    real    :: t_old = 0.0d0
    real    :: t_new = 0.0d0

    integer :: geometry = -100
    real    :: domain(LOW:HIGH, MDIM) = 0.0d0
    integer :: domainBC(LOW:HIGH, MDIM) = PERIODIC 
    real    :: deltas(1:MDIM) = 0.0d0
    real    :: x_expected = 0.0d0
    real    :: y_expected = 0.0d0
    real    :: z_expected = 0.0d0
    integer :: max_level = -1

    type(block_iterator_t) :: itor
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

    integer :: rank = -1
    integer :: ilev = 0
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    integer :: var = 0

    rank = amrex_parallel_myproc()
 
    write(*,*)
    call cpu_time(t_old)

    !!!!! CONFIRM PROPER COORDINATE SYSTEM
    ! Dimensionality
    if (amrex_spacedim /= NDIM) then
        write(*,*) "Wrong dimensionality - ", amrex_spacedim, ' != ', NDIM
        write(*,*) "Recompile AMReX with correct dimensionality"
        write(*,*)
        stop
    end if
    n_tests = n_tests + 1
 
    ! Physical domain
    call Grid_getGeometry(geometry)
    call assertEqual(geometry, CARTESIAN, "Incorrect coordinate system type")

    call Grid_getDomainBoundBox(domain)
#if NDIM == 1
    call assertEqual(domain(LOW,  1), XMIN_EX, "Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, 1), XMAX_EX, "Incorrect high X-coordinate")
    call assertEqual(domain(LOW,  2), 0.0d0,   "Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, 2), 0.0d0,   "Incorrect high Y-coordinate")
    call assertEqual(domain(LOW,  3), 0.0d0,   "Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, 3), 0.0d0,   "Incorrect high Z-coordinate")
#elif NDIM == 2
    call assertEqual(domain(LOW,  1), XMIN_EX, "Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, 1), XMAX_EX, "Incorrect high X-coordinate")
    call assertEqual(domain(LOW,  2), YMIN_EX, "Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, 2), YMAX_EX, "Incorrect high Y-coordinate")
    call assertEqual(domain(LOW,  3), 0.0d0,   "Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, 3), 0.0d0,   "Incorrect high Z-coordinate")
#elif NDIM == 3 
    call assertEqual(domain(LOW,  1), XMIN_EX, "Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, 1), XMAX_EX, "Incorrect high X-coordinate")
    call assertEqual(domain(LOW,  2), YMIN_EX, "Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, 2), YMAX_EX, "Incorrect high Y-coordinate")
    call assertEqual(domain(LOW,  3), ZMIN_EX, "Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, 3), ZMAX_EX, "Incorrect high Z-coordinate")
#endif

    !!!!! CONFIRM PROPER BC
    call Grid_getDomainBC(domainBC)
    call assertEqual(domainBC(LOW,  IAXIS), XL_BC_EX, "Incorrect X-left BC")
    call assertEqual(domainBC(HIGH, IAXIS), XH_BC_EX, "Incorrect X-right BC")
    call assertEqual(domainBC(LOW,  JAXIS), YL_BC_EX, "Incorrect Y-left BC")
    call assertEqual(domainBC(HIGH, JAXIS), YH_BC_EX, "Incorrect Y-right BC")
    call assertEqual(domainBC(LOW,  KAXIS), ZL_BC_EX, "Incorrect Z-left BC")
    call assertEqual(domainBC(HIGH, KAXIS), ZH_BC_EX, "Incorrect Z-right BC")

    !!!!! CONFIRM PROPER REFINEMENT
    call Grid_getMaxRefinement(max_level, mode=1)
    call assertEqual(max_level, MAXLEVEL_EX, "Incorrect maximum refine level")

    do ilev = 1, max_level
        call Grid_getDeltas(ilev, deltas)

        x_expected = XDELTA_EX / 2.0d0**(ilev - 1)
        y_expected = YDELTA_EX / 2.0d0**(ilev - 1)
        z_expected = ZDELTA_EX / 2.0d0**(ilev - 1)
#if NDIM == 1
        call assertEqual(deltas(1), x_expected, "Incorrect high X-coordinate")
        call assertEqual(deltas(2), 0.0d0,      "Incorrect high Y-coordinate")
        call assertEqual(deltas(3), 0.0d0,      "Incorrect high Z-coordinate")
#elif NDIM == 2
        call assertEqual(deltas(1), x_expected, "Incorrect high X-coordinate")
        call assertEqual(deltas(2), y_expected, "Incorrect high Y-coordinate")
        call assertEqual(deltas(3), 0.0d0,      "Incorrect high Z-coordinate")
#elif NDIM == 3 
        call assertEqual(deltas(1), x_expected, "Incorrect high X-coordinate")
        call assertEqual(deltas(2), y_expected, "Incorrect high Y-coordinate")
        call assertEqual(deltas(3), z_expected, "Incorrect high Z-coordinate")
#endif
    end do

    ! TODO: Once unittest is refining mesh, check Grid_getMaxRefinement
    ! with mode that checks actual number of levels in use

    !!!!! CONFIRM PROPER BLOCK/CELL STRUCTURE
    ! Walk across all blocks to test and collect info
    n_blocks = 0
    itor = block_iterator_t(LEAF)
    call itor%blkMetaData(block)
    xBlkMin = block%limits(LOW, 1)
    xBlkMax = block%limits(HIGH, 1)
    yBlkMin = block%limits(LOW, 2)
    yBlkMax = block%limits(HIGH, 2)
    zBlkMin = block%limits(LOW, 3)
    zBlkMax = block%limits(HIGH, 3)
    do while (itor%is_valid())
        n_blocks = n_blocks + 1
        call itor%blkMetaData(block)

        ! DEVNOTE: Should we leave this unittest with simple data
        ! that does not refine so that testing the block structure is easy?
        ! All blocks on coarsest level since no refining
        call assertEqual(block%level, 1, "Incorrect block level")

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

    ! Confirm proper number of blocks and cells
    ! FIXME: This only works for one processor.  Need a reduction here.
#if NDIM == 1
    call assertEqual(n_blocks, NXBLK_EX, &
                     "Incorrect total number of blocks")

    call assertEqual(xBlkMin, 0, "Incorrect origin X-coordinate")
    ! DEVNOTE: FIXME Shouldn't these be zero as well?
    call assertEqual(yBlkMin, 1, "Incorrect origin Y-coordinate")
    call assertEqual(zBlkMin, 1, "Incorrect origin Z-coordinate")

    call assertEqual(xBlkMax + 1, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    call assertEqual(yBlkMax, 1, "More than one cell along Y-axis")
    call assertEqual(zBlkMax, 1, "More than one cell along Z-axis")
#elif NDIM == 2
    call assertEqual(n_blocks, NXBLK_EX*NYBLK_EX, &
                     "Incorrect total number of blocks")

    call assertEqual(xBlkMin, 0, "Incorrect origin X-coordinate")
    call assertEqual(yBlkMin, 0, "Incorrect origin Y-coordinate")
    ! DEVNOTE: FIXME Shouldn't these be zero as well?
    call assertEqual(zBlkMin, 1, "Incorrect origin Z-coordinate")

    call assertEqual(xBlkMax + 1, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    call assertEqual(yBlkMax + 1, NYCELL_EX, &
                     "Incorrect total number of cells along Y-axis")
    call assertEqual(zBlkMax, 1, "More than one cell along Z-axis")
#elif NDIM == 3
    call assertEqual(n_blocks, NXBLK_EX*NYBLK_EX*NZBLK_EX, &
                     "Incorrect total number of blocks")

    call assertEqual(xBlkMin, 0, "Incorrect origin X-coordinate")
    call assertEqual(yBlkMin, 0, "Incorrect origin Y-coordinate")
    call assertEqual(zBlkMin, 0, "Incorrect origin Z-coordinate")

    call assertEqual(xBlkMax + 1, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    call assertEqual(yBlkMax + 1, NYCELL_EX, &
                     "Incorrect total number of cells along Y-axis")
    call assertEqual(zBlkMax + 1, NZCELL_EX, &
                     "Incorrect total number of cells along Z-axis")
#endif

    !!!!! CONFIRM REFINEMENT SETUP
    ! TODO: Get nrefs from AMReX

    !!!!! CONFIRM MPI SETUP
    call assertEqual(rank, gr_meshMe, "AMReX/FLASH ranks are different")

    !!!!! CONFIRM PROPER INITIAL CONDITIONS
    itor = block_iterator_t(LEAF)
    do while (itor%is_valid())
        call itor%blkMetaData(block)
        call Grid_getBlkPtr(block, solnData)

        associate(lo => block%limits(LOW, :), &
                  hi => block%limits(HIGH, :))
            do         k = lo(3), hi(3)
                do     j = lo(2), hi(2)
                    do i = lo(1), hi(1)
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
 
    !!!!! OUTPUT RESULTS
    ! DEVNOTE: reduction to collect number of fails?
    if (rank == MASTER_PE) then
        write(*,*)
        if (n_failed == 0) then
            write(*,*) "SUCCESS - ", &
                       (n_tests - n_failed), "/", n_tests, ' passed' 
        else 
            write(*,*) "FAILURE - ", &
                       (n_tests - n_failed), "/", n_tests, ' passed'
        end if
        write(*,*)
        write(*,*) 'Walltime = ', (t_new - t_old), ' s'
        write(*,*)
    end if

contains

    subroutine assertEqualInt(a, b, msg)
        implicit none

        integer,      intent(IN) :: a
        integer,      intent(IN) :: b
        character(*), intent(IN) :: msg

        character(256) :: buffer = ""

        if (a /= b) then
            write(buffer,'(A,I5,A,I5)') msg, a, " != ", b
            write(*,*) TRIM(ADJUSTL(buffer))
            n_failed = n_failed + 1
        end if
        n_tests = n_tests + 1
    end subroutine assertEqualInt

    subroutine assertEqualReal(a, b, msg)
        implicit none

        real,         intent(IN) :: a
        real,         intent(IN) :: b
        character(*), intent(IN) :: msg

        character(256) :: buffer = ""

        if (a /= b) then
            write(buffer,'(A,F15.8,A,F15.8)') msg, a, " != ", b
            write(*,*) TRIM(ADJUSTL(buffer))
            n_failed = n_failed + 1
        end if
        n_tests = n_tests + 1
    end subroutine assertEqualReal

    subroutine assert_almost_equal(a, b, prec, msg)
        implicit none

        real,         intent(IN) :: a
        real,         intent(IN) :: b
        real,         intent(IN) :: prec
        character(*), intent(IN) :: msg

        character(256) :: buffer = ""

        if (ABS(b - a) > prec) then
            write(buffer,'(A,F15.8,A,F15.8)') msg, a, " != ", b
            write(*,*) TRIM(ADJUSTL(buffer))
            n_failed = n_failed + 1
        end if
        n_tests = n_tests + 1
    end subroutine assert_almost_equal

end subroutine Driver_evolveFlash


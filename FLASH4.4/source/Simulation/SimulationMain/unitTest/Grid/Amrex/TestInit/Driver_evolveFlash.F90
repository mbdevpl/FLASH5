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
                                      Grid_getMaxRefinement
    use Grid_data,             ONLY : gr_nBlockX, gr_nBlockY, gr_nBlockZ, &
                                      gr_meshMe

    implicit none

    interface assertEqual
        procedure :: assertEqualInt
        procedure :: assertEqualReal
    end interface assertEqual

    ! These values should match those values in flash.par
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
    integer,  parameter :: MAXLEVEL_EX = 4

    integer :: n_tests = 0
    integer :: n_failed = 0
    real    :: t_old = 0.0d0
    real    :: t_new = 0.0d0

    integer :: geometry = -100
    real    :: domain(LOW:HIGH, MDIM) = 0.0d0
    real    :: deltas(1:MDIM) = 0.0d0
    real    :: x_expected = 0.0d0
    real    :: y_expected = 0.0d0
    real    :: z_expected = 0.0d0
    integer :: max_level = -1

    integer :: rank = -1
    integer :: ilev = 0

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

    !!!!! CONFIRM PROPER CELL/BLOCK STRUCTURE
    ! Cells
    ! TODO: Get values from AMReX

    ! Refinement levels
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

    ! Blocks
    ! TODO: Get values from AMReX rather than from FLASH
    call assertEqual(gr_nBlockX, NXBLK_EX, "Wrong number of blocks along X")
    call assertEqual(gr_nBlockY, NYBLK_EX, "Wrong number of blocks along Y")
    call assertEqual(gr_nBlockZ, NZBLK_EX, "Wrong number of blocks along Z")

    !!!!! CONFIRM REFINEMENT SETUP
    ! TODO: Get nrefs from AMReX

    !!!!! CONFIRM MPI SETUP
    call assertEqual(rank, gr_meshMe, "AMReX/FLASH ranks are different")

    ! TODO: Once unittest is refining mesh, check Grid_getMaxRefinement
    ! with mode that checks actual number of levels in use

    call cpu_time (t_new)
 
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


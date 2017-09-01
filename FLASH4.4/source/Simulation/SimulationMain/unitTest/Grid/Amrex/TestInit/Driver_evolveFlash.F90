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
!!
!!  Simple stripped down version for testing single units.
!!
!! NOTES
!!  This simulation *must* be configured with at least the following
!!  2D run:
!!     ./setup -auto -2d -nxb=8 -nyb=4 
!!              unitTest/Grid/Amrex/TestInit 
!!             -unit=IO/IOMain/hdf5/serial/AM
!!  3D run:
!!     ./setup -auto -3d -nxb=8 -nyb=4 -nzb=2
!!              unitTest/Grid/Amrex/TestInit 
!!             -unit=IO/IOMain/hdf5/serial/AM
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Driver_evolveFlash()
    use amrex_fort_module,     ONLY : wp => amrex_real, &
                                      amrex_spacedim
    use amrex_amrcore_module,  ONLY : amrex_get_finest_level, &
                                      amrex_get_boxarray
    use amrex_box_module,      ONLY : amrex_print, amrex_box
    use amrex_parallel_module, ONLY : amrex_parallel_myproc, &
                                      amrex_parallel_nprocs

    use Grid_interface,        ONLY : Grid_getDomainBoundBox
    use Grid_data,             ONLY : gr_nBlockX, gr_nBlockY, gr_nBlockZ, &
                                      gr_meshMe

    implicit none
 
    ! These values should match those values in flash.par
    integer,  parameter :: NXCELL_EX  = 64
    integer,  parameter :: NYCELL_EX  = 64
    integer,  parameter :: NZCELL_EX  =  4
    integer,  parameter :: NXBLK_EX   =  8
    integer,  parameter :: NYBLK_EX   = 16
    integer,  parameter :: NZBLK_EX   =  2
    real(wp), parameter :: XMIN_EX    = -1.00_wp
    real(wp), parameter :: XMAX_EX    =  2.00_wp
    real(wp), parameter :: YMIN_EX    = -1.50_wp
    real(wp), parameter :: YMAX_EX    =  4.50_wp
    real(wp), parameter :: ZMIN_EX    =  0.50_wp
    real(wp), parameter :: ZMAX_EX    =  0.75_wp

    character(256) :: msg = ""
    integer        :: n_tests = 0
    integer        :: n_failed = 0
    real(wp)       :: t_old = 0.0_wp
    real(wp)       :: t_new = 0.0_wp

    real           :: domain(LOW:HIGH, MDIM)

    integer  :: rank = -1

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
 
    ! TODO: Check coordinate system type from AMReX

    ! Physical domain
    call Grid_getDomainBoundBox(domain)
    call assert_almost_equal(domain(LOW, 1), XMIN_EX, 0.0_wp, &
                             "Incorrect low X-coordinate")
    call assert_almost_equal(domain(HIGH, 1), XMAX_EX, 0.0_wp, &
                             "Incorrect high X-coordinate")
#if NDIM >= 2
    call assert_almost_equal(domain(LOW, 2), YMIN_EX, 0.0_wp, &
                             "Incorrect low Y-coordinate")
    call assert_almost_equal(domain(HIGH, 2), YMAX_EX, 0.0_wp, &
                             "Incorrect high Y-coordinate")
#if NDIM == 3 
    call assert_almost_equal(domain(LOW, 3), ZMIN_EX, 0.0_wp, &
                             "Incorrect low Z-coordinate")
    call assert_almost_equal(domain(HIGH,3), ZMAX_EX, 0.0_wp, &
                             "Incorrect high Z-coordinate")
#endif
#endif

    !!!!! CONFIRM PROPER CELL/BLOCK STRUCTURE
    ! Cells
    ! TODO: Get values from AMReX

    ! Deltas
    ! TODO: Get values from AMReX

    ! Blocks
    ! TODO: Get values from AMReX rather than from FLASH
    call assert_equal(gr_nBlockX, NXBLK_EX, "Wrong number of blocks along X")
    call assert_equal(gr_nBlockY, NYBLK_EX, "Wrong number of blocks along Y")
    call assert_equal(gr_nBlockZ, NZBLK_EX, "Wrong number of blocks along Z")

    !!!!! CONFIRM REFINEMENT SETUP
    ! TODO: Get max refinement level from AMReX
    ! TODO: Get nrefs from AMReX

    !!!!! CONFIRM MPI SETUP
    call assert_equal(rank, gr_meshMe, "AMReX/FLASH ranks are different")

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

    subroutine assert_equal(a, b, msg)
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
    end subroutine assert_equal

    subroutine assert_almost_equal(a, b, prec, msg)
        implicit none

        real(wp),     intent(IN) :: a
        real(wp),     intent(IN) :: b
        real(wp),     intent(IN) :: prec
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


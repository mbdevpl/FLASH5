!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestRefine/Driver_evolveFlash
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
!!     ./setup -auto -2d -nxb=8 -nyb=8
!!              unitTest/Grid/Amrex/TestRefine
!!             +noio -index-reorder
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Driver_evolveFlash()
    use amrex_fort_module,     ONLY : amrex_spacedim
    use amrex_amrcore_module,  ONLY : amrex_get_finest_level

    use Grid_interface,        ONLY : Grid_getDomainBoundBox, &
                                      Grid_getDeltas, &
                                      Grid_updateRefinement
    use Grid_data,             ONLY : gr_iguard, gr_jguard, gr_kguard, &
                                      gr_meshMe
    use block_iterator,        ONLY : block_iterator_t
    use block_metadata,        ONLY : block_metadata_t, bmd_print

    implicit none

    interface assertEqual
        procedure :: assertEqualInt
        procedure :: assertEqualReal
    end interface assertEqual

    integer :: n_tests
    integer :: n_failed
    real    :: t_old
    real    :: t_new

    real :: domain(LOW:HIGH, 1:MDIM)
    real :: deltas(1:MDIM)

    type(block_iterator_t) :: itor
    type(block_metadata_t) :: block

    integer :: level
    integer :: finest_level

    write(*,*)
    n_tests = 0
    n_failed = 0
    t_old = 0.0d0
    t_new = 0.0d0
    call cpu_time(t_old)

    !!!!! CONFIRM PROPER DIIMENSIONALITY
    if (amrex_spacedim /= 2) then
        write(*,*) "Wrong dimensionality - ", amrex_spacedim, ' != ', 2
        write(*,*) "Recompile AMReX with correct dimensionality"
        write(*,*)
        stop
    end if
    n_tests = n_tests + 1
 
    !!!!! CONFIRM PROPER SETUP
    ! 16x16 domain with dx = dy = (1.0 - 0.0)/16
    call assertEqual(NXB, 8, "Incorrect initial number of cells/block along X")
    call assertEqual(NYB, 8, "Incorrect initial number of cells/block along Y")
    call assertEqual(NZB, 1, "Incorrect initial number of cells/block along Z")

    domain = 0.0d0
    call Grid_getDomainBoundBox(domain)
    call assertEqual(domain(LOW,  IAXIS), 0.0d0, "Incorrect Xi-coord")
    call assertEqual(domain(LOW,  JAXIS), 0.0d0, "Incorrect Yi-coord")
    call assertEqual(domain(LOW,  KAXIS), 0.0d0, "Incorrect Zi-coord")
    call assertEqual(domain(HIGH, IAXIS), 1.0d0, "Incorrect Xf-coord")
    call assertEqual(domain(HIGH, JAXIS), 1.0d0, "Incorrect Yf-coord")
    call assertEqual(domain(HIGH, KAXIS), 0.0d0, "Incorrect Zf-coord")

    deltas = 0.0d0
    call Grid_getDeltas(1, deltas)
    call assertEqual(deltas(IAXIS), 0.0625d0, "Invalid dx at coarse level")
    call assertEqual(deltas(JAXIS), deltas(IAXIS), "dy != dx at coarse")
    call assertEqual(deltas(KAXIS), 0.0d0, "dz != 0 at coarse")

    write(*,*) "LEAF OUTPUT SECTION"
    write(*,*) "-------------------------------------------------------"
    finest_level = amrex_get_finest_level() + 1
    do level = 1, finest_level 
        write(*,*) "LEVEL = ", level
        itor = block_iterator_t(LEAF, level=level)
        do while (itor%is_valid())
            call itor%blkMetaData(block)
            write(*,*) block%limits(LOW, :)
            write(*,*) block%limits(HIGH, :)
            write(*,*) "---------"

            call itor%next()
        end do
    end do

    !!!!! CONFIRM INITIAL REFINEMENT
    ! Started with 2x2 block structure and refined according to initial data
    ! using unittests own gr_markRefineDerefine callback with AMReX

    !!!!! OUTPUT RESULTS
    ! DEVNOTE: reduction to collect number of fails?
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


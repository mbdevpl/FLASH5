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
!!  3D run:
!!     ./setup -auto -3d -nxb=4 -nyb=4 -nzb=4
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
    
    use Grid_interface,        ONLY : Grid_fillGuardCells
    use amrex_interfaces,      ONLY : gr_getFinestLevel
    use gr_physicalMultifabs,  ONLY : unk
    use block_metadata,        ONLY : block_metadata_t

    implicit none

    interface assertEqual
        procedure :: assertEqualInt
        procedure :: assertEqualReal
    end interface assertEqual

    integer, parameter :: FINEST_LEVEL_EX = 2

    integer :: n_tests
    integer :: n_failed
    real    :: t_old
    real    :: t_new

    integer :: rank
    integer :: finest_level
    integer :: lev
    integer :: i, j, k
    integer :: i_bc, j_bc, k_bc
    integer :: var
    integer :: bdd

    type(amrex_mfiter)     :: mfi
    type(amrex_box)        :: bx
    type(block_metadata_t) :: blockDesc

    real, contiguous, pointer :: initData(:, :, :, :)

    n_tests = 0
    n_failed = 0
    t_old = 0.0d0
    t_new = 0.0d0

    rank = amrex_parallel_myproc()
 
    write(*,*)
    call cpu_time(t_old)

    !!!!! CONFIRM PROPER COORDINATE SYSTEM DIMENSIONALITY
    if (amrex_spacedim /= NDIM) then
        write(*,*) "Wrong dimensionality - ", amrex_spacedim, ' != ', NDIM
        write(*,*) "Recompile AMReX with correct dimensionality"
        write(*,*)
        stop
    end if
    n_tests = n_tests + 1

    !!!!! CONFIRM THAT GC ARE ZERO AFTER INIT 
    call gr_getFinestLevel(finest_level)
    call assertEqual(finest_level, FINEST_LEVEL_EX, "Incorrect finest level")

    do lev = 1, finest_level
        call amrex_mfiter_build(mfi, unk(lev-1), tiling=.FALSE.)
        do while (mfi%next())
            bx = mfi%tilebox()

            ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
            ! Level must be 1-based index and limits/limitsGC must be 1-based also
            ! DEVNOTE: Should we use gr_[ijk]guard here?
            blockDesc%level = lev
            ! DEVNOTE: TODO Get grid_index from mfi
            blockDesc%grid_index = -1
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
                ! DEVNOT: TODO Use Grid_getBlkPtr once we have correct grid_index
                initData(loGC(1):, loGC(2):, loGC(3):, 1:) => unk(lev-1)%dataptr(mfi)

                do     k = loGC(KAXIS), hiGC(KAXIS)
                  do   j = loGC(JAXIS), hiGC(JAXIS)
                    do i = loGC(IAXIS), hiGC(IAXIS)
                      do var=UNK_VARS_BEGIN, UNK_VARS_END
                        if (      (lo(IAXIS) <= i) .AND. (i <= hi(IAXIS)) &
                            .AND. (lo(JAXIS) <= j) .AND. (j <= hi(JAXIS)) &
                            .AND. (lo(KAXIS) <= k) .AND. (k <= hi(KAXIS))) then
                            ! Interior has data
                            call assertEqual(initData(i, j, k, var), &
                                             DBLE((i + j + k) * var), &
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
            ! DEVNOTE: Should we use gr_[ijk]guard here?
            blockDesc%level = lev
            ! DEVNOTE: TODO Get grid_index from mfi
            blockDesc%grid_index = -1
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
                ! DEVNOT: TODO Use Grid_getBlkPtr once we have correct grid_index
                initData(loGC(1):, loGC(2):, loGC(3):, 1:) => unk(lev-1)%dataptr(mfi)

                bdd = 8 * 2**(lev - 1)
                do     k = loGC(KAXIS), hiGC(KAXIS)
                  do   j = loGC(JAXIS), hiGC(JAXIS)
                    do i = loGC(IAXIS), hiGC(IAXIS)
                      do var=UNK_VARS_BEGIN, UNK_VARS_END
                        ! Account for periodic BC
                        i_bc = i
                        j_bc = j
                        k_bc = k
                        if (i_bc > bdd) then
                            i_bc = i_bc - bdd 
                        else if (i_bc < 1) then
                            i_bc = i_bc + bdd 
                        end if

                        if (j_bc > bdd) then
                            j_bc = j_bc - bdd 
                        else if (j_bc < 1) then
                            j_bc = j_bc + bdd 
                        end if

                        if (k_bc > bdd) then
                            k_bc = k_bc - bdd
                        else if (k_bc < 1) then
                            k_bc = k_bc + bdd 
                        end if
                        
                        if (      (lo(IAXIS) <= i) .AND. (i <= hi(IAXIS)) &
                            .AND. (lo(JAXIS) <= j) .AND. (j <= hi(JAXIS)) &
                            .AND. (lo(KAXIS) <= k) .AND. (k <= hi(KAXIS))) then
                            ! Interior has data
                            call assertEqual(initData(i, j, k, var), &
                                             DBLE((i + j + k) * var), &
                                             "Bad data")
                        else
                          ! Spot check a few values for now
                          ! DEV: TODO Improve this check or sufficient?
                          if (lev == 2) then
                            if((i == -1) .AND. (j == -1)) then
                              call assertEqual(initData(i, j, k, var), &
                                               17.0d0*var, &
                                               "Incorrect fine GC value")
                            else if ((i == 0) .AND. (j == 5)) then
                              call assertEqual(initData(i, j, k, var), &
                                               11.75d0*var, &
                                               "Incorrect fine GC value")
                            else if ((i == 0) .AND. (j == 6)) then
                              call assertEqual(initData(i, j, k, var), &
                                               12.25d0*var, &
                                               "Incorrect fine GC value")
                            else if ((i == 3) .AND. (j == 9)) then
                              call assertEqual(initData(i, j, k, var), &
                                               7.5d0*var, &
                                               "Incorrect fine GC value")
                            else if ((i == 3) .AND. (j == 10)) then
                              call assertEqual(initData(i, j, k, var), &
                                               8.0d0*var, &
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

    subroutine assertTrue(a, msg)
        implicit none

        logical,      intent(IN) :: a
        character(*), intent(IN) :: msg

        character(256) :: buffer = ""
        
        if (.NOT. a) then
            write(buffer,'(A)') msg
            write(*,*) TRIM(ADJUSTL(buffer))
            n_failed = n_failed + 1
        end if
        n_tests = n_tests + 1
    end subroutine assertTrue

    subroutine assertFalse(a, msg)
        implicit none

        logical,      intent(IN) :: a
        character(*), intent(IN) :: msg

        character(256) :: buffer = ""
        
        if (a) then
            write(buffer,'(A)') msg
            write(*,*) TRIM(ADJUSTL(buffer))
            n_failed = n_failed + 1
        end if
        n_tests = n_tests + 1
    end subroutine assertFalse

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


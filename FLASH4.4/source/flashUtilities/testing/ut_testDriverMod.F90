!!****if* source/flashUtilities/testing/ut_testDriverMod
!!
!! NAME
!!
!!  ut_testDriverMod
!!
!! SYNOPSIS
!!
!!  use ut_testDriverMod
!!
!! DESCRIPTION
!!
!!  A module that encapsulates common routines and variables that can be used in
!!  unittests.  Typically,  this module is used in the unittests version of 
!!  Driver_evolveFlash as follows
!!
!!      <declare variables and setup for test>
!!      ...
!!      call start_test_run
!!      ...
!!      call assertFalse(didAbort, "Action aborted")
!!      call assertAlmostEqual(x, 1.1, 1.0d-15, "Incorrect x value")
!!      ...
!!      call finish_test_run
!!      ...
!!      <clean-up>
!!
!! NOTES
!!  
!!  The finish_test_run prints results to standard output.
!!
!!***

#include "constants.h"

module ut_testDriverMod
    implicit none
    private

    integer, save :: my_n_tests = 0
    integer, save :: my_n_failed = 0
    real,    save :: my_t_start = 0.0d0
    logical, save :: is_testing = .FALSE.

    interface assertEqual
        procedure :: assertEqualInt
        procedure :: assertEqualReal
    end interface assertEqual

    interface assertSetEqual
        procedure :: assertSetEqual2dIntArray
    end interface assertSetEqual

    public :: start_test_run
    public :: finish_test_run

    public :: assertTrue
    public :: assertFalse
    public :: assertEqual
    public :: assertSetEqual
    public :: assertAlmostEqual

contains

    subroutine start_test_run()
        use Driver_Interface, ONLY : Driver_abortFlash

        if (is_testing) then
            call Driver_abortFlash("[start_test_run] Already testing")
        end if

        is_testing = .TRUE.

        my_t_start = 0.0d0
        call cpu_time(my_t_start)
    end subroutine start_test_run

    subroutine finish_test_run
        use Grid_data,        ONLY : gr_meshMe
        use Driver_Interface, ONLY : Driver_abortFlash

        real :: my_t_end
        real :: my_walltime

        if (.NOT. is_testing) then
            call Driver_abortFlash("[finish_test_run] Not testing yet")
        end if

        is_testing = .FALSE.
        
        call cpu_time(my_t_end)
        my_walltime = my_t_end - my_t_start

        ! DEV: TODO reduction to collect number of tests/fails/max walltime?
        if (gr_meshMe == MASTER_PE) then
            write(*,*)
            if (my_n_failed == 0) then
                write(*,*) "SUCCESS - ", &
                           (my_n_tests - my_n_failed), "/", my_n_tests, ' passed' 
            else 
                write(*,*) "FAILURE - ", &
                           (my_n_tests - my_n_failed), "/", my_n_tests, ' passed'
            end if
            write(*,*)
            write(*,*) 'Walltime = ', my_walltime, ' s'
            write(*,*)
        end if
    end subroutine finish_test_run

    subroutine assertTrue(a, msg)
        logical,      intent(IN) :: a
        character(*), intent(IN) :: msg

        character(256) :: buffer = ""
        
        if (.NOT. a) then
            write(buffer,'(A)') msg
            write(*,*) TRIM(ADJUSTL(buffer))
            my_n_failed = my_n_failed + 1
        end if
        my_n_tests = my_n_tests + 1
    end subroutine assertTrue

    subroutine assertFalse(a, msg)
        logical,      intent(IN) :: a
        character(*), intent(IN) :: msg

        character(256) :: buffer = ""
        
        if (a) then
            write(buffer,'(A)') msg
            write(*,*) TRIM(ADJUSTL(buffer))
            my_n_failed = my_n_failed + 1
        end if
        my_n_tests = my_n_tests + 1
    end subroutine assertFalse

    subroutine assertEqualInt(a, b, msg)
        integer,      intent(IN) :: a
        integer,      intent(IN) :: b
        character(*), intent(IN) :: msg

        character(256) :: buffer = ""

        if (a /= b) then
            write(buffer,'(A,I5,A,I5)') msg, a, " != ", b
            write(*,*) TRIM(ADJUSTL(buffer))
            my_n_failed = my_n_failed + 1
        end if
        my_n_tests = my_n_tests + 1
    end subroutine assertEqualInt

    subroutine assertEqualReal(a, b, msg)
        real,         intent(IN) :: a
        real,         intent(IN) :: b
        character(*), intent(IN) :: msg

        character(256) :: buffer = ""

        if (a /= b) then
            write(buffer,'(A,F15.8,A,F15.8)') msg, a, " != ", b
            write(*,*) TRIM(ADJUSTL(buffer))
            my_n_failed = my_n_failed + 1
        end if
        my_n_tests = my_n_tests + 1
    end subroutine assertEqualReal

    subroutine assertAlmostEqual(a, b, prec, msg)
        real,         intent(IN) :: a
        real,         intent(IN) :: b
        real,         intent(IN) :: prec
        character(*), intent(IN) :: msg

        character(256) :: buffer = ""

        if (ABS(b - a) > prec) then
            write(buffer,'(A,F15.8,A,F15.8)') msg, a, " != ", b
            write(*,*) TRIM(ADJUSTL(buffer))
            my_n_failed = my_n_failed + 1
        end if
        my_n_tests = my_n_tests + 1
    end subroutine assertAlmostEqual

    subroutine assertSetEqual2dIntArray(A, B, msg)
        integer,      intent(IN) :: A(:, :)
        integer,      intent(IN) :: B(:, :)
        character(*), intent(IN) :: msg

        logical        :: in_set
        logical        :: failed
        integer        :: j, k

        my_n_tests = my_n_tests + 1

        ! Confirm A subset of B
        failed = .FALSE.
        do j = 1, SIZE(A, 1)
            in_set = .FALSE.
            do k = 1, SIZE(B, 1)
                if (ALL(A(j, :) == B(k, :))) then
                    in_set = .TRUE.
                    exit
                end if
            end do

            if (.NOT. in_set) then
                write(*,*) msg, " - ", A(j, :), " of A not in B"
                failed = .TRUE.
            end if
        end do

        ! Confirm B subset of A
        do j = 1, SIZE(B, 1)
            in_set = .FALSE.
            do k = 1, SIZE(A, 1)
                if (ALL(B(j, :) == A(k, :))) then
                    in_set = .TRUE.
                    exit
                end if
            end do

            if (.NOT. in_set) then
                write(*,*) msg, " - ", B(j, :), " of B not in A"
                failed = .TRUE.
            end if
        end do

        if (failed) then
            my_n_failed = my_n_failed + 1
        end if
    end subroutine assertSetEqual2dIntArray

end module ut_testDriverMod


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
!!  The driver for a toy version of a full FLASH simulation that allows users
!!  to control and therefore experiment with how AMReX mesh refinement and
!!  derefinement is working as implemented in FLASH.  The code advances the data
!!  in UNK in time manually.  At each step, the code sets all data in UNK to
!!  zero except for possibly at a few points, whose non-zero values define
!!  what level of refinement must be achieved in the blocks that contain them.  
!!
!!  This simulation serves as a form for manually testing appropriate refinement
!!  with AMReX (See documentation in folder).  Ideally, the leaf blocks will be
!!  automatically verified at each step.
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
#include "sim_constants.h"

subroutine Driver_evolveFlash()
    use amrex_fort_module,     ONLY : amrex_spacedim
    use amrex_box_module,      ONLY : amrex_box
    use amrex_multifab_module, ONLY : amrex_mfiter, &
                                      amrex_mfiter_build, &
                                      amrex_mfiter_destroy

    use Grid_interface,        ONLY : Grid_getDomainBoundBox, &
                                      Grid_getDeltas, &
                                      Grid_updateRefinement, &
                                      Grid_getBlkPtr, Grid_releaseBlkPtr
    use Grid_data,             ONLY : gr_meshMe, gr_lRefineMax, gr_maxRefine
    use amrex_interfaces,      ONLY : gr_getFinestLevel
    use sim_interface,         ONLY : sim_advance
    use block_metadata,        ONLY : block_metadata_t
    use gr_physicalMultifabs,  ONLY : unk

    implicit none

    interface assertEqual
        procedure :: assertEqualInt
        procedure :: assertEqualReal
    end interface assertEqual

    integer :: n_tests
    integer :: n_failed
    real    :: t_old
    real    :: t_new

    real    :: domain(LOW:HIGH, 1:MDIM)
    real    :: deltas(1:MDIM)
    integer :: finest_level

    real :: points(3, 1:NDIM)
    real :: values(3)

    ! DEV: TODO Get rid of hardcoded max levels here and elsewhere
    integer :: block_count(4)
    integer :: block_count_ex(4)
    integer :: lev, i, j

    real, contiguous, pointer :: solnData(:,:,:,:)
    type(block_metadata_t)    :: blockDesc
    type(amrex_mfiter)        :: mfi
    type(amrex_box)           :: bx

    write(*,*)
    n_tests = 0
    n_failed = 0
    t_old = 0.0d0
    t_new = 0.0d0
    call cpu_time(t_old)

    !!!!! CONFIRM PROPER DIMENSIONALITY
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

    call assertEqual(4, gr_lRefineMax, "Incorrect max number of levels")
    call assertEqual(gr_maxRefine, gr_lRefineMax, "gr_maxRefine != gr_lRefineMax")

    call sim_printLeaves("LEAVES AFTER DATA INIT & REGRID", block_count)
    block_count_ex = [0, 15, 4, 0]
    do lev = 1, SIZE(block_count)
        call assertEqual(block_count(lev), block_count_ex(lev), "Wrong # of levels")
    end do

    !!!!! CONFIRM INITIAL REFINEMENT
    ! Started with 2x2 block structure and refined according to initial data
    ! using unittests own gr_markRefineDerefine callback with AMReX
    call gr_getFinestLevel(finest_level)
    call assertEqual(3, finest_level, "Incorrect finest level after init")

    !!!!! STEP 1/2 - CONFIRM DEREFINEMENT GLOBALLY TO LEVEL 1
    points(:, :) = 0.0d0
    values(:) = 0.0d0
    call sim_advance(1, points, values, &
                     "SETTING ALL DATA TO ZERO AT ALL LEVELS", &
                     "LEAVES AFTER ZEROING ALL DATA & REGRID", &
                     block_count)
    call gr_getFinestLevel(finest_level)
    call assertEqual(1, finest_level, "Incorrect finest level")
    
    block_count_ex = [4, 0, 0, 0]
    do lev = 1, SIZE(block_count)
        call assertEqual(block_count(lev), block_count_ex(lev), "Wrong # of levels")
    end do

    !!!!! STEP 3/4 - CONFIRM REFINEMENT GLOBALLY TO LEVEL 2
    ! Corner cell with periodic BC
    points(:, :) = 0.0d0
    values(:) = 0.0d0
    points(1, :) = [0.99d0, 0.01d0]
    values(1) = REFINE_TO_L2
    call sim_advance(3, points, values, &
                     "SETTING CORNER CELL ONLY FOR LEVEL 2", &
                     "LEAVES AFTER DATA AT CORNER CELL", &
                     block_count)
    call gr_getFinestLevel(finest_level)
    call assertEqual(2, finest_level, "Incorrect finest level")

    block_count_ex = [0, 16, 0, 0]
    do lev = 1, SIZE(block_count)
        call assertEqual(block_count(lev), block_count_ex(lev), "Wrong # of levels")
    end do

    !!!!! STEP 5/6 - CONFIRM LEVEL 2 ONLY ON LOWER-RIGHT
    ! Single point not in corner cell
    points(:, :) = 0.0d0
    values(:) = 0.0d0
    points(1, :) = [0.9d0, 0.1d0]
    values(1) = REFINE_TO_L2 
    call sim_advance(5, points, values, &
                     "SETTING SINGLE CELL ONLY FOR LEVEL 2", &
                     "LEAVES AFTER LEVEL 2 DATA AT SINGLE CELL", &
                     block_count)
    call gr_getFinestLevel(finest_level)
    call assertEqual(2, finest_level, "Incorrect finest level")

    block_count_ex = [3, 4, 0, 0]
    do lev = 1, SIZE(block_count)
        call assertEqual(block_count(lev), block_count_ex(lev), "Wrong # of levels")
    end do

    !!!!! STEP 7/8 - REFINE TO LEVEL 3 ON POINT
    ! Same point but maximize refinement.  However, refinement 
    ! can only increase one level with each advance.
    values(1) = REFINE_TO_L5
    call sim_advance(7, points, values, &
                     "SETTING SINGLE CELL ONLY FOR LEVEL 4", &
                     "LEAVES AFTER ONLY GETTING TO L3 AT SINGLE CELL", &
                     block_count)
    call gr_getFinestLevel(finest_level)
    call assertEqual(3, finest_level, "Incorrect finest level")

    block_count_ex = [0, 15, 4, 0]
    do lev = 1, SIZE(block_count)
        call assertEqual(block_count(lev), block_count_ex(lev), "Wrong # of levels")
    end do

    ! During this step, gr_remakeLevelCallback called on level 2
    ! and then gr_makeFineLevelFromCoarseCallback created level 3 from level 2
    ! Confirm that data is correct in multifabs at all levels
    do lev = 1, finest_level
        call amrex_mfiter_build(mfi, unk(lev-1), tiling=.FALSE.)
        do while(mfi%next())
            bx = mfi%tilebox()

            ! DEVNOTE: TODO Simulate block until we have a natural iterator for FLASH
            ! Level must be 1-based index and limits/limitsGC must be 1-based also
            ! DEVNOTE: Should we use gr_[ijk]guard here?
            blockDesc%level = lev
            blockDesc%grid_index = -1
            blockDesc%limits(LOW,  :) = 1
            blockDesc%limits(HIGH, :) = 1
            blockDesc%limits(LOW,  1:NDIM) = bx%lo(1:NDIM) + 1
            blockDesc%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1
            blockDesc%limitsGC(LOW,  :) = 1
            blockDesc%limitsGC(HIGH, :) = 1
            blockDesc%limitsGC(LOW,  1:NDIM) = blockDesc%limits(LOW,  1:NDIM) - NGUARD
            blockDesc%limitsGC(HIGH, 1:NDIM) = blockDesc%limits(HIGH, 1:NDIM) + NGUARD

            associate (lo => blockDesc%limits(LOW,  :), &
                       hi => blockDesc%limits(HIGH, :), &
                     loGC => blockDesc%limitsGC(LOW,  :), &
                     hiGC => blockDesc%limitsGC(HIGH, :))
                ! Makes this 1-based cell indexing
                solnData(loGC(1):, loGC(2):, loGC(3):, 1:) => unk(lev-1)%dataptr(mfi)

                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        if ((lev == 1) .AND. (i == 15) .AND. (j == 2)) then
                            ! Level one not affected by remake/makeFine
                            call assertEqual(solnData(i, j, 1, 1), &
                                             values(1) / 64.0d0, &
                                             "Wrong data")
                        else if ((lev == 2) .AND. (i == 29) .AND. (j == 4)) then
                            ! Level two data copied into new UNK multifab with remake
                            call assertEqual(solnData(i, j, 1, 1), &
                                             values(1) / 16.0d0, &
                                             "Wrong data")
                        else if (lev == 3) then
                            ! Level three data from level two
                            if      ((i == 57) .AND. (j == 7)) then
                                call assertEqual(solnData(i, j, 1, 1), &
                                                 values(1) / 16.0d0, &
                                                 "Wrong data")
                            else if ((i == 58) .AND. (j == 7)) then
                                call assertEqual(solnData(i, j, 1, 1), &
                                                 values(1) / 16.0d0, &
                                                 "Wrong data")
                            else if ((i == 57) .AND. (j == 8)) then
                                call assertEqual(solnData(i, j, 1, 1), &
                                                 values(1) / 16.0d0, &
                                                 "Wrong data")
                            else if ((i == 58) .AND. (j == 8)) then
                                call assertEqual(solnData(i, j, 1, 1), &
                                                 values(1) / 16.0d0, &
                                                 "Wrong data")
                            end if
                        else
                            call assertEqual(solnData(i, j, 1, 1), 0.0d0, "Data no zero")
                        end if
                    end do
                end do

                nullify(solnData)
            end associate

        end do
        call amrex_mfiter_destroy(mfi)
    end do

    !!!!! STEP 9/10 - ADVANCE WITH NO CHANGE TO ACHIEVE LEVEL 4
    call sim_advance(9, points, values, &
                     "NO DATA CHANGE - LET IT REFINE TO LEVEL 4", &
                     "LEAVES CONSECUTIVE STEPS TO L4 AT SINGLE CELL", &
                     block_count)
    call gr_getFinestLevel(finest_level)
    call assertEqual(4, finest_level, "Incorrect finest level")

    block_count_ex = [0, 12, 15, 4]
    do lev = 1, SIZE(block_count)
        call assertEqual(block_count(lev), block_count_ex(lev), "Wrong # of levels")
    end do

    !!!!! STEP 11/12 - ADVANCE WITH NO CHANGE AND CONFIRM NO CHANGE
    ! We should be limited to refinement up to level 4
    call sim_advance(11, points, values, &
                     "NO DATA CHANGE -  STUCK AT REFINEMENT LEVEL 4", &
                     "LEAVES CONSECUTIVE STEPS TO L4 AT SINGLE CELL", &
                     block_count)
    call gr_getFinestLevel(finest_level)
    call assertEqual(4, finest_level, "Incorrect finest level")

    block_count_ex = [0, 12, 15, 4]
    do lev = 1, SIZE(block_count)
        call assertEqual(block_count(lev), block_count_ex(lev), "Wrong # of levels")
    end do

    !!!!! STEP 13-16 - ADD ONE MORE LEVEL 4 POINT
    points(:, :) = 0.0d0
    values(:) = 0.0d0
    points(1, :) = [0.9d0,  0.1d0]
    points(2, :) = [0.29d0, 0.58d0]
    values(1) = REFINE_TO_L4 
    values(2) = REFINE_TO_L4
    call sim_advance(13, points, values, &
                     "SETTING SECOND LEVEL 4 CELL", &
                     "LEAVES AFTER SECOND LEVEL 4 DATA", &
                     block_count)
    call gr_getFinestLevel(finest_level)
    call assertEqual(4, finest_level, "Incorrect finest level")

    call sim_advance(15, points, values, &
                     "SETTING SECOND LEVEL 4 CELL", &
                     "LEAVES AFTER SECOND LEVEL 4 DATA", &
                     block_count)
    call gr_getFinestLevel(finest_level)
    call assertEqual(4, finest_level, "Incorrect finest level")

    block_count_ex = [0, 8, 30, 8]
    do lev = 1, SIZE(block_count)
        call assertEqual(block_count(lev), block_count_ex(lev), "Wrong # of levels")
    end do

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


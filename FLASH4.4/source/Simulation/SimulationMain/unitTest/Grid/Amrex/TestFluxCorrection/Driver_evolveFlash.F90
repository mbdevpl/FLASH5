!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestFluxCorrection/Driver_evolveFlash
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
!!  The driver for a toy version of a full FLASH simulation that tests if the
!!  AMReX implementation works well when coupled with the FLASH's Grid unit
!!  interface for controlling flux correction.
!!  
!!  This test refines the mesh at initialization so that level 1 has no leaf
!!  blocks and level three covers only a few level-2 blocks.  Several level 3
!!  blocks are placed with an edge at the domain boundary to test correct 
!!  detection of fine-coarse boundaries with and without periodic BC.
!!
!! NOTES
!!  This simulation *must* be configured with at least the following
!!  2D run:
!!     ./setup -auto -2d -nxb=8 -nyb=8
!!              unitTest/Grid/Amrex/TestFluxCorrection
!!             +noio -index-reorder
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Driver_evolveFlash()
    use amrex_fort_module,     ONLY : amrex_spacedim
    use amrex_box_module,      ONLY : amrex_box
    use amrex_multifab_module, ONLY : amrex_mfiter, &
                                      amrex_mfiter_build, &
                                      amrex_mfiter_destroy

    use Grid_interface,        ONLY : Grid_getFluxPtr, Grid_releaseFluxPtr, &
                                      Grid_zeroFluxRegister, &
                                      Grid_addFineToFluxRegister, &
                                      Grid_conserveFluxes
    use gr_amrexInterface,     ONLY : gr_getFinestLevel
    use gr_iterator,           ONLY : gr_iterator_t, &
                                      build_iterator, destroy_iterator
    use block_metadata,        ONLY : block_metadata_t
    use ut_testDriverMod

    implicit none

    type(gr_iterator_t)    :: itor
    type(block_metadata_t) :: block
    real, pointer          :: fluxDataX(:, :, :, :) => null() 
    real, pointer          :: fluxDataY(:, :, :, :) => null()
    real, pointer          :: fluxDataZ(:, :, :, :) => null()

    integer :: finest_level
    integer :: lev, i, j

    !!!!! CONFIRM PROPER DIMENSIONALITY
    write(*,*)
    if (amrex_spacedim /= 2) then
        write(*,*) "Wrong dimensionality - ", amrex_spacedim, ' != ', 2
        write(*,*) "Recompile AMReX with correct dimensionality"
        write(*,*)
        stop
    end if

    call start_test_run

    !!!!! CONFIRM PROPER SETUP
    ! 16x16 domain with dx = dy = (1.0 - 0.0)/16
    call assertEqual(NXB, 8, "Incorrect initial number of cells/block along X")
    call assertEqual(NYB, 8, "Incorrect initial number of cells/block along Y")
    call assertEqual(NZB, 1, "Incorrect initial number of cells/block along Z")

    call gr_getFinestLevel(finest_level)
    call assertEqual(3, finest_level, "Incorrect finest level")

    !!!!! WRITE FLUX DATA TO ALL BLOCKS AT ALL POPULATED LEVELS
    ! Write constant flux to each level that has blocks
    do lev = 1, finest_level
        call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)

        do while (itor%is_valid())
            call itor%blkMetaData(block)

            call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
            
            call assertTrue( associated(fluxDataX), "flux X Ptr should be set")
            call assertTrue( associated(fluxDataY), "flux Y Ptr should be set")
            call assertFalse(associated(fluxDataZ), "flux Z Ptr should be NULL")

            call assertEqual(SIZE(fluxDataX, 1), 9, "Incorrect width of X block")
            call assertEqual(SIZE(fluxDataX, 2), 8, "Incorrect height of X block")
            call assertEqual(SIZE(fluxDataY, 1), 8, "Incorrect width of Y block")
            call assertEqual(SIZE(fluxDataY, 2), 9, "Incorrect height of Y block")

            fluxDataX(:, :, :, :) =  lev
            fluxDataY(:, :, :, :) = -lev

            call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

            call itor%next()
        end do

        call destroy_iterator(itor)
    end do

    ! Confirm proper writing of flux data 
    do lev = 1, finest_level
        call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)

        do while (itor%is_valid())
            call itor%blkMetaData(block)

            call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
            
            associate(lo => block%limits(LOW,  :), &
                      hi => block%limits(HIGH, :))
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        call assertEqual(DBLE( lev), fluxDataX(i, j, 1, 1), &
                                         "Incorrect X flux data on level")
                        call assertEqual(DBLE(-lev), fluxDataY(i, j, 1, 1), &
                                         "Incorrect Y flux data on level")
                    end do
                end do
            end associate

            call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
            
            call itor%next()
        end do
        
        call destroy_iterator(itor)
    end do

    ! Run operations to be tested
    do lev = finest_level, 1, -1
        ! Conserve flux at fine-coarse boundaries
        if (lev < finest_level) then
            call Grid_conserveFluxes(ALLDIR, lev)
        end if

        ! Save corrected flux for next level
        if (lev > 1) then
            call Grid_zeroFluxRegister(lev)
            call Grid_addFineToFluxRegister(lev)
        end if
    end do

    ! No change to flux on level 3 (all leaf blocks)
    lev = 3
    call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)
    do while (itor%is_valid())
        call itor%blkMetaData(block)

        call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        associate(lo => block%limits(LOW,  :), &
                  hi => block%limits(HIGH, :))
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    call assertEqual(DBLE( lev), fluxDataX(i, j, 1, 1), &
                                     "Incorrect X flux data on level 3")
                    call assertEqual(DBLE(-lev), fluxDataY(i, j, 1, 1), &
                                     "Incorrect Y flux data on level 3")
                end do
            end do
        end associate

        call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
        
        call itor%next()
    end do
    call destroy_iterator(itor)
 
    ! Change to flux on level 2 (fine-coarse boundaries)
    lev = 2
    call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)
    do while (itor%is_valid())
        call itor%blkMetaData(block)

        call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        associate(lo => block%limits(LOW,  :), &
                  hi => block%limits(HIGH, :))
            ! Check X-face fluxes
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)+1
                    ! Parent Block (25,1) to (32,8)
                    ! Left face is at fine-fine boundary and right face
                    ! is on outflow boundary.  Therefore, no flux overwriting
                    ! on either.
                    if      (     ((i ==  9) .AND. (j ==  9)) &
                             .OR. ((i ==  9) .AND. (j == 10)) &
                             .OR. ((i ==  9) .AND. (j == 11)) &
                             .OR. ((i ==  9) .AND. (j == 12)) &
                             .OR. ((i ==  9) .AND. (j == 13)) &
                             .OR. ((i ==  9) .AND. (j == 14)) &
                             .OR. ((i ==  9) .AND. (j == 15)) &
                             .OR. ((i ==  9) .AND. (j == 16)) &
                             .OR. ((i == 17) .AND. (j ==  9)) &
                             .OR. ((i == 17) .AND. (j == 10)) &
                             .OR. ((i == 17) .AND. (j == 11)) &
                             .OR. ((i == 17) .AND. (j == 12)) &
                             .OR. ((i == 17) .AND. (j == 13)) &
                             .OR. ((i == 17) .AND. (j == 14)) &
                             .OR. ((i == 17) .AND. (j == 15)) &
                             .OR. ((i == 17) .AND. (j == 16)) ) then
                        ! Parent Block (9,9) to (16,16)
                        ! No X faces on domain.  Both X faces are at fine-coarse
                        ! boundary
                        call assertEqual(DBLE(lev+1), fluxDataX(i, j, 1, 1), &
                                         "Incorrect X flux data on level 2")
                    else if (     ((i == 17) .AND. (j == 25)) &
                             .OR. ((i == 17) .AND. (j == 26)) &
                             .OR. ((i == 17) .AND. (j == 27)) &
                             .OR. ((i == 17) .AND. (j == 28)) &
                             .OR. ((i == 17) .AND. (j == 29)) &
                             .OR. ((i == 17) .AND. (j == 30)) &
                             .OR. ((i == 17) .AND. (j == 31)) &
                             .OR. ((i == 17) .AND. (j == 32)) &
                             .OR. ((i == 25) .AND. (j == 25)) &
                             .OR. ((i == 25) .AND. (j == 26)) &
                             .OR. ((i == 25) .AND. (j == 27)) &
                             .OR. ((i == 25) .AND. (j == 28)) &
                             .OR. ((i == 25) .AND. (j == 29)) &
                             .OR. ((i == 25) .AND. (j == 30)) &
                             .OR. ((i == 25) .AND. (j == 31)) &
                             .OR. ((i == 25) .AND. (j == 32)) ) then
                        ! Parent Block (17,25) to (24,32)
                        ! No X faces on domain boundary.  Both X faces are at
                        ! fine-coarse boundary
                        call assertEqual(DBLE(lev+1), fluxDataX(i, j, 1, 1), &
                                         "Incorrect X flux data on level 2")
                    else if (     ((i == 17) .AND. (j == 1)) &
                             .OR. ((i == 17) .AND. (j == 2)) &
                             .OR. ((i == 17) .AND. (j == 3)) &
                             .OR. ((i == 17) .AND. (j == 4)) &
                             .OR. ((i == 17) .AND. (j == 5)) &
                             .OR. ((i == 17) .AND. (j == 6)) &
                             .OR. ((i == 17) .AND. (j == 7)) &
                             .OR. ((i == 17) .AND. (j == 8)) ) then
                        ! Parent Block (17,1) to (24,8)
                        ! No X faces on domain boundary.  Left face is at F-C
                        ! boundary.  Right face is not.
                        call assertEqual(DBLE(lev+1), fluxDataX(i, j, 1, 1), &
                                         "Incorrect X flux data on level 2")
                    else 
                        call assertEqual(DBLE(lev  ), fluxDataX(i, j, 1, 1), &
                                         "Incorrect X flux data on level 2")
                    end if
                end do
            end do

            ! Check Y-face fluxes
            do     j = lo(JAXIS), hi(JAXIS)+1
                do i = lo(IAXIS), hi(IAXIS)
                    if      (     ((i ==  9) .AND. (j ==  9)) &
                             .OR. ((i == 10) .AND. (j ==  9)) &
                             .OR. ((i == 11) .AND. (j ==  9)) &
                             .OR. ((i == 12) .AND. (j ==  9)) &
                             .OR. ((i == 13) .AND. (j ==  9)) &
                             .OR. ((i == 14) .AND. (j ==  9)) &
                             .OR. ((i == 15) .AND. (j ==  9)) &
                             .OR. ((i == 16) .AND. (j ==  9)) &
                             .OR. ((i ==  9) .AND. (j == 17)) &
                             .OR. ((i == 10) .AND. (j == 17)) &
                             .OR. ((i == 11) .AND. (j == 17)) &
                             .OR. ((i == 12) .AND. (j == 17)) &
                             .OR. ((i == 13) .AND. (j == 17)) &
                             .OR. ((i == 14) .AND. (j == 17)) &
                             .OR. ((i == 15) .AND. (j == 17)) &
                             .OR. ((i == 16) .AND. (j == 17)) ) then
                        ! Parent Block (9,9) to (16,16)
                        ! No Y faces on domain.  Both Y faces are at fine-coarse
                        ! boundary 
                        call assertEqual(DBLE(-lev-1), fluxDataY(i, j, 1, 1), &
                                         "Incorrect Y flux data on level 2")
                    else if (     ((i == 17) .AND. (j == 25)) &
                             .OR. ((i == 18) .AND. (j == 25)) &
                             .OR. ((i == 19) .AND. (j == 25)) &
                             .OR. ((i == 20) .AND. (j == 25)) &
                             .OR. ((i == 21) .AND. (j == 25)) &
                             .OR. ((i == 22) .AND. (j == 25)) &
                             .OR. ((i == 23) .AND. (j == 25)) &
                             .OR. ((i == 24) .AND. (j == 25)) ) then
                        ! Parent Block (17,25) to (24,32)
                        ! Top Y-face is at periodic boundary with level 3 block
                        ! on other side.  Bottom face is at F-C boundary
                        call assertEqual(DBLE(-lev-1), fluxDataY(i, j, 1, 1), &
                                         "Incorrect Y flux data on level 2")
                    else if (     ((i == 17) .AND. (j == 9)) &
                             .OR. ((i == 18) .AND. (j == 9)) &
                             .OR. ((i == 19) .AND. (j == 9)) &
                             .OR. ((i == 20) .AND. (j == 9)) &
                             .OR. ((i == 21) .AND. (j == 9)) &
                             .OR. ((i == 22) .AND. (j == 9)) &
                             .OR. ((i == 23) .AND. (j == 9)) &
                             .OR. ((i == 24) .AND. (j == 9)) ) then
                        ! Parent Block (17,1) to (24,8)
                        ! Bottom Y-face is at periodic boundary with level 3 block
                        ! on other side.  Top face is at F-C boundary
                        call assertEqual(DBLE(-lev-1), fluxDataY(i, j, 1, 1), &
                                         "Incorrect Y flux data on level 2")
                    else if (     ((i == 25) .AND. (j == 1)) &
                             .OR. ((i == 26) .AND. (j == 1)) &
                             .OR. ((i == 27) .AND. (j == 1)) &
                             .OR. ((i == 28) .AND. (j == 1)) &
                             .OR. ((i == 29) .AND. (j == 1)) &
                             .OR. ((i == 30) .AND. (j == 1)) &
                             .OR. ((i == 31) .AND. (j == 1)) &
                             .OR. ((i == 32) .AND. (j == 1)) &
                             .OR. ((i == 25) .AND. (j == 9)) &
                             .OR. ((i == 26) .AND. (j == 9)) &
                             .OR. ((i == 27) .AND. (j == 9)) &
                             .OR. ((i == 28) .AND. (j == 9)) &
                             .OR. ((i == 29) .AND. (j == 9)) &
                             .OR. ((i == 30) .AND. (j == 9)) &
                             .OR. ((i == 31) .AND. (j == 9)) &
                             .OR. ((i == 32) .AND. (j == 9)) ) then
                        ! Parent Block (25,1) to (32,8)
                        ! Bottom Y-face is at periodic boundary with level 2 block
                        ! on other side.  Top face is at F-C boundary
                        call assertEqual(DBLE(-lev-1), fluxDataY(i, j, 1, 1), &
                                         "Incorrect Y flux data on level 2")
                    else if (     ((i == 25) .AND. (j == 33)) &
                             .OR. ((i == 26) .AND. (j == 33)) &
                             .OR. ((i == 27) .AND. (j == 33)) &
                             .OR. ((i == 28) .AND. (j == 33)) &
                             .OR. ((i == 29) .AND. (j == 33)) &
                             .OR. ((i == 30) .AND. (j == 33)) &
                             .OR. ((i == 31) .AND. (j == 33)) &
                             .OR. ((i == 32) .AND. (j == 33)) ) then
                        ! Leaf Block (25,25) to (32,32)
                        ! Top face on domain boundary with periodic BC and with
                        ! level 3 block on other side
                        call assertEqual(DBLE(-lev-1), fluxDataY(i, j, 1, 1), &
                                         "Incorrect Y flux data on level 2")
                    else 
                        call assertEqual(DBLE(-lev  ), fluxDataY(i, j, 1, 1), &
                                         "Incorrect Y flux data on level 2")
                    end if
                end do
            end do
        end associate

        call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        call itor%next()
    end do
    call destroy_iterator(itor)
    
    ! No change to flux on level 1 (no leaf blocks)
    lev = 1
    call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)
    do while (itor%is_valid())
        call itor%blkMetaData(block)

        call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        associate(lo => block%limits(LOW,  :), &
                  hi => block%limits(HIGH, :))
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    call assertEqual(DBLE( lev), fluxDataX(i, j, 1, 1), &
                                     "Incorrect X flux data on level 1")
                    call assertEqual(DBLE(-lev), fluxDataY(i, j, 1, 1), &
                                     "Incorrect Y flux data on level 1")
                end do
            end do
        end associate

        call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        call itor%next()
    end do
    call destroy_iterator(itor)

    call finish_test_run

end subroutine Driver_evolveFlash

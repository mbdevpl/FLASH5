!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestFluxCorrection2/Driver_evolveFlash
!!
!! NAME
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!  The driver for a toy version of a full FLASH simulation that tests if the
!!  AMReX implementation works well when coupled with the FLASH's Grid unit
!!  interface for controlling flux correction.  In particular, it is testing if
!!  the ability to compute flux errors in the flux registers is correctly 
!!  implemented.
!!  
!!  This test refines the mesh at initialization so that level 1 has no leaf
!!  blocks and level three covers only one level-2 block.
!!
!! NOTES
!!  This simulation *must* be configured with at least the following
!!  2D run:
!!     ./setup -auto -2d -nxb=8 -nyb=8
!!              unitTest/Grid/Amrex/TestFluxCorrection2
!!             +noio -index-reorder
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Driver_evolveFlash()
    use amrex_fort_module,     ONLY : amrex_spacedim

    use Grid_interface,        ONLY : Grid_getFluxPtr, Grid_releaseFluxPtr, &
                                      Grid_zeroFluxRegister, &
                                      Grid_addFineToFluxRegister, &
                                      Grid_addCoarseToFluxRegister, &
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
    write(*,*) "Initialization Test Phase"
    write(*,*) "------------------------------------------------------------"
    call assertEqual(NXB, 8, "Incorrect initial number of cells/block along X")
    call assertEqual(NYB, 8, "Incorrect initial number of cells/block along Y")
    call assertEqual(NZB, 1, "Incorrect initial number of cells/block along Z")

    call gr_getFinestLevel(finest_level)
    call assertEqual(3, finest_level, "Incorrect finest level")

    !!!!! ZERO OUT ALL DATA & FLUX REGISTERS AND CONFIRM ZERO DATA IN REGISTERS
    call Grid_zeroFluxData

    ! Zero all flux register data and transfer to flux data structure
    do lev = 1, finest_level
        if (lev > 1) then
            call Grid_zeroFluxRegister(lev)
        end if

        if (lev > finest_level) then
            call Grid_conserveFluxes(ALLDIR, lev)
        end if
    end do

    ! Zero flux at all faces
    call build_iterator(itor, ALL_BLKS, tiling=.FALSE.)
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

        associate(lo => block%limits(LOW,  :), &
                  hi => block%limits(HIGH, :))
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    call assertEqual(DBLE(0.0), fluxDataX(i, j, 1, 1), &
                                     "Incorrect X flux data on level")
                    call assertEqual(DBLE(0.0), fluxDataY(i, j, 1, 1), &
                                     "Incorrect Y flux data on level")
                end do
            end do
        end associate
        call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        call itor%next()
    end do
    call destroy_iterator(itor)

    !!!!! ADD ONLY COARSE FLUX TO REGISTERS
    write(*,*) "Coarse Flux ONLY Test Phase"
    write(*,*) "------------------------------------------------------------"
    do lev = 2, finest_level
        call Grid_zeroFluxRegister(lev)
    end do

    do lev = 1, finest_level
        call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)
        do while (itor%is_valid())
            call itor%blkMetaData(block)

            call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
            associate(lo => block%limits(LOW,  :), &
                      hi => block%limits(HIGH, :))
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)+1
                        fluxDataX(i, j, 1, 1) = j
                    end do
                end do
                do     j = lo(JAXIS), hi(JAXIS)+1
                    do i = lo(IAXIS), hi(IAXIS)
                        fluxDataY(i, j, 1, 1) = -i
                    end do
                end do
            end associate
            call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

            call itor%next()
        end do
        call destroy_iterator(itor)

        if (lev < finest_level) then
            call Grid_addCoarseToFluxRegister(lev)
        end if
    end do

    call Grid_zeroFluxData
    do lev = 1, finest_level-1
        call Grid_conserveFluxes(ALLDIR, lev)
    end do

    ! No Fine/Coarse boundaries on coarsest level
    lev = 1
    call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)
    do while (itor%is_valid())
        call itor%blkMetaData(block)

        call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
        associate(lo => block%limits(LOW,  :), &
                  hi => block%limits(HIGH, :))
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    call assertEqual(DBLE(0.0), fluxDataX(i, j, 1, 1), &
                                     "Incorrect X flux data on level 1")
                    call assertEqual(DBLE(0.0), fluxDataY(i, j, 1, 1), &
                                     "Incorrect Y flux data on level 1")
                end do
            end do
        end associate
        call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        call itor%next()
    end do
    call destroy_iterator(itor)

    ! Second layer has one leaf block on top
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
                    if (     ((i ==  9) .AND. (j ==  9)) &
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
                        call assertEqual(DBLE(j), fluxDataX(i, j, 1, 1), &
                                         "Incorrect X flux data on level 2")
                    else
                        call assertEqual(DBLE(0.0), fluxDataX(i, j, 1, 1), &
                                         "Incorrect X flux data on level 2")
                    end if
                end do
            end do

            ! Check Y-face fluxes
            do     j = lo(JAXIS), hi(JAXIS)+1
                do i = lo(IAXIS), hi(IAXIS)
                    if (     ((i ==  9) .AND. (j ==  9)) &
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
                        call assertEqual(DBLE(-i), fluxDataY(i, j, 1, 1), &
                                         "Incorrect Y flux data on level 2")
                    else
                        call assertEqual(DBLE(0.0), fluxDataY(i, j, 1, 1), &
                                         "Incorrect Y flux data on level 2")
                    end if
                end do
            end do
        end associate
        call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        call itor%next()
    end do
    call destroy_iterator(itor)

    ! No Fine/Coarse boundaries on finest level
    lev = 3
    call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)
    do while (itor%is_valid())
        call itor%blkMetaData(block)
        call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
        associate(lo => block%limits(LOW,  :), &
                  hi => block%limits(HIGH, :))
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    call assertEqual(DBLE(0.0), fluxDataX(i, j, 1, 1), &
                                     "Incorrect X flux data on level 3")
                    call assertEqual(DBLE(0.0), fluxDataY(i, j, 1, 1), &
                                     "Incorrect Y flux data on level 3")
                end do
            end do
        end associate
        call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        call itor%next()
    end do
    call destroy_iterator(itor)

    !!!!! CREATE FLUX ERROR IN REGISTERS
    write(*,*) "Flux Error Test Phase"
    write(*,*) "------------------------------------------------------------"
    do lev = 2, finest_level
        call Grid_zeroFluxRegister(lev)
    end do

    do lev = finest_level, 1, -1
        call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)
        do while (itor%is_valid())
            call itor%blkMetaData(block)

            call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
            associate(lo => block%limits(LOW,  :), &
                      hi => block%limits(HIGH, :))
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)+1
                        fluxDataX(i, j, 1, 1) = j
                    end do
                end do
                do     j = lo(JAXIS), hi(JAXIS)+1
                    do i = lo(IAXIS), hi(IAXIS)
                        fluxDataY(i, j, 1, 1) = -i
                    end do
                end do
            end associate
            call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

            call itor%next()
        end do
        call destroy_iterator(itor)

        ! Construct flux error in flux register whose coarse flux
        ! is associated with current level
        if (lev < finest_level) then
            call Grid_addFineToFluxRegister(  lev+1, coefficient= 1.0)
            call Grid_addCoarseToFluxRegister(lev  , coefficient=-1.0)
        end if
    end do

    call Grid_zeroFluxData
    do lev = 1, finest_level-1
        call Grid_conserveFluxes(ALLDIR, lev)
    end do

    ! No Fine/Coarse boundaries on coarsest level
    lev = 1
    call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)
    do while (itor%is_valid())
        call itor%blkMetaData(block)
        call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
        associate(lo => block%limits(LOW,  :), &
                  hi => block%limits(HIGH, :))
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    call assertEqual(DBLE(0.0), fluxDataX(i, j, 1, 1), &
                                     "Incorrect X flux data on level 1")
                    call assertEqual(DBLE(0.0), fluxDataY(i, j, 1, 1), &
                                     "Incorrect Y flux data on level 1")
                end do
            end do
        end associate
        call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        call itor%next()
    end do
    call destroy_iterator(itor)

    ! Second layer had one block refined 
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
                    if (     ((i ==  9) .AND. (j ==  9)) &
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
                        ! Coarse side has length A => Fine side has A/2
                        ! For coarse cell with y-index j, the flux is 
                        !    jA
                        ! The total flux in the two adjacent fine cells is
                        !    2*j*A/2 + (2*j - 1)*A/2
                        !
                        ! Therefore, the flux difference is (2*j - 1)*A/2 and
                        ! the flux density is (2*j - 1)/2.
                        call assertEqual(DBLE(0.5*(2.0*j - 1)), fluxDataX(i, j, 1, 1), &
                                         "Incorrect X flux data on level 2")
                    else
                        call assertEqual(DBLE(0.0), fluxDataX(i, j, 1, 1), &
                                         "Incorrect X flux data on level 2")
                    end if
                end do
            end do

            ! Check Y-face fluxes
            do     j = lo(JAXIS), hi(JAXIS)+1
                do i = lo(IAXIS), hi(IAXIS)
                    if (     ((i ==  9) .AND. (j ==  9)) &
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
                        ! See X-face comments
                        call assertEqual(DBLE(-0.5*(2.0*i - 1)), fluxDataY(i, j, 1, 1), &
                                         "Incorrect Y flux data on level 2")
                    else
                        call assertEqual(DBLE(0.0), fluxDataY(i, j, 1, 1), &
                                         "Incorrect Y flux data on level 2")
                    end if
                end do
            end do
        end associate
        call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        call itor%next()
    end do
    call destroy_iterator(itor)

    ! No Fine/Coarse boundaries on finest level
    lev = 3
    call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)
    do while (itor%is_valid())
        call itor%blkMetaData(block)
        call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
        associate(lo => block%limits(LOW,  :), &
                  hi => block%limits(HIGH, :))
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    call assertEqual(DBLE(0.0), fluxDataX(i, j, 1, 1), &
                                     "Incorrect X flux data on level 1")
                    call assertEqual(DBLE(0.0), fluxDataY(i, j, 1, 1), &
                                     "Incorrect Y flux data on level 1")
                end do
            end do
        end associate
        call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        call itor%next()
    end do
    call destroy_iterator(itor)

    !!!!! USE FINE/COARSE MULTIPLICATIVE FACTORS TO CANCEL FLUX
    write(*,*) "Flux Cancellation Test Phase"
    write(*,*) "------------------------------------------------------------"
    do lev = 2, finest_level
        call Grid_zeroFluxRegister(lev)
    end do

    do lev = finest_level, 1, -1
        call build_iterator(itor, ALL_BLKS, lev, tiling=.FALSE.)
        do while (itor%is_valid())
            call itor%blkMetaData(block)

            call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
            fluxDataX(:, :, :, :) =  lev
            fluxDataY(:, :, :, :) = -lev
            call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

            call itor%next()
        end do
        call destroy_iterator(itor)

        ! Construct flux error in flux register whose coarse flux
        ! is associated with current level
        if (lev < finest_level) then
            call Grid_addFineToFluxRegister(  lev+1, coefficient=-DBLE(lev)/DBLE(lev+1))
            call Grid_addCoarseToFluxRegister(lev  , coefficient=1.0)
        end if
    end do

    call Grid_zeroFluxData
    do lev = 1, finest_level-1
        call Grid_conserveFluxes(ALLDIR, lev)
    end do

    ! Zero flux everywhere
    call build_iterator(itor, ALL_BLKS, tiling=.FALSE.)
    do while (itor%is_valid())
        call itor%blkMetaData(block)
        call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
        associate(lo => block%limits(LOW,  :), &
                  hi => block%limits(HIGH, :))
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    call assertEqual(DBLE(0.0), fluxDataX(i, j, 1, 1), &
                                     "Incorrect X flux data on level")
                    call assertEqual(DBLE(0.0), fluxDataY(i, j, 1, 1), &
                                     "Incorrect Y flux data on level")
                end do
            end do
        end associate
        call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

        call itor%next()
    end do
    call destroy_iterator(itor)

    call finish_test_run

end subroutine Driver_evolveFlash


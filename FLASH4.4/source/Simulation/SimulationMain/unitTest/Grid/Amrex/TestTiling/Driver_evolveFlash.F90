!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestTiling/Driver_evolveFlash
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
!!  This unit test exercises the flash_iterator and flash_tile objects both with
!!  and without tiling enabled to confirm that iterating over tiles in accord
!!  with tiling runtime parameters is correct.  In addition, it confirms that
!!  tiling related features of the flash_tile_t objects are correct.
!!  
!! NOTES
!!  This simulation *must* be configured with at least the following
!!  2D run:
!!     ./setup -auto -2d -nxb=8 -nyb=8 
!!              unitTest/Grid/Amrex/TestTiling
!!              +noio -index-reorder
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Driver_evolveFlash()
    use amrex_box_module,      ONLY : amrex_box

    use flash_iterator,        ONLY : flash_iterator_t
    use flash_tile,            ONLY : flash_tile_t
    use Grid_interface,        ONLY : Grid_getTileIterator, &
                                      Grid_releaseTileIterator
    use Grid_data,             ONLY : gr_enableTiling, &
                                      gr_tileSize
    use ut_testDriverMod

    implicit none

    integer :: cnt

    type(flash_iterator_t) :: itor
    type(flash_tile_t)     :: tileDesc
    type(flash_tile_t)     :: encBlk

    type(amrex_box) :: tileBox
    type(amrex_box) :: encBox

    real, pointer :: solnData(:,:,:,:)
    integer       :: bnd(1:4)

    integer :: i, j, k, var
    integer :: axis

    nullify(solnData)

    call start_test_run

    ! Confirm contents of parfile that are hardcoded into the expected results
    ! and behavior of the tests
    call assertTrue(gr_enableTiling, "Tiling not enabled")
    call assertEqual(gr_tileSize(IAXIS), 4, "Incorrect X tile size")
    call assertEqual(gr_tileSize(JAXIS), 2, "Incorrect Y tile size")
    call assertEqual(gr_tileSize(KAXIS), 1, "Incorrect Z tile size")

    call assertEqual(NXB, 8, "Invalid number of cells/block in X")
    call assertEqual(NYB, 8, "Invalid number of cells/block in Y")
    call assertEqual(NZB, 1, "Invalid number of cells/block in Z")

    call assertEqual(NGUARD, 4, "Invalid number of guardcells")

    ! Get total blocks by explicitly turning off tiling at iterator creation
    cnt = 0
    call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.)
    do while(itor%isValid())
        call itor%currentTile(tileDesc)

        associate(lo      => tileDesc%limits(LOW,  :), &
                  hi      => tileDesc%limits(HIGH, :), &
                  loGC    => tileDesc%grownLimits(LOW,  :), &
                  hiGC    => tileDesc%grownLimits(HIGH, :), &
                  blkloGC => tileDesc%blkLimitsGC(LOW,  :), &
                  blkhiGC => tileDesc%blkLimitsGC(HIGH, :))
            ! Confirm that the tile we are given has the size of a block
            ! (the trivial tiling)
            call assertEqual(NXB, hi(IAXIS) - lo(IAXIS) + 1, "Invalid tile x-length")
            call assertEqual(NYB, hi(JAXIS) - lo(JAXIS) + 1, "Invalid tile y-length")
            call assertEqual(NZB, hi(KAXIS) - lo(KAXIS) + 1, "Invalid tile z-length")

            ! Confirm that grownLimits has the size of the block interior + GC halo 
            call assertEqual(NXB+K1D*2*NGUARD, hiGC(IAXIS) - loGC(IAXIS) + 1, &
                             "Invalid tileGC x-length")
            call assertEqual(NYB+K2D*2*NGUARD, hiGC(JAXIS) - loGC(JAXIS) + 1, &
                             "Invalid tileGC y-length")
            call assertEqual(NZB+K3D*2*NGUARD, hiGC(KAXIS) - loGC(KAXIS) + 1, &
                             "Invalid tileGC z-length")

            ! When the tile is a block, then the grown tile is the
            ! interior + GC halo
            call assertTrue(ALL(loGC == blkLoGC), &
                            "blkLimitsGC low != grownLimits low for block")
            call assertTrue(ALL(hiGC == blkHiGC), &
                            "blkLimitsGC high != grownLimits high for block")
        end associate

        ! Confirm that when the tile is a block that the enclosing block is the
        ! tile itself
        encBlk = tileDesc%enclosingBlock()
        call assertEqual(encBlk%grid_index, tileDesc%grid_index, &
                         "Incorrect enclosing block grid index")
        call assertEqual(encBlk%tile_index, 0, &
                         "Incorrect tile index")
        call assertEqual(encBlk%tile_index, tileDesc%tile_index, &
                         "Incorrect enclosing block tile index")
 
        associate(lo      => encBlk%limits(LOW,  :), &
                  hi      => encBlk%limits(HIGH, :), &
                  loGC    => encBlk%grownLimits(LOW,  :), &
                  hiGC    => encBlk%grownLimits(HIGH, :), &
                  blkloGC => encBlk%blkLimitsGC(LOW,  :), &
                  blkhiGC => encBlk%blkLimitsGC(HIGH, :))
            call assertTrue(ALL(lo == tileDesc%limits(LOW,  :)), &
                            "Bad enclosed block limits LOW")
            call assertTrue(ALL(hi == tileDesc%limits(HIGH, :)), &
                            "Bad enclosed block limits HIGH")
            call assertTrue(ALL(loGC == tileDesc%grownLimits(LOW,  :)), &
                            "Bad enclosed block grownLimits LOW")
            call assertTrue(ALL(hiGC == tileDesc%grownLimits(HIGH, :)), &
                            "Bad enclosed block grownLimits HIGH")
            call assertTrue(ALL(blkloGC == tileDesc%blkLimitsGC(LOW,  :)), &
                            "Bad enclosed block blkLimitsGC LOW")
            call assertTrue(ALL(blkhiGC == tileDesc%blkLimitsGC(HIGH, :)), &
                            "Bad enclosed block blkLimitsGC HIGH")
        end associate

        ! Confirm that the UNK data pointer is for the FAB of the enclosing
        ! block and including the block's GC
        call tileDesc%getDataPtr(solnData, CENTER)
        call assertTrue(associated(solnData), "Unable to get cell-centered data")
        bnd = lbound(solnData)
        call assertTrue( ALL(bnd(1:MDIM) == tileDesc%blkLimitsGC(LOW,  1:MDIM)), &
                         "lbound of data pointer is wrong size")
        bnd = ubound(solnData)
        call assertTrue( ALL(bnd(1:MDIM) == tileDesc%blkLimitsGC(HIGH, 1:MDIM)), &
                         "ubound of data pointer is wrong size")
        call tileDesc%releaseDataPtr(solnData, CENTER)

        ! The flux multifabs do not contain GC data.  Confirm that the FAB
        ! has the size of the enclosing block rather than the size of the
        ! tile
        call tileDesc%getDataPtr(solnData, FLUXX)
        call assertTrue(associated(solnData), "Unable to get FACEX data")
        bnd = ubound(solnData) - lbound(solnData) + 1
        call assertEqual(NXB+1, bnd(IAXIS), "Wrong flux FAB X size")
        call assertEqual(NYB,   bnd(JAXIS), "Wrong flux FAB Y size")
        call assertEqual(NZB,   bnd(KAXIS), "Wrong flux FAB Z size")
        call tileDesc%releaseDataPtr(solnData, FLUXX)

        cnt = cnt + 1

        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)
    call assertEqual(cnt, 2, "Incorrect number of blocks")

    ! *YES* TILING ITERATION
    ! -----------------------------------------------------------------------
    cnt = 0
    call Grid_getTileIterator(itor, ALL_BLKS, tiling=.TRUE.)
    do while(itor%isValid())
        call itor%currentTile(tileDesc)

        associate(lo      => tileDesc%limits(LOW,  :), &
                  hi      => tileDesc%limits(HIGH, :), &
                  loGC    => tileDesc%grownLimits(LOW,  :), &
                  hiGC    => tileDesc%grownLimits(HIGH, :), &
                  blkloGC => tileDesc%blkLimitsGC(LOW,  :), &
                  blkhiGC => tileDesc%blkLimitsGC(HIGH, :))
            ! Confirm appropriate size of tile (the interior)
            do axis = 1, MDIM
                call assertEqual(gr_tileSize(axis), hi(axis) - lo(axis) + 1, &
                                 "Invalid tile length")
            end do 

            ! Confirm that the grown tile has been grown appropriately
            ! Two blocks in x direction
            if      (lo(IAXIS) == 1) then
                call assertEqual(-3, loGC(IAXIS), "Invalid grown tile lo x")
                call assertEqual( 4, hiGC(IAXIS), "Invalid grown tile hi x")
            else if (lo(IAXIS) == 5) then
                call assertEqual( 5, loGC(IAXIS), "Invalid grown tile lo x")
                call assertEqual(12, hiGC(IAXIS), "Invalid grown tile hi x")
            else if (lo(IAXIS) == 9) then
                call assertEqual( 5, loGC(IAXIS), "Invalid grown tile lo x")
                call assertEqual(12, hiGC(IAXIS), "Invalid grown tile hi x")
            else if (lo(IAXIS) == 13) then
                call assertEqual(13, loGC(IAXIS), "Invalid grown tile lo x")
                call assertEqual(20, hiGC(IAXIS), "Invalid grown tile hi x")
            end if

            if      (lo(JAXIS) == 1) then
                call assertEqual(-3, loGC(JAXIS), "Invalid grown tile lo y")
                call assertEqual( 2, hiGC(JAXIS), "Invalid grown tile hi y")
            else if (lo(JAXIS) == 3) then
                call assertEqual( 3, loGC(JAXIS), "Invalid grown tile lo y")
                call assertEqual( 4, hiGC(JAXIS), "Invalid grown tile hi y")
            else if (lo(JAXIS) == 5) then
                call assertEqual( 5, loGC(JAXIS), "Invalid grown tile lo y")
                call assertEqual( 6, hiGC(JAXIS), "Invalid grown tile hi y")
            else if (lo(JAXIS) == 7) then
                call assertEqual( 7, loGC(JAXIS), "Invalid grown tile lo y")
                call assertEqual(12, hiGC(JAXIS), "Invalid grown tile hi y")
            end if

            call assertEqual(1, loGC(KAXIS), "Invalid grown tile lo z")
            call assertEqual(1, hiGC(KAXIS), "Invalid grown tile hi z")

            ! The size of blkLimitsGC should be that for a block
            call assertEqual(NXB+K1D*2*NGUARD, blkhiGC(IAXIS) - blkloGC(IAXIS) + 1, "Invalid tile x-length")
            call assertEqual(NYB+K2D*2*NGUARD, blkhiGC(JAXIS) - blkloGC(JAXIS) + 1, "Invalid tile y-length")
            call assertEqual(NZB+K3D*2*NGUARD, blkhiGC(KAXIS) - blkloGC(KAXIS) + 1, "Invalid tile z-length")
        end associate

        encBlk = tileDesc%enclosingBlock()

        tileBox = amrex_box(tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :))
        encBox  = amrex_box(encBlk%limits(LOW, :), encBlk%limits(HIGH, :))
        call assertTrue(encBox%contains(tileBox), &
                        "tile not contained in enclosing block")

        associate(lo      => encBlk%limits(LOW,  :), &
                  hi      => encBlk%limits(HIGH, :), &
                  loGC    => encBlk%grownLimits(LOW,  :), &
                  hiGC    => encBlk%grownLimits(HIGH, :), &
                  blkloGC => encBlk%blkLimitsGC(LOW,  :), &
                  blkhiGC => encBlk%blkLimitsGC(HIGH, :))
            call assertEqual(NXB, hi(IAXIS) - lo(IAXIS) + 1, &
                             "Enclosing block has wrong X length")
            call assertEqual(NYB, hi(JAXIS) - lo(JAXIS) + 1, &
                             "Enclosing block has wrong Y length")
            call assertEqual(NZB, hi(KAXIS) - lo(KAXIS) + 1, &
                             "Enclosing block has wrong Z length")

            call assertEqual(NXB+K1D*2*NGUARD, hiGC(IAXIS) - loGC(IAXIS) + 1, &
                             "Enclosing block has wrong X length")
            call assertEqual(NYB+K2D*2*NGUARD, hiGC(JAXIS) - loGC(JAXIS) + 1, &
                             "Enclosing block has wrong Y length")
            call assertEqual(NZB+K3D*2*NGUARD, hiGC(KAXIS) - loGC(KAXIS) + 1, &
                             "Enclosing block has wrong Z length")

            ! If it is a block, then blkLimitsGC = grownLimits
            call assertTrue(ALL(encBlk%blkLimitsGC(LOW,  :) == encBlk%grownLimits(LOW, :)), &
                            "enclosing blocks blkLimitsGC low != grownLimits low")
            call assertTrue(ALL(encBlk%blkLimitsGC(HIGH, :) == encBlk%grownLimits(HIGH, :)), &
                            "enclosing blocks blkLimitsGC high != grownLimits high")

            ! The blkLimitsGC of the block should be grownLimits of enclosing block
            call assertTrue(ALL(tileDesc%blkLimitsGC(LOW,  :) == encBlk%grownLimits(LOW, :)), &
                            "tile blkLimitsGC low != enc block grownLimits low")
            call assertTrue(ALL(tileDesc%blkLimitsGC(HIGH, :) == encBlk%grownLimits(HIGH, :)), &
                            "tile blkLimitsGC high != enc block grownLimits high")
        end associate

        ! Confirm that the UNK data pointer is for the FAB of the enclosing
        ! block and including the block's GC.  It should not be set relative to
        ! the tile size
        call tileDesc%getDataPtr(solnData, CENTER)
        call assertTrue(associated(solnData), "Unable to get cell-centered data")
        bnd = lbound(solnData)
        call assertTrue( ALL(bnd(1:MDIM) == tileDesc%blkLimitsGC(LOW,  1:MDIM)), &
                         "lbound of data pointer is wrong size")
        bnd = ubound(solnData)
        call assertTrue( ALL(bnd(1:MDIM) == tileDesc%blkLimitsGC(HIGH, 1:MDIM)), &
                         "ubound of data pointer is wrong size")
        call tileDesc%releaseDataPtr(solnData, CENTER)

        ! The flux multifabs do not contain GC data.  Confirm that the FAB
        ! has the size of the enclosing block rather than the size of the
        ! tile
        call tileDesc%getDataPtr(solnData, FLUXX)
        call assertTrue(associated(solnData), "Unable to get FACEX data")
        bnd = ubound(solnData) - lbound(solnData) + 1
        call assertEqual(NXB+1, bnd(IAXIS), "Wrong flux FAB X size")
        call assertEqual(NYB,   bnd(JAXIS), "Wrong flux FAB Y size")
        call assertEqual(NZB,   bnd(KAXIS), "Wrong flux FAB Z size")
        call tileDesc%releaseDataPtr(solnData, FLUXX)

        cnt = cnt + 1

        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)
    call assertEqual(cnt, 16, "Incorrect number of tiles")

    ! Confirm proper setting of data
    cnt = 0
    call Grid_getTileIterator(itor, ALL_BLKS, tiling=.TRUE.)
    do while(itor%isValid())
        call itor%currentTile(tileDesc)
        call tileDesc%getDataPtr(solnData, CENTER)

        associate(lo => tileDesc%limits(LOW,  :), &
                  hi => tileDesc%limits(HIGH, :))
            do           var = UNK_VARS_BEGIN, UNK_VARS_END
                do         k = lo(KAXIS), hi(KAXIS)
                    do     j = lo(JAXIS), hi(JAXIS)
                        do i = lo(IAXIS), hi(IAXIS)
                            call assertEqual(solnData(i, j, k, var), &
                                             DBLE(i + j + k)*var, &
                                             "Bad data")
                        end do
                    end do
                end do
            end do
        end associate

        call tileDesc%releaseDataPtr(solnData, CENTER)

        cnt = cnt + 1

        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)
    call assertEqual(cnt, 16, "Incorrect number of tiles")
 
    ! *NO* TILING ITERATION
    ! -----------------------------------------------------------------------
    ! OVERLOAD RUNTIME PARAMETERS ! This is not an expected use, but it is available
    gr_enableTiling = .FALSE. 

    ! Get total blocks by disabling all tiling.  This should ignore
    ! our tiling parameter value
    cnt = 0
    call Grid_getTileIterator(itor, ALL_BLKS, tiling=.TRUE.)
    do while(itor%isValid())
        cnt = cnt + 1
        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)
    call assertEqual(cnt, 2, "Incorrect number of blocks")

    call finish_test_run

end subroutine Driver_evolveFlash


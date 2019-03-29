!!****if* source/Grid/GridMain/AMR/Amrex/Grid_zeroFluxData
!!
!! NAME
!!  Grid_zeroFluxData
!!
!! SYNOPSIS
!!  call Grid_zeroFluxData
!!
!! DESCRIPTION 
!!  Request that the Grid unit zero all flux data managed by the unit.
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Grid_zeroFluxData()
  use Grid_interface,       ONLY : Grid_getTileIterator, &
                                   Grid_releaseTileIterator
  use Grid_iterator,        ONLY : Grid_iterator_t
  use Grid_tile,            ONLY : Grid_tile_t

  implicit none

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc

  real, pointer :: fluxData(:, :, :, :)
  integer       :: i, j, k, var

  nullify(fluxData)

  if (NFLUXES < 1)   RETURN
 
  call Grid_getTileIterator(itor, ALL_BLKS, tiling=.TRUE.)
  do while (itor%isValid())
    call itor%currentTile(tileDesc)

    associate(lo => tileDesc%limits(LOW,  :), &
              hi => tileDesc%limits(HIGH, :))
       call tileDesc%getDataPtr(fluxData, FLUXX)
       if (associated(fluxData)) then
          do        var = 1,         NFLUXES
             do       k = lo(KAXIS), hi(KAXIS)
                do    j = lo(JAXIS), hi(JAXIS)
                   do i = lo(IAXIS), hi(IAXIS)+1
                      fluxData(i, j, k, var) = 0.0
                   end do
                end do
             end do
          end do
       end if
       call tileDesc%releaseDataPtr(fluxData, FLUXX)

       call tileDesc%getDataPtr(fluxData, FLUXY)
       if (associated(fluxData)) then
          do        var = 1,         NFLUXES
             do       k = lo(KAXIS), hi(KAXIS)
                do    j = lo(JAXIS), hi(JAXIS)+1
                   do i = lo(IAXIS), hi(IAXIS)
                      fluxData(i, j, k, var) = 0.0
                   end do
                end do
             end do
          end do
       end if
       call tileDesc%releaseDataPtr(fluxData, FLUXY)

       call tileDesc%getDataPtr(fluxData, FLUXZ)
       if (associated(fluxData)) then
          do        var = 1,         NFLUXES
             do       k = lo(KAXIS), hi(KAXIS)+1
                do    j = lo(JAXIS), hi(JAXIS)
                   do i = lo(IAXIS), hi(IAXIS)
                      fluxData(i, j, k, var) = 0.0
                   end do
                end do
             end do
          end do
       end if
       call tileDesc%releaseDataPtr(fluxData, FLUXZ)
    end associate

    call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
end subroutine Grid_zeroFluxData


!!****if* source/Grid/GridMain/gr_findMean
!!
!! NAME
!!  gr_findMean
!! 
!! SYNOPSIS
!!  call gr_findMean(integer(in)  :: iSrc,
!!                   integer(in)  :: iType,
!!                   logical(in)  :: bGuardcell,
!!                      real(out) :: mean)
!!
!! DESCRIPTION
!!  Calculates the mean of a function
!!
!!
!! ARGUMENTS
!!  iSrc -- the index (e.g. DENS_VAR) into the unk routine to calculate over
!!  iType -- the type of mean.  Valid values will be  (feel free to implement more)
!!     2 = arithmetic average
!!  bGuardcell -- logical indicating whether guard cells should be included in the calculation
!!  mean -- the requested mean
!!
!!
!!***

!!REORDER(4):solnData

#include "constants.h"
#include "Flash.h"

subroutine gr_findMean(iSrc, iType, bGuardcell, mean)
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : Grid_getTileIterator, &
                               Grid_releaseTileIterator
  use Grid_data,        ONLY : gr_meshComm
  use Grid_iterator,    ONLY : Grid_iterator_t
  use Grid_tile,        ONLY : Grid_tile_t

  implicit none

  integer, intent(in)  :: iSrc
  integer, intent(in)  :: iType
  logical, intent(in)  :: bGuardcell
  real,    intent(out) :: mean

#include "Flash_mpi.h"

  real :: localVolume, localSum, blockSum, sum
  real :: numZonesInv

  real :: blockVolume, cellVolume, volume
  integer :: i, j, k, ierr
  integer :: ili, iui, jli, jui, kli, kui
  integer, dimension(2,MDIM) :: limits
  real, dimension(MDIM) :: blockSize
  real, dimension(:,:,:,:), pointer :: solnData
  integer :: nxbBlock, nybBlock, nzbBlock
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc
!!==============================================================================

  nullify(solnData)

  mean = 0.0

  if (iType /= 2) then
     call Driver_abortFlash("[gr_findMean] Can only do arithmetic mean with iType=1")
  else if (bGuardCell) then
     call Driver_abortFlash("[gr_findMean] Inclusion of guardcells not yet coded")
  end if

  localVolume = 0.
  localSum = 0.

  call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
  do while (itor%isValid())
     call itor%currentTile(tileDesc)

     limits = tileDesc%limits
     call tileDesc%getDataPtr(solnData, CENTER)

     ! DEVNOTE: This appears to be for Cartesian only blocks/cells.
     !          Do we need to calculate physical volume or could cell/block
     !          count could be used?
     call tileDesc%physicalSize(blockSize)

     blockVolume = blockSize(IAXIS)
     nxbBlock = limits(HIGH,IAXIS) - limits(LOW,IAXIS) + 1
     numZonesInv = 1. / real(nxbBlock)

     if (NDIM >= 2) then
        blockVolume = blockVolume * blockSize(JAXIS)
        nybBlock = limits(HIGH,JAXIS) - limits(LOW,JAXIS) + 1
        numZonesInv = numZonesInv / real(nybBlock)
     end if

     if (NDIM == 3) then
        blockVolume = blockVolume * blockSize(KAXIS)
        nzbBlock = limits(HIGH,KAXIS) - limits(LOW,KAXIS) + 1
        numZonesInv = numZonesInv / real(nzbBlock)
     end if

     cellVolume = blockVolume * numZonesInv
     localVolume = localVolume + blockVolume
     blockSum = 0.

     ili = limits(LOW,IAXIS)
     iui = limits(HIGH,IAXIS)
     jli = limits(LOW,JAXIS)
     jui = limits(HIGH,JAXIS)
     kli = limits(LOW,KAXIS)
     kui = limits(HIGH,KAXIS)

     do       k = kli, kui
        do    j = jli, jui
           do i = ili, iui
              blockSum = blockSum + solnData(iSrc,i,j,k)
           enddo
        enddo
     enddo

     localSum = localSum + blockSum * cellVolume

     call tileDesc%releaseDataPtr(solnData, CENTER)

     call itor%next()
  enddo
  call Grid_releaseTileIterator(itor)

  call mpi_allreduce ( localSum, sum, 1, FLASH_REAL, & 
       MPI_SUM, gr_meshComm, ierr )
  call mpi_allreduce ( localVolume, volume, 1, FLASH_REAL, & 
       MPI_SUM, gr_meshComm, ierr )

  mean = sum / volume
end subroutine gr_findMean


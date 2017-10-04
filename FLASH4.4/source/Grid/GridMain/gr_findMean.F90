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

subroutine gr_findMean(iSrc, iType, bGuardcell, mean)
  
  use Driver_interface, ONLY: Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, &
                             Grid_getBlkPtr, Grid_releaseBlkPtr
  use Grid_data, ONLY : gr_meshComm
  use block_iterator, ONLY : block_iterator_t, destroy_iterator
  use block_metadata, ONLY : block_metadata_t

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: iSrc, iType
  logical, intent(in) :: bGuardcell
  real, intent(out) :: mean

  real :: localVolume, localSum, blockSum, bvol,  sum
  real :: numZonesInv

  real :: blockVolume, cellVolume, volume
  integer :: i, j, k, ierr
  integer :: ili, iui, jli, jui, kli, kui
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(MDIM) :: blockSize
  real, dimension(:,:,:,:), pointer :: solnData
  integer :: nxbBlock, nybBlock, nzbBlock
  type(block_iterator_t) :: itor
  type(block_metadata_t) :: block
!!==============================================================================

  mean = 0.0

  if (iType /= 2) then
     call Driver_abortFlash("[gr_findMean] Can only do arithmetic mean with iType=1")
  end if
  if (bGuardCell) then
     call Driver_abortFlash("[gr_findMean] Inclusion of guardcells not yet coded")
  end if

  localVolume = 0.
  localSum = 0.

  itor = block_iterator_t(LEAF)
  do while (itor%is_valid())
     call itor%blkMetaData(block)
     
     blkLimits   = block%limits
     blkLimitsGC = block%limitsGC
     call Grid_getBlkPtr(block, solnData)

     ! DEVNOTE: This appears to be for Cartesian only blocks/cells.
     !          Do we need to calculate physical volume or could cell/block
     !          count could be used?
     call Grid_getBlkPhysicalSize(block, blockSize)

     blockVolume = blockSize(IAXIS)
     nxbBlock = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
     numZonesInv = 1. / real(nxbBlock)

     if (NDIM >= 2) then
        blockVolume = blockVolume * blockSize(JAXIS)
        nybBlock = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
        numZonesInv = numZonesInv / real(nybBlock)
     end if

     if (NDIM == 3) then
        blockVolume = blockVolume * blockSize(KAXIS)
        nzbBlock = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1
        numZonesInv = numZonesInv / real(nzbBlock)
     end if


     cellVolume = blockVolume * numZonesInv
     localVolume = localVolume + blockVolume
     blockSum = 0.

     if (bGuardcell) then
        ili = blkLimitsGC(LOW,IAXIS)
        iui = blkLImitsGC(HIGH,IAXIS)
        jli = blkLimitsGC(LOW,JAXIS)
        jui = blkLImitsGC(HIGH,JAXIS)
        kli = blkLimitsGC(LOW,KAXIS)
        kui = blkLImitsGC(HIGH,KAXIS)
     else
        ili = blkLimits(LOW,IAXIS)
        iui = blkLImits(HIGH,IAXIS)
        jli = blkLimits(LOW,JAXIS)
        jui = blkLImits(HIGH,JAXIS)
        kli = blkLimits(LOW,KAXIS)
        kui = blkLImits(HIGH,KAXIS)
     endif

     do k = kli, kui
        do j = jli, jui
           do i = ili, iui
              blockSum = blockSum + solnData(iSrc,i,j,k)
           enddo
        enddo
     enddo

     localSum = localSum + blockSum * cellVolume

     call Grid_releaseBlkPtr(block, solnData)
     nullify(solnData)

     call itor%next()
  enddo
#if defined(__GFORTRAN__) && (__GNUC__ <= 4)
  call destroy_iterator(itor)
#endif

  call mpi_allreduce ( localSum, sum, 1, FLASH_REAL, & 
       MPI_SUM, gr_meshComm, ierr )
  call mpi_allreduce ( localVolume, volume, 1, FLASH_REAL, & 
       MPI_SUM, gr_meshComm, ierr )

  mean = sum / volume

  return

end subroutine gr_findMean


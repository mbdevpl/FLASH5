!!F****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/Grid_updateSolidBodyForces
!!
!! NAME
!!  Grid_updateSolidBodyForces
!!
!! SYNOPSIS
!!
!!  Grid_updateSolidBodyForces(integer, INTENT(in)    :: blkID,
!!                             real, dimension(NPART_PROPS), INTENT(inout) :: particledata)
!!  
!! DESCRIPTION 
!!  
!!  The velocity update and forcing routine
!!
!!  Overview of the algoritm
!!
!! ARGUMENTS
!!
!!  blkID:  the local block ID
!!
!!  particledata : data corresponding to the particle
!! 
!!  
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_updateSolidBodyForces(ibd,p,blkID,particleData)

  use Grid_data, ONLY : gr_meshMe
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox, &
                             Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr

  use ib_interface, only : ib_forcing

  implicit none
  include "Flash_mpi.h"

!! Argument list 
  integer, intent(IN) :: blkID,ibd,p
  real, dimension(NPART_PROPS), intent(INOUT) :: particleData
!-------------------------------------------------------

! Call the ib_forcing subroutine to get the forces on each particle
call ib_forcing(ibd,p,blkID,particleData)
  
  return
end subroutine Grid_updateSolidBodyForces

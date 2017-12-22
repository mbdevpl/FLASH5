!!****if* source/physics/Hydro/HydroMain/Hydro_loop5Body
!!
!!
!! NAME
!!
!!  Hydro_loop5Body
!!
!!
!! SYNOPSIS
!!
!!  Hydro_loop5Body(integer(IN) :: blockCount, 
!!        integer(IN) :: blockList(blockCount)
!!        real(IN)    :: timeEndAdv,
!!        real(IN)    :: dt,
!!        real(IN)    :: dtOld,
!!        integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!! 
!!  Performs physics update in a directionally unsplit fashion.
!!
!!  The blockList and blockCount arguments tell this routine on 
!!  which blocks and on how many to operate.  blockList is an 
!!  integer array of size blockCount that contains the local 
!!  block numbers of blocks on which to advance.
!!
!!  dt gives the timestep through which this update should advance,
!!  and timeEndAdv tells the time that this update will reach when
!!  it finishes.  dtOld gives the previously taken timestep.
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - end time
!!  dt         - timestep
!!  dtOld      - old timestep
!!  sweepOrder - argument for the unsplit scheme, for unsplit Hydro this
!!               just a dummy variable to be consistent with the API.
!!
!!***

#include "constants.h"

Subroutine Hydro_loop5Body(block, blkLimitsGC, Uin, blkLimits, Uout,  del, timeEndAdv, dt,  dtOld)
  use block_metadata, ONLY : block_metadata_t
                        
  implicit none

  real,    INTENT(IN) :: timeEndAdv, dt, dtOld

  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits,blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: Uout
  real,dimension(MDIM),intent(IN) :: del
  real, pointer, dimension(:,:,:,:) :: Uin
  type(block_metadata_t),intent(IN) :: block

  return
End Subroutine Hydro_loop5Body

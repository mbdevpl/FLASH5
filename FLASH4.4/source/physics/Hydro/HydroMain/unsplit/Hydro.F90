!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro
!!
!!
!! NAME
!!
!!  Hydro
!!
!!
!! SYNOPSIS
!!
!!  Hydro(integer(IN) :: blockCount, 
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
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a toplayer stub function
!!
!!***

Subroutine Hydro(Uin,Uout,block,timeEndAdv,dt,dtold,sweepDummy)

  use Hydro_data,       ONLY : hy_useHydro, hy_gpotAlreadyUpToDate
  use hy_uhd_interface, ONLY : hy_uhd_unsplit
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use block_metadata,   ONLY : block_metadata_t

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, INTENT(IN) :: sweepOrder
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld
  type(block_metadata_t), INTENT(IN) :: block

  real, pointer, dimension(:,:,:,:) :: Uout
  real, pointer, dimension(:,:,:,:) :: Uin

  read,dimension(MDIM) :: del
  integer,dimension(LOW:HIGH,MDIM)::blkLimits,blkLimitsGC
  
  hy_gpotAlreadyUpToDate = .FALSE. ! reset this flag, may be set .TRUE. below if warranted.

  if (.not. hy_useHydro) return

  call Timers_start("hydro_unsplit")
  
  del=block%del
  call hy_uhd_unsplit(Uin,blkLimitGC,&
                      Uout,blkLimits,&
                      del,dt, dtOld )

  call Timers_stop("hydro_unsplit")


End Subroutine Hydro

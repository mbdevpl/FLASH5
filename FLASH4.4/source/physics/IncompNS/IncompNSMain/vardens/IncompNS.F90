!!****if* source/physics/IncompNS/IncompNSMain/vardens/IncompNS
!!
!!
!! NAME
!!
!!  IncompNS
!!
!!
!! SYNOPSIS
!!
!!  IncompNS(integer(IN) :: blockCount, 
!!      integer(IN) :: blockList(blockCount)
!!      real(IN)    :: timeEndAdv,
!!      real(IN)    :: dt,
!!      real(IN)    :: dtOld,
!!      integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!! 
!!  Performs INS timestep advancement.
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
!!  numProcs   - total number of processors
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - dummy consistent with toplayer stub function
!!  dt         - timestep
!!  dtOld      - dummy consistent with toplayer stub function
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a toplayer stub function
!!
!!***

subroutine IncompNS( blockCount, blockList, &
                     timeEndAdv, dt,  dtOld,&
                     sweepOrder)

  use ins_interface, ONLY : ins_ab2rk3, ins_ab2rk3_VD
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Driver_data,   ONLY : dr_nstep

  use IncompNS_data, ONLY : ins_nstep,ins_meshMe

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, INTENT(IN)    :: sweepOrder
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT) :: blockList(MAXBLOCKS) !blockCount)
  real,    INTENT(IN)    :: timeEndAdv, dt, dtOld
  integer :: x, i

  ins_nstep = dr_nstep

  call Timers_start("ins_ab2rk3")

  !call ins_ab2rk3(  blockCount, blockList, dt)
  !call ins_ab2rk3(  blockCount, blockList, timeEndAdv, dt)
  if (ins_meshMe .eq. 0) print*,"Using Variable Density INS (ins_ab2rk3_VD)..."
  call ins_ab2rk3_VD(  blockCount, blockList, timeEndAdv, dt)

  call Timers_stop("ins_ab2rk3")


end subroutine IncompNS

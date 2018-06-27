!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/ImBound
!!
!!
!! NAME
!!
!!  ImBound
!!
!!
!! SYNOPSIS
!!
!!  ImBound(blockCount,blockList,dt)
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***

#include "Flash.h"

subroutine ImBound(blockCount, blockList, dt, forcflag)

  use Driver_data, only : dr_simTime
  use Driver_interface, only : Driver_abortFlash
  use Grid_interface, only : Grid_sbSelectMaster
  use gr_sbInterface, only : gr_sbCreateParticles,gr_sbGetProcBlock, &
                             gr_sbSendPosn, gr_sbDistributedForces,  &
                             gr_sbFinalize

  use ImBound_data, only : ib_dt,ib_BlockMarker_flag

  use gr_sbData, only : gr_sbNumBodies,gr_sbBodyInfo,gr_sbFirstCall

  use Logfile_interface,ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop

#ifdef FORCES_FROM_VOL_INTEGRALS
  use sm_assemble_interface, only : sm_assemble_FluidVolIntegrals
#endif
  implicit none
#include "ImBound.h"
#include "constants.h"

  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
  real, INTENT(IN) :: dt
  integer, INTENT(IN) :: forcflag
  !! -----------------------------------------------------


  ib_dt = dt

  select case (forcflag)
  case(FORCE_FLOW) 

  ! Set Refinement to False in all blocks
  ib_BlockMarker_flag = .FALSE.

  call Logfile_stamp( 'Entering Forcing' , '[ImBound]')

  ! Flush Particles.
  call gr_sbFinalize()

  ! First, select Master:
  !call Timers_start("Grid_sbSelectMaster")
  !call Grid_sbSelectMaster()
  !call Timers_stop("Grid_sbSelectMaster")

  call Timers_start("gr_sbCreateParticles")
  call gr_sbCreateParticles()
  call Timers_stop("gr_sbCreateParticles")

  call Timers_start("gr_sbGetProcBlock")
  call gr_sbGetProcBlock()
  call Timers_stop("gr_sbGetProcBlock")

  call Timers_start("gr_sbSendPosn")
  call gr_sbSendPosn()
  call Timers_stop("gr_sbSendPosn")

  gr_sbFirstCall = CONSTANT_ZERO

  call Logfile_stamp( 'After Forcing' , '[ImBound]')

  case(COMPUTE_FORCES)

  call gr_sbDistributedForces()

#ifdef FORCES_FROM_VOL_INTEGRALS
  call sm_assemble_FluidVolIntegrals()
#endif  
  case default

  call Driver_abortFlash("ImBound : forcflag doen not correspond to any available option.") 
  
  end select

  return
end subroutine ImBound


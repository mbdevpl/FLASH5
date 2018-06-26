!!****if* source/physics/ImBound/ImBoundMain/LagForce/parallel/ImBound_init
!!
!! NAME
!!
!!  ImBound_init
!!
!!
!! SYNOPSIS
!!
!!  call ImBound_init(restart)
!!  
!! VARIABLES
!!
!! restart = restart flag, logical.
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine ImBound_init(restart)

  use Driver_interface, ONLY: Driver_getMype, Driver_getComm
  use ImBound_data, ONLY: ib_meshMe, ib_globalMe, ib_meshComm, ib_maxIterForcing
  use gr_ptVPData, ONLY : gr_ptVPBufferFactor
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
#include "constants.h"
  logical, INTENT(IN) :: restart

  call Driver_getMype(MESH_COMM,ib_meshMe)
  call Driver_getComm(MESH_COMM,ib_meshComm)
  call Driver_getMype(GLOBAL_COMM,ib_globalMe)

  ! Set Buffer factor for Virtual Particle creation to 1.5:
  gr_ptVPBufferFactor = 1.5

  ! Get the maximum iterations of the IB forcing
  call RuntimeParameters_get("ib_maxIterForcing", ib_maxIterForcing)

  if(ib_meshMe == MASTER_PE) write(*,*) 'Maximum iterations in IB forcing', &
                        &    ib_maxIterForcing 

  if(ib_globalMe==MASTER_PE)print*,'Immersed_Boundaries initialized'
  return

end subroutine ImBound_init

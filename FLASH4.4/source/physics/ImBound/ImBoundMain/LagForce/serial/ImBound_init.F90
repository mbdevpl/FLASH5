!!****if* source/physics/ImBound/ImBoundMain/LagForce/serial/ImBound_init
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
!! ARGUMENTS
!!
!!  restart - restart flag, logical.
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
  use ImBound_data, ONLY: ib_meshMe, ib_globalMe, ib_meshComm
  implicit none
#include "constants.h"
  logical, INTENT(IN) :: restart

  call Driver_getMype(MESH_COMM,ib_meshMe)
  call Driver_getComm(MESH_COMM,ib_meshComm)
  call Driver_getMype(GLOBAL_COMM,ib_globalMe)
  if(ib_globalMe==MASTER_PE)print*,'Immersed_Boundaries initialized'
  return

end subroutine ImBound_init

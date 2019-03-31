!!****f* source/physics/Hydro/Hydro
!!
!!
!! NAME
!!
!!  Hydro
!!
!! SYNOPSIS
!!
!!  Hydro( real(IN)    :: timeEndAdv, 
!!         real(IN)    :: dt, 
!!         real(IN)    :: dtOld, 
!!         integer(IN) :: sweepOrder )
!!
!! DESCRIPTION
!!
!!  Perform a 1, 2, or 3D hydro update.  This version handles
!!  directionally split hydro schemes.  The input
!!  parameter sweepOrder determines the ordering of sweep
!!  directions.  For example, in 3d, SWEEP_XYZ means to perform 
!!  1d sweeps first in the x direction, then the y, then z direction, 
!!  while SWEEP_ZYX means to performs the 1d sweeps in the opposite order.
!!
!!  dt gives the timestep through which this update should advance,
!!  and timeEndAdv tells the time that this update will reach when
!!  it finishes.  dtOld gives the previously taken timestep.  
!!
!! ARGUMENTS
!!
!!  timeEndAdv -  end time
!!  dt -          timestep
!!  dtOld -       old timestep
!!  sweepOrder -  direction of hydro sweep, can be: SWEEP_XYZ or SWEEP_ZYX
!!                as defined in  constants.h
!!
!!
!!***


subroutine Hydro(simTime, dt, dtOld, sweeporder)
implicit none
#include "Flash.h"
#include "constants.h"
  
  real,    INTENT(IN) :: simTime, dt, dtOld
  integer, optional, intent(IN) :: sweeporder
end subroutine Hydro

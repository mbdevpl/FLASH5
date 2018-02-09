!!****if* source/Simulation/SimulationMain/Sedov/Simulation_adjustEvolution
!!
!! NAME
!!  Simulation_adjustEvolution
!!
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: blkcnt,
!!                              integer(IN) :: blklst(blkcnt),
!!                              integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  
!! ARGUMENTS
!!  blkcnt - number of blocks
!!  blklst - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***

#include "Flash.h"

subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)
  use Simulation_interface, ONLY: Simulation_computeAnalytical
  use Grid_interface,   ONLY : Grid_dump
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Driver_data, ONLY: dr_simTime, dr_initialSimTime
  implicit none

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real, intent(in) :: dt
  real, intent(in) :: stime

#if NDIM==1
  real    :: tnew
  integer :: lb

  call Timers_start("adjustEvo")

  if (dr_simTime .LE. dr_initialSimTime) then
     do lb = 1, blkcnt
        call Grid_dump((/DENS_VAR,PRES_VAR,VELX_VAR,EINT_VAR/),4, blklst(lb),gcell=.FALSE.)
     end do
  end if


!!$  ! We want the NEW time that we are advancing TO.
!!$  ! This is based on where Simulation_adjustEvolution is called from Driver_evolveFlash:
!!$  ! dr_simTime has not been updated yet at this point.
!!$  !!DEV: Check whether this needs adjustments for non-"Unsplit" Driver_evolveFlash.F90 !
!!$  !!DEV: Check whether this needs adjustments for when STS is used!
!!$  tnew = stime+dt
!!$
!!$  do lb = 1, blkcnt
!!$     call Simulation_computeAnalytical(blklst(lb), tnew)
!!$  end do
  call Timers_stop("adjustEvo")
#endif

end subroutine Simulation_adjustEvolution

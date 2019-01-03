!!****if* source/physics/Hydro/localAPI/hy_advanceBlk
!!
!! NAME
!!  hy_advanceBlk
!!
!! SYNOPSIS
!!  hy_advanceBlk(flash_tile_t(IN) :: tileDesc, 
!!                real, pointer    :: Uin,
!!                real, pointer    :: Uout,
!!                real(IN)         :: timeEndAdv,
!!                real(IN)         :: dt,
!!                real(IN)         :: dtOld,
!!                integer(IN)      :: sweepOrder)
!!
!! DESCRIPTION
!!  Performs physics update in a directionally unsplit fashion on the given
!!  tile.
!!
!!  dt gives the timestep through which this update should advance,
!!  and timeEndAdv tells the time that this update will reach when
!!  it finishes.  dtOld gives the previously taken timestep.
!!
!! ARGUMENTS
!!  tileDesc   - the tile on which to advance the solution
!!  Uin        - pointer to original solution data
!!  Uout       - pointer to advanced solution
!!  timeEndAdv - end time
!!  dt         - timestep
!!  dtOld      - old timestep
!!  sweepOrder - argument for the unsplit scheme, for unsplit Hydro this
!!               just a dummy variable to be consistent with the API.
!!
!!***

#include "UHD.h"

Subroutine hy_advanceBlk(tileDesc, Uin, Uout, timeEndAdv, dt, dtOld, sweepOrder)
  use flash_tile,       ONLY : flash_tile_t

  implicit none

  type(flash_tile_t), intent(IN)         :: tileDesc
  real,                          pointer :: Uout(:,:,:,:)
  real,                          pointer :: Uin(:,:,:,:)
  real,               intent(IN)         :: timeEndAdv
  real,               intent(IN)         :: dt
  real,               intent(IN)         :: dtOld
  integer,            intent(IN)         :: sweeporder

  RETURN
End Subroutine hy_advanceBlk


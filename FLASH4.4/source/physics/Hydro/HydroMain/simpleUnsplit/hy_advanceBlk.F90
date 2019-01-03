!!****if* source/physics/Hydro/HydroMain/simpleUnsplit/hy_advanceBlk
!!
!! NAME
!!  hy_advanceBlk
!!
!! DESCRIPTION
!!  Refer to stub for detailed documentation.
!!
!!***

#include "UHD.h"

Subroutine hy_advanceBlk(tileDesc, Uin, Uout, timeEndAdv, dt, dtOld, sweepOrder)
  use Hydro_data,       ONLY : hy_useHydro, &
                               hy_riemannSolver, &
                               hy_eosModeAfter
  use hy_interface,     ONLY : hy_hllUnsplit, hy_llfUnsplit
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface,    ONLY : Eos_wrapped
  use flash_tile,       ONLY : flash_tile_t

  implicit none

  type(flash_tile_t), intent(IN)         :: tileDesc
  real,                          pointer :: Uout(:,:,:,:)
  real,                          pointer :: Uin(:,:,:,:)
  real,               intent(IN)         :: timeEndAdv
  real,               intent(IN)         :: dt
  real,               intent(IN)         :: dtOld
  integer,            intent(IN)         :: sweeporder

  real :: deltas(1:MDIM)

  if (.not. hy_useHydro) return 

  call tileDesc%deltas(deltas)

  select case (hy_riemannSolver)
  case(HLL)
     call hy_hllUnsplit(tileDesc%limits, Uin, lbound(Uin), Uout, deltas, dt)
  ! DEV: FIXME LLF has not been updated yet
!  case(LLF)
!     call hy_llfUnsplit(tileDesc%limits, Uin, Uout, deltas, dt)
  case default
     call Driver_abortFlash("[hy_advanceBlk]: Unknown Riemann solver")
  end select

  call Eos_wrapped(hy_eosModeAfter, tileDesc%limits, Uout)

End Subroutine hy_advanceBlk


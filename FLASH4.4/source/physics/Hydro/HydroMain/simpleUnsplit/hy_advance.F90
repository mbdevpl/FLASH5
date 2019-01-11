!!****if* source/physics/Hydro/HydroMain/simpleUnsplit/hy_advance
!!
!! NAME
!!  hy_advance
!!
!! DESCRIPTION
!!  Refer to stub for detailed documentation.
!!
!!***

#include "UHD.h"

subroutine hy_advance(simTime, dt, dtOld)
  use Grid_interface,   ONLY : Grid_getTileIterator, &
                               Grid_releaseTileIterator
  use Timers_interface, ONLY : Timers_start, &
                               Timers_stop
  use Hydro_data,       ONLY : hy_useHydro, &
                               hy_riemannSolver, &
                               hy_eosModeAfter
  use hy_interface,     ONLY : hy_hllComputeFluxes, &
                               hy_hllUpdateSolution
  use Eos_interface,    ONLY : Eos_wrapped
  use flash_iterator,   ONLY : flash_iterator_t
  use flash_tile,       ONLY : flash_tile_t

  implicit none

  real, intent(IN) :: simTime
  real, intent(IN) :: dt
  real, intent(IN) :: dtOld
 
  real, pointer :: Uout(:,:,:,:)
  real, pointer :: Uin(:,:,:,:)
  real, pointer :: flX(:,:,:,:)
  real, pointer :: flY(:,:,:,:)
  real, pointer :: flZ(:,:,:,:)

  type(flash_iterator_t) :: itor
  type(flash_tile_t)     :: tileDesc

  real    :: deltas(1:MDIM)

  nullify(Uin)
  nullify(Uout)
  nullify(flX)
  nullify(flY)
  nullify(flZ)

  if (.not. hy_useHydro) return 

  call Timers_start("loop1")
  call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     call tileDesc%deltas(deltas)
     call tileDesc%getDataPtr(flX, FLUXX)
     call tileDesc%getDataPtr(flY, FLUXY)
     call tileDesc%getDataPtr(flZ, FLUXZ)
     call tileDesc%getDataPtr(Uin, CENTER)

     select case (hy_riemannSolver)
     case(HLL)
        call hy_hllComputeFluxes(tileDesc%limits, &
                                 Uin, lbound(Uin), &
                                 flX, flY, flZ, lbound(flX), &
                                 deltas, dt)
     case default
        call Driver_abortFlash("[hy_advance]: Unknown Riemann solver")
     end select

     call tileDesc%releaseDataPtr(Uin, CENTER)
     call tileDesc%releaseDataPtr(flX, FLUXX)
     call tileDesc%releaseDataPtr(flY, FLUXY)
     call tileDesc%releaseDataPtr(flZ, FLUXZ)

     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     call tileDesc%deltas(deltas)

     call tileDesc%getDataPtr(flX, FLUXX)
     call tileDesc%getDataPtr(flY, FLUXY)
     call tileDesc%getDataPtr(flZ, FLUXZ)
     call tileDesc%getDataPtr(Uin, CENTER)
     Uout => Uin

     select case (hy_riemannSolver)
     case(HLL)
        call hy_hllUpdateSolution(tileDesc%limits, &
                                  Uin, lbound(Uin), Uout, &
                                  flX, flY, flZ, lbound(flX), &
                                  deltas, dt)
     case default
        call Driver_abortFlash("[hy_advance]: Unknown Riemann solver")
     end select

     call Eos_wrapped(hy_eosModeAfter, tileDesc%limits, Uout)

     call tileDesc%releaseDataPtr(Uin, CENTER)
     call tileDesc%releaseDataPtr(flX, FLUXX)
     call tileDesc%releaseDataPtr(flY, FLUXY)
     call tileDesc%releaseDataPtr(flZ, FLUXZ)
     nullify(Uout)

     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
  call Timers_stop("loop1")

#ifdef DEBUG_DRIVER
  print*, 'return from Hydro/MHD timestep'  ! DEBUG
  print*,'returning from hydro myPE=',dr_globalMe
#endif

end subroutine hy_advance


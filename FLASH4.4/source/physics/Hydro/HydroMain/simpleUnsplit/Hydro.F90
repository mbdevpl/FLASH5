!!****f* source/physics/Hydro/HydroMain/simpleUnsplit/Hydro
!!
!! NAME
!!  Hydro
!!
!! DESCRIPTION
!!  Please see stub version of this file for detailed information.
!!
!!  For the simpleUnsplit implementation of Hydro, we compute all fluxes in one
!!  pass and then update solution based on solution of last iteration and these
!!  fluxes.
!!
!!  This two-loop structure allows for saving the fluxes as an intermediate
!!  result in a different Grid-owned data structure and then updating the solution
!!  in-place.  This is necessitated by tiling so that we aren't overwriting a
!!  tile's worth of the solution data before computing the fluxes on neighboring
!!  tiles.
!!
!!***

#define DEBUG_GRID_GCMASK

#include "UHD.h"

subroutine Hydro(simTime, dt, dtOld)
  use Grid_interface,    ONLY : Grid_fillGuardCells, &
                                Grid_getTileIterator, &
                                Grid_releaseTileIterator
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Timers_interface,  ONLY : Timers_start, &
                                Timers_stop
  use Eos_interface,     ONLY : Eos_wrapped
  use Hydro_interface,   ONLY : Hydro_prepareBuffers, &
                                Hydro_freeBuffers
  use Hydro_data,        ONLY : hy_useHydro, &
                                hy_riemannSolver, &
                                hy_eosModeAfter, &
                                hy_gcMaskSize, &
                                hy_gcMask
  use hy_interface,      ONLY : hy_hllComputeFluxes, &
                                hy_hllUpdateSolution
  use flash_iterator,    ONLY : flash_iterator_t
  use flash_tile,        ONLY : flash_tile_t

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld
  
  real, pointer :: Uout(:,:,:,:)
  real, pointer :: Uin(:,:,:,:)
  real, pointer :: flX(:,:,:,:)
  real, pointer :: flY(:,:,:,:)
  real, pointer :: flZ(:,:,:,:)

  type(flash_iterator_t) :: itor
  type(flash_tile_t)     :: tileDesc

  real :: deltas(1:MDIM)

#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif

  nullify(Uin)
  nullify(Uout)
  nullify(flX)
  nullify(flY)
  nullify(flZ)

  if (.NOT. hy_useHydro) RETURN

  call Timers_start("Hydro")

  call Hydro_prepareBuffers()

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(hy_gcMask, .TRUE., '[Hydro]', 'gcNeed')
  end if
#endif
 
  !!!!!----- OBTAIN ALL GC DATA FOR USE WITH STENCIL
  ! Assume that the data is good only on the block interiors
  call Grid_fillGuardCells(CENTER, ALLDIR, &
                           doEos=.TRUE., &
                           maskSize=hy_gcMaskSize, &
                           mask=hy_gcMask, &
                           makeMaskConsistent=.TRUE., &
                           doLogMask=.NOT.gcMaskLogged)

  ! DEV: Should shock detection be done here?

  !!!!!----- COMPUTE FLUXES ON ALL LEAF BLOCKS INTERIORS
  call Timers_start("compute fluxes")
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
  call Timers_stop("compute fluxes")

  !!!!!----- COMPUTE SOLUTIONS ON ALL LEAF BLOCK INTERIORS
  call Timers_start("update solution")
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
  call Timers_stop("update solution")

  call Hydro_freeBuffers()

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
#endif

  call Timers_stop("Hydro")

end subroutine Hydro


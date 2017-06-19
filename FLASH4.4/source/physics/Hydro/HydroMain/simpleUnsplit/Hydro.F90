!!****if* source/physics/Hydro/HydroMain/simpleUnsplit/Hydro
!!
!!
!! NAME
!!
!!  Hydro
!!
!!
!! SYNOPSIS
!!
!!  Hydro(integer(IN) :: blockCount, 
!!        integer(IN) :: blockList(blockCount)
!!        real(IN)    :: timeEndAdv,
!!        real(IN)    :: dt,
!!        real(IN)    :: dtOld,
!!        integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!! 
!!  Performs physics update in a directionally unsplit fashion.
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
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - end time
!!  dt         - timestep
!!  dtOld      - old timestep
!!  sweepOrder - argument for the unsplit scheme, for unsplit Hydro this
!!               just a dummy variable to be consistent with the API.
!!
!!***

Subroutine Hydro( del, tileLimits, Uout, timeEndAdv, dt,  dtOld,&
                  sweepOrder)

  use Hydro_data,       ONLY : hy_useHydro, hy_riemannSolver
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getListOfBlocks,    &
                             Grid_getBlkPtr,          &
                             Grid_getBlkIndexLimits,  &
                             Grid_releaseBlkPtr

  use hy_simpleInterface, ONLY : hy_hllUnsplit, hy_llfUnsplit
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_wrapped
  use Hydro_data, ONLY :hy_meshComm
  use Hydro_data, ONLY :hy_gcMaskSize,       &
                         hy_gcMask,           &
                         hy_unsplitEosMode,   &
                         hy_eosModeAfter
                        

  implicit none

#include "UHD.h"

  integer, INTENT(IN) :: sweepOrder
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld

  real,dimension(MDIM),intent(IN) :: del
  real, pointer, dimension(:,:,:,:) :: Uout
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: tileLimits
  real, pointer, dimension(:,:,:,:) :: Uin



  if (.not. hy_useHydro) return 

!!  call Timers_start("hydro_sUnsplit")

  Uin => Uout
  
  select case (hy_riemannSolver)
  case(HLL)
     call hy_hllUnsplit(tileLimits, Uin, lbound(Uin), Uout, del, dt)
  case(LLF)
     call hy_llfUnsplit(tileLimits, Uin, Uout, del, dt)
  case default
     call Driver_abortFlash("Hydro: what?")
  end select
  
  !! Call to Eos - note this is a variant where we pass a buffer not a blockID.
  call Eos_wrapped(hy_eosModeAfter, tileLimits, Uout)


!!$#ifdef DEBUG_GRID_GCMASK
!!$  if (.NOT.gcMaskLogged) then
!!$     gcMaskLogged = .TRUE.
!!$  end if
!!$#endif

!!$  select case (hy_riemannSolver)
!!$  case(HLL)
!!$     call hy_hllUnsplit(blockCount, blockList, dt, dtOld)
!!$  case(LLF)
!!$     call hy_llfUnsplit(blockCount, blockList, dt, dtOld)
!!$  case default
!!$     call Driver_abortFlash("Hydro: what?")
!!$  end select

!!  call Timers_stop("hydro_sUnsplit")


End Subroutine Hydro

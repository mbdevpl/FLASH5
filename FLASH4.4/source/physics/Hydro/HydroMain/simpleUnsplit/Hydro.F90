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

Subroutine Hydro( blockCount, blockList, &
                  timeEndAdv, dt,  dtOld,&
                  sweepOrder)

  use Hydro_data,       ONLY : hy_useHydro, hy_riemannSolver
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_genGetBlkPtr,         &
                             Grid_genReleaseBlkPtr,     &
                             Grid_getBlkData


  use hy_simpleInterface, ONLY : hy_hllUnsplit, hy_llfUnsplit
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_wrapped
  use Hydro_data, ONLY :hy_gcMaskSize,       &
                         hy_gcMask,           &
                         hy_unsplitEosMode,   &
                         hy_eosModeAfter
                        

  implicit none

#include "UHD.h"

  integer, INTENT(IN) :: blockCount,sweepOrder
  integer, INTENT(IN) :: blockList(blockCount)
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld

  real,dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: Uin,Uout
  integer, dimension(LOW:HIGH,MDIM) :: tileLimits ,blkLimitsGC
  integer :: ib,blockID
  logical :: gcMask(hy_gcMaskSize)
#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif


  if (.not. hy_useHydro) return 

  call Timers_start("hydro_sUnsplit")

!!ChageForAMRex -- Here is where we put in the iterator and extract the relevant metadata
!!ChageForAMRex -- from the iterator and then use the case statement to transfer control to the
!!ChageForAMRex -- right implementation.

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(hy_gcMask, .FALSE., '[hy_hllUnsplit]', 'gcNeed')
  end if
#endif

  !! Guardcell filling routine
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=hy_gcMaskSize, mask=hy_gcMask,makeMaskConsistent=.true.,doLogMask=.NOT.gcMaskLogged)

  do ib=1,blockCount

     blockID = blockList(ib)

!!ChageForAMRex -- this information should come from the interator as meta-data
     call Grid_getDeltas(blockID,del)

!!$     dtdx = dt / del(IAXIS)
!!$     if (NDIM > 1) dtdy = dt / del(JAXIS)
!!$     if (NDIM > 2) dtdz = dt / del(KAXIS)

!!ChageForAMRex -- this information should come from the iterator as meta-data
     call Grid_getBlkIndexLimits(blockID,tileLimits,blkLimitsGC)

     call Grid_getBlkPtr(blockID,Uout,CENTER)
     Uin => Uout

     select case (hy_riemannSolver)
     case(HLL)
        call hy_hllUnsplit(tileLimits, Uin, Uout, del, dt)
     case(LLF)
        call hy_llfUnsplit(tileLimits, Uin, Uout, del, dt)
     case default
        call Driver_abortFlash("Hydro: what?")
     end select

        !! Call to Eos - note this is a variant where we pass a buffer not a blockID.
     call Eos_wrapped(hy_eosModeAfter, tileLimits, Uout)
     call Grid_releaseBlkPtr(blockID,Uout,CENTER)

  end do

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
#endif

!!$  select case (hy_riemannSolver)
!!$  case(HLL)
!!$     call hy_hllUnsplit(blockCount, blockList, dt, dtOld)
!!$  case(LLF)
!!$     call hy_llfUnsplit(blockCount, blockList, dt, dtOld)
!!$  case default
!!$     call Driver_abortFlash("Hydro: what?")
!!$  end select

  call Timers_stop("hydro_sUnsplit")


End Subroutine Hydro

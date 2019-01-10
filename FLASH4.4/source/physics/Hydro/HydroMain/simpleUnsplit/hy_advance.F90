!!****if* source/physics/Hydro/HydroMain/simpleUnsplit/hy_advance
!!
!! NAME
!!  hy_advance
!!
!! DESCRIPTION
!!  Refer to stub for detailed documentation.
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine hy_advance(simTime, dt, dtOld)
  use Grid_interface,   ONLY : Grid_getTileIterator, &
                               Grid_releaseTileIterator, &
                               Grid_getMaxRefinement
  use Timers_interface, ONLY : Timers_start, &
                               Timers_stop
  use hy_interface,     ONLY : hy_advanceBlk
  use flash_iterator,   ONLY : flash_iterator_t
  use flash_tile,       ONLY : flash_tile_t

  use gr_physicalMultifabs,  ONLY : unk
  use amrex_multifab_module, ONLY : amrex_multifab, &
                                    amrex_multifab_build, &
                                    amrex_multifab_destroy

  implicit none

  real, intent(IN) :: simTime
  real, intent(IN) :: dt
  real, intent(IN) :: dtOld
 
  type(amrex_multifab) :: tmp

  real, pointer :: Uout(:,:,:,:)
  real, pointer :: Uin(:,:,:,:)

  type(flash_iterator_t) :: itor
  type(flash_tile_t)     :: tileDesc

  integer :: level
  integer :: finest_level

  nullify(Uin)
  nullify(Uout)

  call Timers_stop("loop1")
  call Grid_getMaxRefinement(finest_level, mode=3)
  do level = 1, finest_level
     call amrex_multifab_build(tmp, &
                               unk(level-1)%ba, &
                               unk(level-1)%dm, &
                               NUNK_VARS, NGUARD)

     call Grid_getTileIterator(itor, LEAF, level=level, tiling=.TRUE.)
     do while(itor%isValid())
        call itor%currentTile(tileDesc)

        call tileDesc%getDataPtr(Uin, CENTER)
        associate(lo => lbound(Uin))
           Uout(lo(1):, lo(2):, lo(3):, lo(4):) => tmp%dataPtr(tileDesc%grid_index)
        end associate

        call hy_advanceBlk(tileDesc, Uin, Uout, simTime, dt, dtOld, SWEEP_ALL)

        call tileDesc%releaseDataPtr(Uin, CENTER)
        nullify(Uout)

        call itor%next()
     end do
     call Grid_releaseTileIterator(itor)

     call unk(level-1)%copy(tmp, UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, 0)
     call amrex_multifab_destroy(tmp)
  end do
  call Timers_stop("loop1")

#ifdef DEBUG_DRIVER
  print*, 'return from Hydro/MHD timestep'  ! DEBUG
  print*,'returning from hydro myPE=',dr_globalMe
#endif

end subroutine hy_advance


#include "Flash.h"

subroutine Hydro_prepareBuffers
  use Grid_interface,  ONLY : Grid_getListOfBlocks
  use hy_memInterface, ONLY :  hy_memAllocScratch
  use Hydro_data, ONLY : hy_fullRiemannStateArrays
  implicit none

#include "constants.h"
#include "UHD.h"

#ifndef FLASH_GRID_ANYAMREX

  integer :: blockList(MAXBLOCKS)
  integer :: blockCount


  call Grid_getListOfBlocks(LEAF, blockList, blockCount)


  ! Independently of whether hy_fullRiemannStateArrays is set:
  call hy_memAllocScratch(SCRATCH_CTR,HY_VAR1_SCRATCHCTR_VAR,2, 0,0,0, &
       blockList(1:blockCount) )

  if (hy_fullRiemannStateArrays) then
     call hy_memAllocScratch(SCRATCH_FACEX,&
       min(HY_P01_FACEXPTR_VAR,HY_N01_FACEXPTR_VAR),&
       2*HY_SCRATCH_NUM, &
       0,1,0, &
       blockList(1:blockCount) )
#if NDIM > 1
     call hy_memAllocScratch(SCRATCH_FACEY,&
       min(HY_P01_FACEYPTR_VAR,HY_N01_FACEYPTR_VAR),&
       2*HY_SCRATCH_NUM,&
       0,1,0, &
       blockList(1:blockCount) )
#endif
#if NDIM > 2
     call hy_memAllocScratch(SCRATCH_FACEZ,&
       min(HY_P01_FACEZPTR_VAR,HY_N01_FACEZPTR_VAR),&
       2*HY_SCRATCH_NUM,&
       0,1,0, &
       blockList(1:blockCount) )
#endif
  end if

#endif

end subroutine Hydro_prepareBuffers

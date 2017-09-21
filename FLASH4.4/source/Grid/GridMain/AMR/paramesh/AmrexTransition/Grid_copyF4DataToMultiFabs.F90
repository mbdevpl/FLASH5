subroutine Grid_copyF4DataToMultiFabs(gds, phi, nodetype, reverse)
  use Grid_interface,      ONLY : Grid_getLocalNumBlks, &
                                  Grid_getListOfBlocks, &
                                  Grid_getBlkPtr,&
                                  Grid_releaseBlkPtr

  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_meshComm
  use tree,      ONLY : lrefine_max, lrefine

  use amrex_multifab_module
  use amrex_distromap_module
  use amrex_boxarray_module
  use block_1lev_iterator, ONLY : block_1lev_iterator_t
  use block_metadata, ONLY : block_metadata_t

#include "Flash_mpi_implicitNone.fh"
#include "Flash.h"
#include "constants.h"

  type(amrex_multifab),OPTIONAL,intent(INOUT) :: phi(:)
  integer,intent(IN),OPTIONAL :: gds
  integer,intent(IN),OPTIONAL :: nodetype
  logical,intent(IN),OPTIONAL :: reverse

  type(block_1lev_iterator_t) :: itor
  type(block_metadata_t) :: blockDesc
  integer   :: gdsr
  integer   :: maxLev
  integer   :: myNodetype, listNodetype
  integer   :: localNumBlocks, globalNumBlocks
  logical   :: doReverse
  logical,save :: didWarn1 = .FALSE.
  real,POINTER,dimension(:,:,:,:) :: pf, pa

!!$  integer, dimension(LOW:HIGH,MDIM) :: tileLimits,blkLimitsGC, cornerID
  integer, dimension(LOW:HIGH,MDIM) :: cornerID
  integer, dimension(LOW:HIGH,NDIM) :: tileLimits
  integer, dimension(:,:,:),allocatable :: dimLimits
  integer, dimension(:,:,:),allocatable :: locLim
  integer, parameter :: limSize = 2*NDIM !number of integers in (LO:HI,:NDIM)
  integer, dimension(MDIM) :: stride
  integer, dimension(:), allocatable :: oldp,newp
  integer, dimension(:), allocatable :: procMap
  integer :: blockCount
  integer,dimension(:),allocatable :: blks
  integer,dimension(:),allocatable :: recvcount,displs
  integer :: offsLoc, offsGlob, levelBlocks, n

  integer :: ierr

  integer:: ib, blockID, level,ilev, iproc
  integer:: ibLoc, blkLev
  integer,allocatable :: lbpl( :) ! array for local blocks-per-level
  integer,allocatable :: gbpl( :) ! array for global blocks-per-level
  integer,allocatable :: disp( :) ! array for displacement
  integer,allocatable :: bpl(:,:) ! array for blocks-per-level (and per proc)

  if (present(gds)) then
     gdsr = gds
  else
     gdsr = CENTER
  end if
  if (present(nodetype)) then
     myNodetype = nodetype
  else
     myNodetype = LEAF
  end if
  if (myNodetype == ALL_BLKS) then
     listNodetype = REFINEMENT
  else
     listNodetype = myNodetype
  end if

  doReverse = .FALSE.
  if (present(reverse)) then
     doReverse = reverse
  end if

  if (present(phi)) then
     if (.NOT.didWarn1) &
          print*,'Grid_copyF4DataToMultiFabs: WARNING - argument phi is currently ignored.'
     didWarn1 = .TRUE.
  end if


  maxLev = lrefine_max

  call Grid_getLocalNumBlks(localNumBlocks)
  allocate(blks(localNumBlocks))

  do level = 1,maxLev
     call Grid_getListOfBlocks(listNodetype,blks,blockCount,level)
     if (blockCount > 0) then
        ib = 0
!!$        999 format('Grid_copyF4DataToMultiFabs[',A1,']: create iter for level',A3,':')
!!$        if (.NOT.doReverse) then
!!$           print 999,'>',level
!!$        else
!!$           print 999,'<',level
!!$        end if
        itor = block_1lev_iterator_t(myNodetype, level)
        do while (itor%is_valid())
           call itor%blkMetaData(blockDesc)
           ib = ib + 1
           do while(lrefine(blks(ib)).NE.level)
              ib = ib + 1
           end do
!!$           print*,ib,gr_meshMe,blks(ib),blockDesc%level,blockDesc%grid_index
           call Grid_getBlkPtr(blks(ib),pf,gds)
           call Grid_getBlkPtr(blockDesc,pa,gds)
           if (ANY(SHAPE(pf).NE.SHAPE(pa))) then
              print*,'lbound(pf),ubound(pf):',lbound(pf),ubound(pf)
              print*,'lbound(pa),ubound(pa):',lbound(pa),ubound(pa)
              call Driver_abortFlash("SHAPE SHIFTERS!")
           end if
           if (doReverse) then
              pf = pa           !Do the actual copying! (Flash4 data <- multifab array
           else
              pa = pf           !Do the actual copying! (Flash4 data -> multifab array
           end if
           call Grid_releaseBlkPtr(blks(ib),pf,gds)
           call Grid_releaseBlkPtr(blockDesc,pa,gds)
           call itor%next()
        end do
!!$        print*,'Grid_copyF4DataToMultiFabs: destroying iter for level',level,'...'
        call itor%destroy_iterator()
     end if
  end do

  deallocate(blks)

end subroutine Grid_copyF4DataToMultiFabs

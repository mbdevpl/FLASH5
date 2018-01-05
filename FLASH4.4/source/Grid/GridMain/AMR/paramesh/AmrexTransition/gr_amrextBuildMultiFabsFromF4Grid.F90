subroutine gr_amrextBuildMultiFabsFromF4Grid(gds, maxLev, nodetype)
  use Grid_interface,      ONLY : Grid_getLocalNumBlks, &
                                  Grid_getListOfBlocks, &
                                  Grid_getBlkIndexLimits, &
                                  Grid_getBlkCornerID, &
                                  Grid_updateRefinement,&
                                  Grid_fillGuardCells,&
                                  Grid_getDeltas,&
                                  Grid_getBlkPtr,&
                                  Grid_releaseBlkPtr,&
                                  Grid_getMaxRefinement

  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_meshComm
  use gr_physicalMultifabs, ONLY : Unk

  use amrex_multifab_module
  use amrex_distromap_module
  use amrex_boxarray_module

#include "Flash_mpi_implicitNone.fh"
#include "Flash.h"
#include "constants.h"

  integer,intent(IN) :: gds
!!$  type(amrex_multifab),intent(INOUT) :: phi_mf(0:)
  type(amrex_multifab),POINTER :: phi_mf(:)
  integer,intent(IN) :: maxLev
  integer,intent(IN),OPTIONAL :: nodetype

  integer   :: myNodetype
  integer   :: localNumBlocks, globalNumBlocks

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

  type(amrex_distromap) :: dm
  type(amrex_boxarray)  :: ba

  if (gds==CENTER) then
     phi_mf => Unk
  else                          ! for now
     call Driver_abortFlash('gr_amrextBuildMultiFabsFromF4Grid: gds must be CENTER!')
  end if

  if (present(nodetype)) then
     myNodetype = nodetype
  else
     myNodetype = LEAF
  end if

  call Grid_getLocalNumBlks(localNumBlocks)
  allocate(blks(localNumBlocks))

  allocate(lbpl(                   maxLev)) !local blocks per level
  allocate(gbpl(                   maxLev)) !global blocks per level
!!$  allocate(bpl(0:gr_meshNumProcs-1,maxLev)) !(global) blocks per level (and proc)
  allocate(bpl(maxLev,gr_meshNumProcs))
  lbpl(:) = 0
  bpl(:,:) = 0

  call Grid_getListOfBlocks(myNodetype,blks,blockCount)
  ibLoc = 0
  ! count how many blocks we have at each refinement level
  do ib=1,blockCount
     blockID=blks(ib)
     call Grid_getBlkRefineLevel(blockID,blkLev)
     ibLoc = lbpl(blkLev) + 1
     lbpl(blkLev) = ibLoc

#ifdef DEBUG_GRID
     print*,'ib,blockID,blkLev,ibLoc:',ib,blockID,blkLev,ibLoc
#endif
  end do

  allocate(disp(maxLev))
  if (maxLev>0) disp(1) = 0
  do ilev = 2,maxLev
     disp(ilev) = disp(ilev-1) + lbpl(ilev-1)
  end do

  allocate(locLim(LOW:HIGH,NDIM,blockCount))
  allocate(oldp(blockCount))
  allocate(newp(blockCount))
  lbpl(:) = 0
  tileLimits(:,:) = 1
  do ib=1,blockCount
     blockID=blks(ib)
     call Grid_getBlkRefineLevel(blockID,blkLev)
     ibLoc = lbpl(blkLev) + 1
     lbpl(blkLev) = ibLoc
     newp(ib) = disp(blkLev) + ibLoc !new position of block ib
#ifdef DEBUG_GRID
     print*,'ib,blockID,blkLev,ibLoc,disp(blkLev):',ib,blockID,blkLev,ibLoc,disp(blkLev)
#endif
     oldp(newp(ib)) = ib             !old position (in blks) of newp(ib)
     call Grid_getBlkCornerID(blockID,cornerID(LOW,:),stride,cornerID(HIGH,:))
     tileLimits(LOW ,:NDIM) = (cornerID(LOW ,:NDIM)-1) / stride(:NDIM) + 1
     tileLimits(HIGH,:NDIM) =  cornerID(HIGH,:NDIM) / stride(:NDIM)
     locLim(LOW:HIGH,1:NDIM,newp(ib)) = tileLimits(LOW:HIGH,1:NDIM)
  end do
  deallocate(disp)



  ! spread the information everywhere
  call MPI_Allgather(lbpl, maxLev, FLASH_INTEGER, bpl, maxLev, FLASH_INTEGER, &
       gr_meshComm, ierr)

  gbpl(:) = SUM(bpl,dim=2)

  globalNumBlocks = sum(gbpl)
  allocate(dimLimits(LOW:HIGH,NDIM,globalNumBlocks))
#ifdef DEBUG_GRID
  print*,'SHAPE(dimLimits) is',SHAPE(dimLimits)
#endif

  allocate(recvcount(0:gr_meshNumProcs))
  allocate(displs   (0:gr_meshNumProcs))
  recvcount(:) = 0; displs(:) = 0
  do iproc=1,gr_meshNumProcs
     displs(iproc) = displs(iproc-1) + recvcount(iproc-1)
     do ilev=1,maxLev
        recvcount(iproc) = recvcount(iproc) + bpl(ilev,iproc)
     end do
  end do

  recvcount(:) = limSize * recvcount(:)
  displs(  1:) = limSize * displs(  1:)
!!$  recvcount
!!$  displs

  call MPI_AllgatherV(locLim,    recvcount(gr_meshMe+1), FLASH_INTEGER, &
                      dimLimits, recvcount(1:gr_meshNumProcs), displs(1:gr_meshNumProcs), FLASH_INTEGER, &
                      gr_meshComm, ierr)

  displs(  1:) = displs(  1:) / limSize

  do ib=1,blockCount
     

  end do


  deallocate(newp)
  deallocate(oldp)
  deallocate(locLim)



  do level=1,maxLev
#ifdef DEBUG_GRID
     print*,' ***************   Amrext-BUILD LEVEL', level,'  **********************'
#endif
     levelBlocks = gbpl(level)
     allocate(procMap(levelBlocks))
     allocate(locLim(LOW:HIGH,NDIM,levelBlocks))
     offsLoc = 0
     do iproc=1,gr_meshNumProcs
        offsGlob = displs(iproc) + sum(bpl(1:level-1,iproc))
        n = bpl(level,iproc)
#ifdef DEBUG_GRID
        print*,'level,iproc,bpl(level,iproc):',level,iproc,bpl(level,iproc)
#endif
        if (n>0) then
           locLim(LOW:HIGH,:NDIM,offsLoc+1:offsLoc+n) = dimLimits(:,:NDIM,offsGlob+1:offsGlob+n)
           procMap              (offsLoc+1:offsLoc+n) = iproc - 1
           offsLoc = offsLoc + n
        end if
     end do

     call amrex_distromap_build(dm,procMap(1:levelBlocks))
#ifdef DEBUG_GRID
     call amrex_print(dm)
#endif
     call amrex_boxarray_build(ba,locLim(LOW:HIGH,:NDIM,1:levelBlocks))
#ifdef DEBUG_GRID
     call amrex_print(ba)
#endif
     call amrex_multifab_destroy(phi_mf(level-1))
     call amrex_multifab_build(phi_mf(level-1), ba, dm, NUNK_VARS, ng=NGUARD)

     call amrex_boxarray_destroy(ba)
     call amrex_distromap_destroy(dm)

     deallocate(locLim)
     deallocate(procMap)

  end do

  deallocate(bpl)
  deallocate(dimLimits)
  deallocate(blks)

end subroutine gr_amrextBuildMultiFabsFromF4Grid

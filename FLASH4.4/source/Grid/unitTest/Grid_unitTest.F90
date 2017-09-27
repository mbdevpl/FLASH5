!!****if* source/Grid/unitTest/Grid_unitTest
!!
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  call Grid_unitTest(integer(in):: fileUnit,
!!                     logical(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This unit test exercises the data accessing functions of the Grid unit
!!  The routine has direct access to all the mesh data structures such as 
!!  "unk", "facevarx" etc. It uses the Grid_getBlk/Point/RowData functions 
!!  to fetch some or all of the block data, and then compares it with
!!  the corresponding section of the appropriate array.
!!
!! ARGUMENTS
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!!
!!***

#undef FIXEDBLOCKSIZE

subroutine Grid_unitTest(fileUnit,perfect)


                        
  use physicaldata, ONLY : unk
  use Grid_interface, ONLY : Grid_getBlkData, &
       Grid_getPointData, Grid_getRowData, Grid_getMaxRefinement
  use block_iterator, ONLY : block_iterator_t
  use block_metadata, ONLY : block_metadata_t
  
  implicit none
  
#include "Flash.h"
#include "constants.h"

  
  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors


  real, allocatable, dimension(:,:,:,:) :: dataBlock
  real, pointer :: solnData(:,:,:,:)   

  integer,dimension(MDIM) :: size,startingPos

  integer, dimension(MDIM) :: errorCoords
  integer,dimension(NUNK_VARS) :: vars

  integer,dimension(MDIM)::sizeBlk
  integer,dimension(MDIM-1)::sizePlane
  integer :: sizeRow
  real, allocatable,dimension(:,:,:) :: dataBlk
  real, allocatable,dimension(:,:) :: dataPlane
  real, allocatable :: dataRow(:)
  real :: dataPoint

  integer :: lev, maxLev
  type(block_iterator_t) :: itor
  type(block_metadata_t) :: block
  integer :: ilocal,jlocal,klocal
  integer, dimension(LOW:HIGH,MDIM) :: limGC, lim,blkLimits,blkLimitsGC
  
  integer :: i,j,beg(MDIM),fin(MDIM)
  real :: error

  ! DEVNOTE: Update code to change out paramesh-specific code
  error = 0.0
  call Grid_getMaxRefinement(maxLev,mode=1)

  do lev=1,maxLev
     itor = block_iterator_t(LEAF, level=lev)
     do while(itor%is_valid())
        call itor%blkMetaData(block)
        limGC=block%limitsGC
        lim=block%limits
        blkLimits=block%localLimits
        blkLimitsGC=block%localLimitsGC
        
        startingPos=1
        beg=startingPos
        fin=beg+blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)
        sizeBlk=fin-beg+1

        allocate(dataBlk(sizeBlk(IAXIS), sizeBlk(JAXIS), sizeBlk(KAXIS)))     

        call Grid_getBlkData(block, CENTER, DENS_VAR, EXTERIOR, startingPos, &
             dataBlk, sizeBlk)
        error = error+maxval(abs(dataBlk(:,:,:)-&
             unk(DENS_VAR,beg(IAXIS):fin(IAXIS),beg(JAXIS):fin(JAXIS),beg(KAXIS):fin(KAXIS),block%id)))

        startingPos=1
        beg=startingPos+blkLimits(LOW,:)-1
        fin=beg+blkLimits(HIGH,:)-blkLimits(LOW,:)
        sizeBlk=fin-beg+1

        deallocate(dataBlk)
        allocate(dataBlk(sizeBlk(IAXIS), sizeBlk(JAXIS), sizeBlk(KAXIS)))     
             
        call Grid_getBlkData(block, CENTER, DENS_VAR, INTERIOR, startingPos, &
             dataBlk, sizeBlk)
        error = error+maxval(abs(dataBlk(:,:,:)-&
             unk(DENS_VAR,beg(IAXIS):fin(IAXIS),beg(JAXIS):fin(JAXIS),beg(KAXIS):fin(KAXIS),block%id)))
        deallocate(dataBlk)
        print*,'the cumulative error after exterior blk get is ',error


        !!  Now test Grid_getPointData, 1 var, EXTERIOR
        startingPos(1) = 5
        startingPos(2) = 6
        startingPos(3) = 1

        beg=startingPos
        fin=beg
        call Grid_getPointData(block, CENTER, DENS_VAR, EXTERIOR, startingPos, &
             dataPoint)

        error = error + abs(dataPoint-unk(DENS_VAR,beg(IAXIS),beg(JAXIS),beg(KAXIS),block%id))
        print*,'the cumulative error after exterior point get is ',error

        startingPos(1) = 7
        startingPos(2) = 12
        startingPos(3) = 1
        print *, "testing Grid_getPointData, 1 var, EXTERIOR, from guardcells"
        
        beg=startingPos
        fin=beg

        call Grid_getPointData(block, CENTER, DENS_VAR, EXTERIOR, startingPos, &
             dataPoint)
        
        error = error + abs(dataPoint-unk(DENS_VAR,beg(IAXIS),beg(JAXIS),beg(KAXIS),block%id))

             


        !  Now test Grid_getPointData, 1 var, INTERIOR
        print *, "testing Grid_getPointData, 1 var, INTERIOR"
        
        startingPos(1) = 7
        startingPos(2) = 3
        startingPos(3) = 1
        beg=startingPos+blkLimits(LOW,:)-1
        fin=beg
        
        call Grid_getPointData(block, CENTER, DENS_VAR, INTERIOR, startingPos, &
          dataPoint)
        error = error + abs(dataPoint-unk(DENS_VAR,beg(IAXIS),beg(JAXIS),beg(KAXIS),block%id))

        print*,'the cumulative error is ',error
        call itor%next()     
     end do
  end do
  print*,'the cumulative error is ',error
  print *, dataPoint
  
!!$
!!$
!!$
!!$
!!$
!!$  Now test Grid_getRowData, 1 var, INTERIOR
!!$  print *, "testing Grid_getRowData, 1 var, INTERIOR, I AXIS"
!!$
!!$  startingPos(1) = 3
!!$  startingPos(2) = 1
!!$  startingPos(3) = 4
!!$
!!$  dataSizeRow = gr_ihi-gr_ilo+1-startingPos(1)+1
!!$
!!$  ib=startingPos(1)+grd(1)
!!$  ie=ib+dataSizeRow-1
!!$  jb=startingPos(2)+grd(2)
!!$  je=jb
!!$  kb=startingPos(3)+grd(3)
!!$  ke=kb
!!$
!!$
!!$
!!$
!!$  allocate(datablockRow(dataSizeRow))
!!$    
!!$  do i = 1,localNumblks
!!$     call Grid_getRowData(i, CENTER, DENS_VAR, INTERIOR, IAXIS, startingPos, &
!!$          datablockRow, dataSizeRow)
!!$
!!$     error = error+maxval(abs(datablockRow-unk(DENS_VAR,ib:ie,jb,kb,i)))
!!$     
!!$  end do
!!$  print*,'error is ',error
!!$ 
!!$  deallocate(datablockRow)
!!$
!!$ 
!!$
!!$
!!$  Now test Grid_getRowData, 1 var, INTERIOR JAXIS
!!$  print *, "testing Grid_getRowData, 1 var, INTERIOR, JAXIS"
!!$
!!$  startingPos(1) = 1
!!$  startingPos(2) = 5
!!$  startingPos(3) = 1
!!$
!!$  dataSizeRow = gr_jhi-gr_jlo+1-startingPos(2)+1
!!$
!!$  ib=startingPos(1)+grd(1)
!!$  ie=ib
!!$  jb=startingPos(2)+grd(2)
!!$  je=jb+dataSizeRow-1
!!$  kb=startingPos(3)+grd(3)
!!$  ke=kb
!!$
!!$  allocate(datablockRow(dataSizeRow))
!!$    
!!$
!!$  do i = 1,localNumblks
!!$     call Grid_getRowData(i, CENTER, DENS_VAR, INTERIOR, JAXIS, startingPos, &
!!$           datablockRow, dataSizeRow)
!!$     
!!$     error = error+maxval(abs(datablockRow-unk(DENS_VAR,ib,jb:je,kb,i)))
!!$  end do
!!$  print*,'error is ',error
!!$  
!!$  deallocate(datablockRow)
!!$  
!!$
!!$
!!$
!!$
!!$  startingPos(1) = 2
!!$  startingPos(2) = 4
!!$  startingPos(3) = 6
!!$
!!$  dataSizeRow = gr_ihiGC-startingPos(1)+1
!!$
!!$  ib=startingPos(1)
!!$  ie=ib+dataSizeRow-1
!!$  jb=startingPos(2)
!!$  je=jb
!!$  kb=startingPos(3)
!!$  ke=kb
!!$
!!$  allocate(datablockRow(dataSizeRow))
!!$    
!!$  print *, "testing Grid_getRowData, 1 var, EXTERIOR, IAXIS"
!!$  do i = 1,localNumblks
!!$     call Grid_getRowData(i, CENTER, DENS_VAR, EXTERIOR, IAXIS, &
!!$          startingPos, datablockRow, dataSizeRow)
!!$     error = error+maxval(abs(datablockRow-unk(DENS_VAR,ib:ie,jb,kb,i)))
!!$  end do
!!$  print*,'error is ',error
!!$ 
!!$  deallocate(datablockRow)
  
  return
 
end subroutine Grid_unitTest

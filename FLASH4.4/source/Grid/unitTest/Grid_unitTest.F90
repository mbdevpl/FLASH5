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


                        
  use physicaldata, ONLY : unk,facevarx,facevary,facevarz
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
  real :: dataPoint, dataPoint2
  real,dimension(MDIM) :: point
  integer :: lev, maxLev
  type(block_iterator_t) :: itor
  type(block_metadata_t) :: block
  integer :: ilocal,jlocal,klocal
  integer, dimension(LOW:HIGH,MDIM) :: limGC, lim,blkLimits,blkLimitsGC
  
  integer :: i,j,beg(MDIM),fin(MDIM),iSize,jSize,kSize
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
        
        !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Test for Grid_getBlkData  !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        startingPos=1
        beg=startingPos+blkLimits(LOW,:)-1
        fin=beg+blkLimits(HIGH,:)-blkLimits(LOW,:)
        fin(IAXIS)=fin(IAXIS)+1

        iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
        jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
        kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1


        allocate(dataBlk(lim(LOW,IAXIS):lim(HIGH,IAXIS)+1,&
             lim(LOW,JAXIS):lim(HIGH,JAXIS), &
             lim(LOW,KAXIS):lim(HIGH,KAXIS)))
!        call Grid_getBlkData(block, FACEX, AREA_FACE_VAR, GLOBALIDX1, &
        call Grid_getBlkData(block, CELL_FACEAREA, ILO_FACE,GLOBALIDX1, &
             (/lim(LOW,IAXIS),lim(LOW,JAXIS),lim(LOW,KAXIS)/), &
             dataBlk(lim(LOW,IAXIS):lim(HIGH,IAXIS)+1,&
             lim(LOW,JAXIS):lim(HIGH,JAXIS),  &
             lim(LOW,KAXIS):lim(HIGH,KAXIS)), &
             (/isize+1, jsize, ksize/) )
        error = error+maxval(abs(dataBlk(:,:,:)-&
             facevarx(1,beg(IAXIS):fin(IAXIS),beg(JAXIS):fin(JAXIS),beg(KAXIS):fin(KAXIS),block%id)))

        startingPos=1
        beg=startingPos+blkLimits(LOW,:)-1
        fin=beg+blkLimits(HIGH,:)-blkLimits(LOW,:)

        deallocate(dataBlk)
        print*,'Grid_getBlkData  ::  error after face get is ',error
        if (error > 0) perfect = .FALSE.
             
        error=0.0
        allocate(dataBlk(lim(LOW,IAXIS):lim(HIGH,IAXIS),&
             lim(LOW,JAXIS):lim(HIGH,JAXIS), &
             lim(LOW,KAXIS):lim(HIGH,KAXIS)))
        call Grid_getBlkData(block, CELL_VOLUME, 0, GLOBALIDX1, &
             (/lim(LOW,IAXIS),lim(LOW,JAXIS),lim(LOW,KAXIS)/), &
            dataBlk(lim(LOW,IAXIS):lim(HIGH,IAXIS), &
            lim(LOW,JAXIS):lim(HIGH,JAXIS), &
            lim(LOW,KAXIS):lim(HIGH,KAXIS)), &
            (/isize, jsize, ksize/) )
!!$        call Grid_getBlkData(block, FACEY, 1, INTERIOR, startingPos, &
!!$             dataBlk, sizeBlk)
        error = error+maxval(abs(dataBlk(:,:,:)-&
             unk(1,beg(IAXIS):fin(IAXIS),beg(JAXIS):fin(JAXIS),beg(KAXIS):fin(KAXIS),block%id)))
        deallocate(dataBlk)
        print*,'Grid_getBlkData  ::  the cumulative error after Center get is ',error
        if (error > 0) perfect = .FALSE.

        !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -------x------- !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Test for Grid_getPlaneData  !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        startingPos=1
        beg=startingPos+blkLimits(LOW,:)-1
        fin=beg+blkLimits(HIGH,:)-blkLimits(LOW,:)
        fin(IAXIS)=fin(IAXIS)+1

        iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
        jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
        kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

        error=0.0
        allocate(dataPlane(lim(LOW,IAXIS):lim(HIGH,IAXIS)+1,&
             lim(LOW,JAXIS):lim(HIGH,JAXIS)))
        call Grid_getPlaneData(block, CELL_FACEAREA, ILO_FACE,GLOBALIDX1, XYPLANE, &
             (/lim(LOW,IAXIS),lim(LOW,JAXIS),lim(LOW,KAXIS)/), &
             dataPlane(lim(LOW,IAXIS):lim(HIGH,IAXIS)+1,&
             lim(LOW,JAXIS):lim(HIGH,JAXIS)), &
             (/isize+1, jsize/) )

        error = error+maxval(abs(dataPlane(:,:)-&
             facevarx(1,beg(IAXIS):fin(IAXIS),beg(JAXIS):fin(JAXIS),beg(KAXIS),block%id)))
        
        startingPos=1
        beg=startingPos+blkLimits(LOW,:)-1
        fin=beg+blkLimits(HIGH,:)-blkLimits(LOW,:)

        deallocate(dataPlane)
        print*,'Grid_getPlaneData  ::  error after face get is ',error
        if (error > 0) perfect = .FALSE.
             
        error=0.0
        allocate(dataPlane(lim(LOW,IAXIS):lim(HIGH,IAXIS),&
             lim(LOW,JAXIS):lim(HIGH,JAXIS)))
        call Grid_getPlaneData(block, CELL_VOLUME, 0, GLOBALIDX1, XYPLANE, &
             (/lim(LOW,IAXIS),lim(LOW,JAXIS),lim(LOW,KAXIS)/), &
            dataPlane(lim(LOW,IAXIS):lim(HIGH,IAXIS), &
            lim(LOW,JAXIS):lim(HIGH,JAXIS)), &
            (/isize, jsize/) )
        error = error+maxval(abs(dataPlane(:,:)-&
             unk(1,beg(IAXIS):fin(IAXIS),beg(JAXIS):fin(JAXIS),beg(KAXIS),block%id)))
        deallocate(dataPlane)
        print*,'Grid_getPlaneData  ::  the cumulative error after Center get is ',error
        if (error > 0) perfect = .FALSE.

        !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -------x------- !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Test for Grid_getRowData  !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        startingPos=1
        beg=startingPos+blkLimits(LOW,:)-1
        fin=beg+blkLimits(HIGH,:)-blkLimits(LOW,:)
        fin(IAXIS)=fin(IAXIS)+1

        iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
!         jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
!         kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

        error=0.0
        allocate(dataRow(lim(LOW,IAXIS):lim(HIGH,IAXIS)+1))
        call Grid_getRowData(block, CELL_FACEAREA, ILO_FACE,GLOBALIDX1, IAXIS, &
             (/lim(LOW,IAXIS),lim(LOW,JAXIS),lim(LOW,KAXIS)/), &
             dataRow(lim(LOW,IAXIS):lim(HIGH,IAXIS)+1), &
             (isize+1) )
        error = error+maxval(abs(dataRow(:)-&
             facevarx(1,beg(IAXIS):fin(IAXIS),beg(JAXIS),beg(KAXIS),block%id)))
        
        startingPos=1
        beg=startingPos+blkLimits(LOW,:)-1
        fin=beg+blkLimits(HIGH,:)-blkLimits(LOW,:)

        deallocate(dataRow)
        print*,'Grid_getRowData  ::  error after face get is ',error
        if (error > 0) perfect = .FALSE.
             
        error=0.0
        allocate(dataRow(lim(LOW,IAXIS):lim(HIGH,IAXIS)))
        call Grid_getRowData(block, CELL_VOLUME, 0, GLOBALIDX1, IAXIS, &
             (/lim(LOW,IAXIS),lim(LOW,JAXIS),lim(LOW,KAXIS)/), &
             dataRow(lim(LOW,IAXIS):lim(HIGH,IAXIS)), &
             isize )
        error = error+maxval(abs(dataRow(:)-&
             unk(1,beg(IAXIS):fin(IAXIS),beg(JAXIS),beg(KAXIS),block%id)))
        deallocate(dataRow)
        print*,'Grid_getRowData  ::  the cumulative error after Center get is ',error
        if (error > 0) perfect = .FALSE.

        !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Test for Grid_getPointData  !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        startingPos=1
        beg=startingPos+blkLimits(LOW,:)-1
        fin=beg+blkLimits(HIGH,:)-blkLimits(LOW,:)
        fin(IAXIS)=fin(IAXIS)+1

        error=0.0
!             Grid_getRowData(blockID, gridDataStruct, variable, beginCount, &
!                     row, startingPos, datablock, dataSize)
!             Grid_getPointData(block, gridDataStruct, variable, beginCount, &
!                     position, datablock)
        call Grid_getPointData(block, CELL_FACEAREA, ILO_FACE,GLOBALIDX1, &
             (/lim(LOW,IAXIS),lim(LOW,JAXIS),lim(LOW,KAXIS)/), &
             dataPoint )
        error = error+abs(dataPoint-&
             facevarx(1,beg(IAXIS),beg(JAXIS),beg(KAXIS),block%id))
        
        startingPos=1
        beg=startingPos+blkLimits(LOW,:)-1
        fin=beg+blkLimits(HIGH,:)-blkLimits(LOW,:)

        print*,'Grid_getPointData  ::  error after face get is ',error
        if (error > 0) perfect = .FALSE.
             
        error=0.0
        call Grid_getPointData(block, CELL_VOLUME, 0, GLOBALIDX1, &
             (/lim(LOW,IAXIS),lim(LOW,JAXIS),lim(LOW,KAXIS)/), &
             dataPoint )
        error = error+abs(dataPoint-&
             unk(1,beg(IAXIS),beg(JAXIS),beg(KAXIS),block%id))
        print*,'Grid_getPointData  ::  the cumulative error after Center get is ',error
        if (error > 0) perfect = .FALSE.

        !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -------x------- !!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$        !!  Now test Grid_getPointData, 1 var, EXTERIOR
!!$        startingPos(1) = 5
!!$        startingPos(2) = 6
!!$        startingPos(3) = 1
!!$
!!$        beg=startingPos
!!$        fin=beg
!!$        call Grid_getPointData(block, CENTER, DENS_VAR, EXTERIOR, startingPos, &
!!$             dataPoint)
!!$
!!$        error = error + abs(dataPoint-unk(DENS_VAR,beg(IAXIS),beg(JAXIS),beg(KAXIS),block%id))
!!$        print*,'the cumulative error after exterior point get is ',error
!!$
!!$        startingPos(1) = 7
!!$        startingPos(2) = 12
!!$        startingPos(3) = 1
!!$        print *, "testing Grid_getPointData, 1 var, EXTERIOR, from guardcells"
!!$        
!!$        beg=startingPos
!!$        fin=beg
!!$
!!$        call Grid_getPointData(block, CENTER, DENS_VAR, EXTERIOR, startingPos, &
!!$             dataPoint)
!!$        
!!$        error = error + abs(dataPoint-unk(DENS_VAR,beg(IAXIS),beg(JAXIS),beg(KAXIS),block%id))
!!$
!!$             
!!$
!!$
!!$        !  Now test Grid_getPointData, 1 var, INTERIOR
!!$        print *, "testing Grid_getPointData, 1 var, INTERIOR"
!!$        
!!$        startingPos(1) = 7
!!$        startingPos(2) = 3
!!$        startingPos(3) = 1
!!$        beg=startingPos+blkLimits(LOW,:)-1
!!$        fin=beg
!!$        
!!$        call Grid_getPointData(block, CENTER, DENS_VAR, INTERIOR, startingPos, &
!!$          dataPoint)
!!$        error = error + abs(dataPoint-unk(DENS_VAR,beg(IAXIS),beg(JAXIS),beg(KAXIS),block%id))
!!$
!!$        print*,'the cumulative error is ',error
        call itor%next()     
     end do
  end do
!   print*,'the cumulative error is ',error
  
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

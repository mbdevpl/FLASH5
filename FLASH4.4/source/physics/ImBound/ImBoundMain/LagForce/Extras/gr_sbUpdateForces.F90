!!****if* source/physics/ImBound/ImBoundMain/LagForce/Extras/gr_sbUpdateForces
!!
!! NAME
!!  gr_sbUpdateForces
!!
!! SYNOPSIS
!!
!!  gr_sbUpdateForces()
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!***

#define DO_MARKER_FORCING /* If commented, no IB forcing from Markers */
!#define  INVERSE_GCELL_FILL 
!#define ITER_FORCING /* Iterative procedure in Kermpe and Frohlich 2005, JCP */

#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

subroutine gr_sbUpdateForces
  
  use Grid_interface, ONLY : Grid_updateSolidBodyForces, Grid_getBlkPtr,     &
                             Grid_releaseBlkPtr, Grid_getListOfBlocks,       &
                             Grid_getDeltas, Grid_getBlkCenterCoords,        &
                             Grid_getBlkPhysicalSize,Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells

  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, totalPart, &
                        gr_sbDebug, gr_sbParticleCount, solid_body

  use gr_ptData, ONLY :  gr_ptBlkList, gr_ptBlkCount

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use ImBound_data, only : ib_dt, ib_maxIterForcing

  use gr_ptInterface, ONLY : gr_ptMove, gr_ptSetIndices, gr_ptResetIndices

  use Particles_data, ONLY : pt_indexList, pt_indexCount,pt_maxPerProc,pt_posinitialized

  use ib_interface, only : ib_forceInsideBody

  use Driver_interface, only : Driver_getNStep

  use ins_interface, only : ins_fluxfix

#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs, gr_maxParticlesPerProc
  use tree, only : neigh, lrefine, lrefine_max
#else
  use Grid_data, ONLY : gr_axisComm, gr_exch, gr_gridDataStruct, &
       gr_justExchangedGC,gr_domainBC, &
       gr_offset, gr_meshMe, gr_meshComm, gr_meshNumProcs
#endif


  implicit none
#include "Flash_mpi.h"
  type(solid_body), pointer :: bodyInfo
  real, dimension(MDIM) :: particleposn
  integer :: i,j,k,b,p, gettingFrom, blockID, recvCount
  real :: particleData(NPART_PROPS)
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData,facexData2,faceyData2,facezData2
  integer, save, dimension(MAXBLOCKS) :: listOfBlocks
  integer, save :: count
  integer :: lb,ierr
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: ct,localCount,globalCount,localPart
  logical :: moveDone,coords_in_blks=.true.

  real :: del(MDIM),coord(MDIM),bsize(MDIM)

  integer :: iForceIter, maxForceIter
  character(len=2) :: indForceIter
  character(len=6) :: indNStep
  integer :: NStep

  !! Forcing Test Variables
  real :: Body_Cen(MDIM),Voli
  real :: MomArmi_x, MomArmi_y, MomArmi_z
  real :: Fxtot(gr_sbNumBodies),  Fytot(gr_sbNumBodies), &
          Fxtoti(gr_sbNumBodies),Fytoti(gr_sbNumBodies)
  real :: Momz(gr_sbNumBodies),Momzi(gr_sbNumBodies)                                
  real :: Momx(gr_sbNumBodies),  Momxi(gr_sbNumBodies),  &
          Momy(gr_sbNumBodies),  Momyi(gr_sbNumBodies)
  real :: Fztot(gr_sbNumBodies), Fztoti(gr_sbNumBodies) 

  real :: dx,dy,dz,xcell,ycell,zcell,dxdydz
#ifdef FLASH_GRID_PARAMESH
  integer :: nxc,nyc,nzc
#endif
 	 
  !! Inverse gcfill variables:
#ifndef INVERSE_GCELL_FILL 

#ifdef FLASH_GRID_PARAMESH
  real :: particles2(NPART_PROPS,gr_maxParticlesPerProc)
#else
  real :: particles2(NPART_PROPS,pt_maxPerProc)
#endif

  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)

#else
  integer, dimension(LOW:HIGH,LOW:HIGH) :: istr_A,istr_F,jstr_A,jstr_F,kstr_A,kstr_F
  integer, dimension(LOW:HIGH,LOW:HIGH) :: iend_A,iend_F,jend_A,jend_F,kend_A,kend_F
  integer, dimension(MDIM) :: shf_low, myflg
!  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
  integer :: myaxis,beginDataType,endDataType
  integer, dimension(MDIM, MDIM+1) :: sendRight, sendLeft, recvRight, recvLeft

  real :: forcexData(GRID_IHI_GC+1,GRID_JHI_GC  ,GRID_KHI_GC  ,MAXBLOCKS)
  real :: forceyData(GRID_IHI_GC  ,GRID_JHI_GC+1,GRID_KHI_GC  ,MAXBLOCKS)
  real :: forcezData(GRID_IHI_GC  ,GRID_JHI_GC  ,GRID_KHI_GC+1,MAXBLOCKS)

  ! Guardcell info for Force x:
  real :: iaxis_forcexData(NGUARD       , GRID_JHI_GC  ,GRID_KHI_GC  ,LOW:HIGH, MAXBLOCKS)
  real :: jaxis_forcexData(GRID_IHI_GC+1, NGUARD       ,GRID_KHI_GC  ,LOW:HIGH, MAXBLOCKS)
  real :: kaxis_forcexData(GRID_IHI_GC+1, GRID_JHI_GC  ,NGUARD       ,LOW:HIGH, MAXBLOCKS)

  ! Guardcell info for Force y:
  real :: iaxis_forceyData(NGUARD       , GRID_JHI_GC+1,GRID_KHI_GC  ,LOW:HIGH, MAXBLOCKS)
  real :: jaxis_forceyData(GRID_IHI_GC  , NGUARD       ,GRID_KHI_GC  ,LOW:HIGH, MAXBLOCKS)
  real :: kaxis_forceyData(GRID_IHI_GC  , GRID_JHI_GC+1,NGUARD       ,LOW:HIGH, MAXBLOCKS)

  ! Guardcell info for Force z:
  real :: iaxis_forcezData(NGUARD       , GRID_JHI_GC  ,GRID_KHI_GC+1,LOW:HIGH, MAXBLOCKS)
  real :: jaxis_forcezData(GRID_IHI_GC  , NGUARD       ,GRID_KHI_GC+1,LOW:HIGH, MAXBLOCKS)
  real :: kaxis_forcezData(GRID_IHI_GC  , GRID_JHI_GC  ,NGUARD       ,LOW:HIGH, MAXBLOCKS)
#endif

  real :: FxTotL, FyTotL, FzTotL, MomxL, MomyL, MomzL, &
          Fx, Fy, Fz, FxProc, FyProc, FzProc, MomXProc, MomYProc, MomZProc


  integer, save :: countf = 0
  integer :: countRP, countVP
  character(len=28) :: fileName
  character(len=6) :: index_count,index_proc

  integer :: val
  integer :: Particles_Allbods


#ifdef FLASH_GRID_PARAMESH
  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1
#endif

  Particles_Allbods = 0

  !FzTotL = 0.

!  call Timers_start("update_forces")

  call Timers_start('InitEulerForce')

#ifdef FLASH_GRID_PARAMESH
  call Grid_getListOfBlocks(LEAF, listOfBlocks, count)
  call Grid_getListOfBlocks(LEAF, gr_ptBlkList, gr_ptBlkCount)
#else
  count=CONSTANT_ONE
  listOfBlocks(CONSTANT_ONE)=CONSTANT_ONE
#endif

  ! Set temp Eulerian Force Variable to zero:
#ifdef ITER_FORCING

!  maxForceIter = 3

  do lb=1,count

     blockID = listOfBlocks(lb)

     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     facexData(FORO_FACE_VAR,:,:,:) = 0.
     faceyData(FORO_FACE_VAR,:,:,:) = 0.

#if NDIM == MDIM
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     facezData(FORO_FACE_VAR,:,:,:) = 0.
#endif

     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo

  call MPI_Barrier (gr_meshComm, ierr) 

  do b = 1, gr_sbNumBodies
    bodyInfo => gr_sbBodyInfo(b)
    if(gr_sbParticleCount(b) .gt. 0) then    
    bodyInfo % particles(FULO_PART_PROP,:) = 0. 
    bodyInfo % particles(FVLO_PART_PROP,:) = 0. 
#if NDIM == MDIM
    bodyInfo % particles(FWLO_PART_PROP,:) = 0. 
#endif
    endif
  enddo
!#else
!  maxForceIter = 1
#endif


#ifdef ITER_FORCING
  do iForceIter = 1, ib_maxIterForcing ! iterative forcing
#endif

  Particles_Allbods = 0

  ! Set Eulerian Force Variable to zero:
  do lb=1,count

     blockID = listOfBlocks(lb)

     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     facexData(FORC_FACE_VAR,:,:,:) = 0.
     faceyData(FORC_FACE_VAR,:,:,:) = 0.

#if NDIM == MDIM
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     facezData(FORC_FACE_VAR,:,:,:) = 0.
#endif

     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

  enddo

  call Timers_stop('InitEulerForce') ! Shizhao 

#ifdef DO_MARKER_FORCING
  call Timers_start("Force_particles")

  ! Loop over all bodies
  do b = 1, gr_sbNumBodies
     bodyInfo => gr_sbBodyInfo(b)

#ifndef INVERSE_GCELL_FILL
     call Timers_start("Create_Particles2")
     particles2(:,:) = 0.
#ifdef FLASH_GRID_PARAMESH
     globalCount = gr_maxParticlesPerProc  
     moveDone=.false.
#else
     globalCount=pt_maxPerProc 
#endif

     !! Now get all the indices into the data structure setup right
     !! for the rest of the unit
     call gr_ptSetIndices(pt_indexList, pt_indexCount)

     if (bodyInfo % myPE == bodyInfo % bodyMaster) then

        totalPart = bodyInfo % totalPart

        ct = 0 
        do i = 1, totalPart
           if (int(bodyInfo % particles(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then
              ct= ct+1
              particles2(:,ct) = bodyInfo % particles(:,i)
           endif
        enddo
        
        localCount  = ct

        ! Generate virtual Particles: This call to gr_ptMove is calles with the
        ! only purpose of generating virtual particles for surface markers on the
        ! framework of the Immersed Boundary method. This might have to change if
        ! Different types of particles are used.
#ifdef FLASH_GRID_PARAMESH
        call gr_ptMove(particles2,NPART_PROPS, localCount,globalCount, moveDone)
#else
        call Timers_start('GridMoveParticles')
        pt_posinitialized=.true.
        call Grid_moveParticles(particles2,NPART_PROPS, globalCount, localCount, &
             pt_indexList, pt_indexCount, coords_in_blks)
        call Timers_stop('GridMoveParticles')
#endif
     else

        localCount = gr_sbParticleCount(b)
  
        if (localCount .gt. 0) &
        particles2(:,1:localCount) = bodyInfo%particles(:,1:localCount)

        ! Generate virtual Particles:
#ifdef FLASH_GRID_PARAMESH     
        call gr_ptMove(particles2,NPART_PROPS, localCount,globalCount, moveDone)
#else
        call Timers_start('GridMoveParticles')
        pt_posinitialized=.true.
        call Grid_moveParticles(particles2,NPART_PROPS, globalCount, localCount, &
             pt_indexList, pt_indexCount, coords_in_blks)
        call Timers_stop('GridMoveParticles')
#endif
     endif
     call Timers_stop("Create_Particles2")

     call gr_ptResetIndices(pt_indexList, pt_indexCount)     

#endif  /* #ifndef INVERSE_GCELL_FILL */

!!$     ! Do some tests on virtual particles
!!$     ! count real and virtual Particles
!!$     countRP = 0; countVP = 0;
!!$     write(*,*) gr_meshMe,'LocalCount after gr_ptMove=',localCount
!!$     do i = 1,localCount
!!$        if (int(particles2(TAG_PART_PROP,i)) .eq. 1) then
!!$           countRP = countRP + 1
!!$        elseif (int(particles2(TAG_PART_PROP,i)) .eq. -1) then       
!!$           countVP = countVP + 1
!!$        endif    
!!$     enddo
!!$     write(*,*) gr_meshMe,'Markers=',countRP,'Virtual Particles=',countVP
!!$
!!$
!!$     ! Write Real Particles:
!!$     countf=countf+1
!!$     write(index_count,"(I6.6)") countf
!!$     write(index_proc,"(I6.6)") gr_meshMe
!!$     open(unit=113,file='./IOData/Rpart'//index_count//'.'//index_proc//'.plt',form='formatted')
!!$
!!$     write(113,'(A)') 'VARIABLES = "X" , "Y" , "Z"'
!!$     write(113,'(A,I8,A)')'ZONE I=',countRP,',DATAPACKING = POINT'
!!$
!!$     do i = 1,localCount
!!$        if (int(particles2(TAG_PART_PROP,i)) .eq. 1) then
!!$           write(113,'(3E16.8)') particles2(POSX_PART_PROP,i),particles2(POSY_PART_PROP,i), &
!!$                                 particles2(POSZ_PART_PROP,i)     
!!$        end if
!!$     enddo
!!$     close(113)
!!$
!!$
!!$     ! Write Virtual Particles:
!!$     open(unit=113,file='./IOData/Vpart'//index_count//'.'//index_proc//'.plt',form='formatted')
!!$
!!$     write(113,'(A)') 'VARIABLES = "X" , "Y" , "Z", "blk"'
!!$     write(113,'(A,I8,A)')'ZONE I=',countVP,',DATAPACKING = POINT'
!!$
!!$     do i = 1,localCount
!!$        if (int(particles2(TAG_PART_PROP,i)) .eq. -1) then
!!$           write(113,'(4E16.8)') particles2(POSX_PART_PROP,i),particles2(POSY_PART_PROP,i), &
!!$                                 particles2(POSZ_PART_PROP,i),particles2(BLK_PART_PROP,i)
!!$        end if
!!$     enddo
!!$     close(113)



!!$     ! Write out for Matlab:
!!$     countf=countf+1
!!$     write(index_count,"(I6.6)") countf
!!$     write(index_proc,"(I6.6)") gr_meshMe
!!$     open(unit=113,file='./IOData/Rpart'//index_count//'.'//index_proc//'.dat',form='formatted')
!!$
!!$     do i = 1,localCount
!!$        if (int(particles2(TAG_PART_PROP,i)) .eq. 1) then
!!$           write(113,'(4E16.8)') particles2(POSX_PART_PROP,i),particles2(POSY_PART_PROP,i), &
!!$                                 particles2(NMLX_PART_PROP,i),particles2(NMLY_PART_PROP,i)
!!$        end if
!!$     enddo
!!$     close(113)
!!$
!!$     ! Write Virtual Particles:
!!$     open(unit=113,file='./IOData/Vpart'//index_count//'.'//index_proc//'.dat',form='formatted')
!!$
!!$     do i = 1,localCount
!!$        if (int(particles2(TAG_PART_PROP,i)) .eq. -1) then
!!$           write(113,'(4E16.8)') particles2(POSX_PART_PROP,i),particles2(POSY_PART_PROP,i), &
!!$                                 particles2(NMLX_PART_PROP,i),particles2(NMLY_PART_PROP,i)
!!$        end if
!!$     enddo
!!$     close(113)     
!!$
!!$
!!$     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
!!$        ! Write Lagrangian Mesh:
!!$        open(unit=113,file='./IOData/Mesh'//index_count//'.'//index_proc//'.dat',form='formatted')
!!$
!!$        do i = 1,bodyInfo%NumVertices
!!$           write(113,'(2E16.8)') bodyInfo%xb(i),bodyInfo%yb(i)
!!$        enddo
!!$        close(113)     
!!$     endif


     recvCount = 0
     if (bodyInfo % myPE == bodyInfo % bodyMaster) then

#ifndef INVERSE_GCELL_FILL

        ! Find number of particles in the Master Processor:
        localPart = 0
        do j=1,size(bodyInfo%particlesPerProc,DIM=2)
           if(bodyInfo%particlesPerProc(1,j) .eq. bodyInfo%bodyMaster) then
              localPart = bodyInfo%particlesPerProc(2,j)
           endif
        enddo

        particles2(FUL_PART_PROP,:) = 0.
        particles2(FVL_PART_PROP,:) = 0.
#if NDIM == MDIM
        particles2(FWL_PART_PROP,:) = 0.
#endif
        recvCount = localCount

        do i = 1, localCount

           if (int(particles2(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then
#ifdef FLASH_GRID_PARAMESH
              blockID = int(particles2(BLK_PART_PROP,i))
#else
              blockID=CONSTANT_ONE
#endif
              particleData = particles2(1:NPART_PROPS,i)


#else /* #ifndef INVERSE_GCELL_FILL */

        totalPart = bodyInfo % totalPart
        ct = 0
        do i = 1, totalPart
           if (int(bodyInfo % particles(PROC_PART_PROP,i)) == bodyInfo%bodyMaster) then
              ct= ct+1
           endif
        enddo
        localPart  = ct

        if (localPart .gt. 0) then
           bodyInfo % particles(FUL_PART_PROP,1:totalPart) = 0.
           bodyInfo % particles(FVL_PART_PROP,1:totalPart) = 0.
#if NDIM == MDIM
           bodyInfo % particles(FWL_PART_PROP,1:totalPart) = 0.
#endif
        endif
        localCount = localPart
        recvCount = localCount

        ct = 0
        do i = 1, totalPart ! In inverse GC fill mode particles are not sorted
                            ! in the master, we loop through all of them.   

           if (int(bodyInfo % particles(PROC_PART_PROP,i)) == bodyInfo%bodyMaster) then
#ifdef FLASH_GRID_PARAMESH
              blockID = int(bodyInfo % particles(BLK_PART_PROP,i))
#else
              blockID=CONSTANT_ONE
#endif
              particleData = bodyInfo % particles(1:NPART_PROPS,i)
#endif /* #ifndef INVERSE_GCELL_FILL */

              !call Timers_start("Grid_updateSolidBodyForces")      
              call Grid_updateSolidBodyForces(b,int(particleData(GLOB_PART_PROP)),blockID, particleData)
              !call Timers_stop("Grid_updateSolidBodyForces")

#ifndef INVERSE_GCELL_FILL
              if (i .le. localPart) &
              bodyInfo % particles(1:NPART_PROPS,i) = particleData   
#else
              ct = ct + 1
              bodyInfo % particles(1:NPART_PROPS,ct) = particleData
#endif

           end if
        end do

#ifdef INVERSE_GCELL_FILL
        if (ct .ne. recvCount) write(*,*) 'Error, Particles in master ct=',ct,'not equal to local count=',&
                                           recvCount,', Master Proc=',bodyInfo%bodyMaster 
#endif

     else ! if (bodyInfo % myPE == bodyInfo % bodyMaster) then 
        
        gettingFrom = gr_sbParticleCount(b)

#ifndef INVERSE_GCELL_FILL 

        if (localCount > 0) then
           particles2(FUL_PART_PROP,:) = 0.
           particles2(FVL_PART_PROP,:) = 0.
#if NDIM == MDIM
           particles2(FWL_PART_PROP,:) = 0.
#endif        
           recvCount = localCount           

#else /* #ifndef INVERSE_GCELL_FILL */

        recvCount = 0
        if (gettingFrom .gt. 0) then

           bodyInfo % particles(FUL_PART_PROP,1:gettingFrom) = 0.
           bodyInfo % particles(FVL_PART_PROP,1:gettingFrom) = 0.
#if NDIM == MDIM
           bodyInfo % particles(FWL_PART_PROP,1:gettingFrom) = 0.
#endif           
           recvCount = gettingFrom
           localcount= recvCount
#endif  /* #ifndef INVERSE_GCELL_FILL */

           do p = 1, recvCount

#ifndef INVERSE_GCELL_FILL 
              blockID = int(particles2(BLK_PART_PROP,p))
              particleData =  particles2(1:NPART_PROPS,p)
#else
              blockID = int(bodyInfo % particles(BLK_PART_PROP,p))
              particleData = bodyInfo % particles(1:NPART_PROPS,p)
#endif

              !call Timers_start("Grid_updateSolidBodyForces") 
              call Grid_updateSolidBodyForces(b,p,blockID, particleData)
              !call Timers_stop("Grid_updateSolidBodyForces") 

#ifndef INVERSE_GCELL_FILL
              if (p .le. gettingFrom) &
              bodyInfo % particles(1:NPART_PROPS,p) = particleData 
#else
              bodyInfo % particles(1:NPART_PROPS,p) = particleData
#endif

           enddo

        end if

     end if ! if (bodyInfo % myPE == bodyInfo % bodyMaster) then

!#ifdef INVERSE_GCELL_FILL
!     val = recvCount
!     call MPI_Allreduce(val,recvCount, CONSTANT_ONE, FLASH_INTEGER, &
!     !                   MPI_SUM, MPI_COMM_WORLD, ierr)
!     fz = 0.
!     !if (bodyInfo%bodyMaster .eq. gr_meshMe) then
!     !do i = 1,bodyInfo%totalPart
!     !   if (bodyInfo%particles(PROC_PART_PROP,i) .eq. bodyInfo%bodyMaster) &
!     !   fz = fz + &
!     !   bodyInfo%particles(FWL_PART_PROP,i)*bodyInfo%particles(HL_PART_PROP,i)*bodyInfo%particles(AREA_PART_PROP,i)
!     !enddo
!     !else
!     do i=1,val
!       fz = fz + bodyInfo%particles(FWL_PART_PROP,i)*bodyInfo%particles(HL_PART_PROP,i)*bodyInfo%particles(AREA_PART_PROP,i) 
!     enddo
!     !endif
!     fy = fz
!     call MPI_Allreduce(fy,fz, CONSTANT_ONE, FLASH_REAL,MPI_SUM, MPI_COMM_WORLD, ierr)
!     if ((b .eq. 1081) .and. (gr_meshMe .eq. 1080) ) then
!         do i=1,val
!           write(*,*) i,int(bodyInfo%particles(GLOB_PART_PROP,i)),bodyInfo%particles(FWL_PART_PROP,i), &
!                      bodyInfo%particles(HL_PART_PROP,i)*bodyInfo%particles(AREA_PART_PROP,i)
!           !gr_meshMe,'Body=',b,'Nparticles=',val,fy,fz
!         enddo
!     endif
!     if (gr_meshMe .eq. MASTER_PE) then
!       !write(*,*) gr_meshMe,'Body=',b,'Nparticles=',recvCount,fz
!       FzTotL = FzTotL + Fz
!       !write(*,*) ' ' 
!     endif
!     call mpi_barrier(MPI_COMM_WORLD, ierr)
!#endif

  Particles_Allbods = Particles_Allbods + recvCount

  end do  ! End Loop over Bodies
 call Timers_stop("Force_particles")

 call Timers_start("Barrier_UpForces")
 lb = Particles_Allbods
 call MPI_Allreduce(lb,Particles_Allbods, CONSTANT_ONE, FLASH_INTEGER,MPI_MAX, MPI_COMM_WORLD,ierr)
 call Timers_stop("Barrier_UpForces")

 if(gr_meshMe .eq. MASTER_PE) write(*,*) 'MAX Particles per Proc=',Particles_Allbods !,FzTotL

!  ! Calculate the Total Eulerian Forces for all bodies:
!  ! loop over the leaf blocks
!  Fz=0.;
!  do lb = 1,count
!     blockID = listofblocks(lb)
!     ! Get face data (velocities):
!     call Grid_getBlkPtr(blockID,facexData,FACEX)
!     call Grid_getBlkPtr(blockID,faceyData,FACEY)
!#if NDIM == MDIM
!     call Grid_getBlkPtr(blockID,facezData,FACEZ)
!#endif 
!     ! Get dx,dy,dz
!     call Grid_getDeltas(blockID,del)
!     ! The Block center coord
!     call Grid_getBlkCenterCoords(blockId,coord)
!     ! The Block Physical size
!     call Grid_getBlkPhysicalSize(blockID,bsize)
!     ! Get Blocks internal limits indexes:
!     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!     ! Local cell size:
!     dx = del(IAXIS)
!     dy = del(JAXIS)
!     dz = 1.
!#if NDIM == MDIM
!     dz = del(KAXIS)
!#endif
!     dxdydz=dx*dy*dz
!     ! Calculating the total forces 
!#if NDIM == MDIM
!     ! Z - Direction
!     Fz = SUM (facezData( FORC_FACE_VAR,                &
!          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
!          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
!          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)    &
!          ))*dxdydz + Fz
!#endif
!     ! Release face data (Forces):
!     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
!     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
!#if NDIM == MDIM
!     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
!#endif
!  enddo

!#if NDIM == MDIM
!  ! Fz:
!  FzProc = Fz
!  call MPI_Allreduce(FzProc, FzTotL, 1, FLASH_REAL,&
!       MPI_SUM, MPI_COMM_WORLD, ierr)
!#endif
!  if (gr_meshMe .eq. MASTER_PE)  then
!     write(*,*) ' '
!#if NDIM == MDIM
!     write(*,'(A35,g18.10)') ' Total Euler Fz=',FzTotL
!#endif
!  endif


 call Timers_start("Inverse_GCfill")
  !! Here Inverse guardcell filling for Forcing field:
  !! ------------------------------------------------
#ifdef INVERSE_GCELL_FILL

  ! Set To zero Gcell Force x:
  iaxis_forcexData = 0.
  jaxis_forcexData = 0.
  kaxis_forcexData = 0.

  do lb=1,count
     blockID = listOfBlocks(lb) 
     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)     
     forcexData(:,:,:,blockID)= facexData(FORC_FACE_VAR,:,:,:)
     forceyData(:,:,:,blockID)= faceyData(FORC_FACE_VAR,:,:,:) 
#if NDIM == MDIM
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     forcezData(:,:,:,blockID)= facezData(FORC_FACE_VAR,:,:,:)
#endif
     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo


  ! This is done in 2 steps:
  ! First Cross terms:
  do myaxis = IAXIS,KAXIS
     if((NDIM .ne. MDIM) .and. (myaxis .eq. KAXIS)) cycle

     ! Indexes for data to be transfered:
     ! IAXIS
     istr_F(LOW,LOW:HIGH) = GRID_ILO;   istr_F(HIGH,LOW:HIGH)= GRID_IHI+1
     istr_A(LOW,LOW:HIGH) = GRID_ILO;   istr_A(HIGH,LOW:HIGH)= GRID_IHI+1

     ! JAXIS
     jstr_F(LOW,LOW:HIGH) = GRID_JLO;   jstr_F(HIGH,LOW:HIGH)= GRID_JHI+1
     jstr_A(LOW,LOW:HIGH) = GRID_JLO;   jstr_A(HIGH,LOW:HIGH)= GRID_JHI+1

     ! KAXIS
     kstr_F(LOW,LOW:HIGH) = GRID_KLO;   kstr_F(HIGH,LOW:HIGH)= GRID_KHI+1
     kstr_A(LOW,LOW:HIGH) = GRID_KLO;   kstr_A(HIGH,LOW:HIGH)= GRID_KHI+1

     select case (myaxis)
     case(IAXIS)
       ! Limits LOW and HIGH for for LOW,HIGH boundaries, i.e. istr_F(lim,bound) 
       ! IAXIS
       ! Low Boundary
       istr_F(LOW,LOW)  = NGUARD
       istr_F(HIGH,LOW) = NGUARD+1

       istr_A(LOW,LOW)  = istr_F(LOW,LOW)  + 2
       istr_A(HIGH,LOW) = istr_F(HIGH,LOW) + 2

       ! High Boundary
       istr_F(LOW,HIGH) = GRID_IHI + 1
       istr_F(HIGH,HIGH)= GRID_IHI + 2

       istr_A(LOW,HIGH) = istr_F(LOW,HIGH) - 2
       istr_A(HIGH,HIGH)= istr_F(HIGH,HIGH)- 2

       shf_low(IAXIS) = 0
       shf_low(JAXIS) =-1
       shf_low(KAXIS) =-1

       myflg(IAXIS)   = 1
       myflg(JAXIS)   = 0
       myflg(KAXIS)   = 0

     case(JAXIS)

       ! JAXIS
       ! Low Boundary
       jstr_F(LOW,LOW)  = NGUARD
       jstr_F(HIGH,LOW) = NGUARD+1

       jstr_A(LOW,LOW)  = jstr_F(LOW,LOW)  + 2
       jstr_A(HIGH,LOW) = jstr_F(HIGH,LOW) + 2

       ! High Boundary
       jstr_F(LOW,HIGH) = GRID_JHI + 1
       jstr_F(HIGH,HIGH)= GRID_JHI + 2

       jstr_A(LOW,HIGH) = jstr_F(LOW,HIGH) - 2
       jstr_A(HIGH,HIGH)= jstr_F(HIGH,HIGH)- 2
       
       shf_low(IAXIS) =-1
       shf_low(JAXIS) = 0
       shf_low(KAXIS) =-1

       myflg(IAXIS)   = 0
       myflg(JAXIS)   = 1
       myflg(KAXIS)   = 0

#if NDIM == MDIM
     case(KAXIS)

       ! KAXIS
       ! Low Boundary
       kstr_F(LOW,LOW)  = NGUARD
       kstr_F(HIGH,LOW) = NGUARD+1
       
       kstr_A(LOW,LOW)  = kstr_F(LOW,LOW)  + 2
       kstr_A(HIGH,LOW) = kstr_F(HIGH,LOW) + 2

       ! High Boundary
       kstr_F(LOW,HIGH) = GRID_KHI + 1
       kstr_F(HIGH,HIGH)= GRID_KHI + 2

       kstr_A(LOW,HIGH) = kstr_F(LOW,HIGH) - 2
       kstr_A(HIGH,HIGH)= kstr_F(HIGH,HIGH)- 2

       shf_low(IAXIS) =-1
       shf_low(JAXIS) =-1
       shf_low(KAXIS) = 0

       myflg(IAXIS)   = 0
       myflg(JAXIS)   = 0
       myflg(KAXIS)   = 1
#endif

     end select

     ! First set to zero the AUX_FACE variable, and load the FORC field in
     ! direction myaxis:
     do lb=1,count

        blockID = listOfBlocks(lb)

        ! Get face data (velocities):
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)

        facexData(AUX_FACE_VAR,:,:,:) = 0.
        faceyData(AUX_FACE_VAR,:,:,:) = 0.

#if NDIM == MDIM
        call Grid_getBlkPtr(blockID,facezData,FACEZ)
        facezData(AUX_FACE_VAR,:,:,:) = 0.
#endif

        facexData(FORC_FACE_VAR,:,:,:) = forcexData(:,:,:,blockID)
        faceyData(FORC_FACE_VAR,:,:,:) = forceyData(:,:,:,blockID)
#if NDIM == MDIM
        facezData(FORC_FACE_VAR,:,:,:) = forcezData(:,:,:,blockID)
#endif

        ! Low Boundary:
        facexData(AUX_FACE_VAR, istr_A(LOW,LOW)+shf_low(IAXIS)*myflg(IAXIS):istr_A(HIGH,LOW)+shf_low(IAXIS)*myflg(IAXIS),   &
                                jstr_A(LOW,LOW)+shf_low(IAXIS)*myflg(JAXIS):jstr_A(HIGH,LOW)+shf_low(IAXIS)*myflg(JAXIS),   &
                                kstr_A(LOW,LOW)+shf_low(IAXIS)*myflg(KAXIS):kstr_A(HIGH,LOW)+shf_low(IAXIS)*myflg(KAXIS)) = &
        facexData(FORC_FACE_VAR,istr_F(LOW,LOW)+shf_low(IAXIS)*myflg(IAXIS):istr_F(HIGH,LOW)+shf_low(IAXIS)*myflg(IAXIS),   &
                                jstr_F(LOW,LOW)+shf_low(IAXIS)*myflg(JAXIS):jstr_F(HIGH,LOW)+shf_low(IAXIS)*myflg(JAXIS),   &
                                kstr_F(LOW,LOW)+shf_low(IAXIS)*myflg(KAXIS):kstr_F(HIGH,LOW)+shf_low(IAXIS)*myflg(KAXIS)) 


        faceyData(AUX_FACE_VAR, istr_A(LOW,LOW)+shf_low(JAXIS)*myflg(IAXIS):istr_A(HIGH,LOW)+shf_low(JAXIS)*myflg(IAXIS),   &
                                jstr_A(LOW,LOW)+shf_low(JAXIS)*myflg(JAXIS):jstr_A(HIGH,LOW)+shf_low(JAXIS)*myflg(JAXIS),   &
                                kstr_A(LOW,LOW)+shf_low(JAXIS)*myflg(KAXIS):kstr_A(HIGH,LOW)+shf_low(JAXIS)*myflg(KAXIS)) = &
        faceyData(FORC_FACE_VAR,istr_F(LOW,LOW)+shf_low(JAXIS)*myflg(IAXIS):istr_F(HIGH,LOW)+shf_low(JAXIS)*myflg(IAXIS),   &
                                jstr_F(LOW,LOW)+shf_low(JAXIS)*myflg(JAXIS):jstr_F(HIGH,LOW)+shf_low(JAXIS)*myflg(JAXIS),   &
                                kstr_F(LOW,LOW)+shf_low(JAXIS)*myflg(KAXIS):kstr_F(HIGH,LOW)+shf_low(JAXIS)*myflg(KAXIS))


#if NDIM == MDIM         
        facezData(AUX_FACE_VAR, istr_A(LOW,LOW)+shf_low(KAXIS)*myflg(IAXIS):istr_A(HIGH,LOW)+shf_low(KAXIS)*myflg(IAXIS),   &
                                jstr_A(LOW,LOW)+shf_low(KAXIS)*myflg(JAXIS):jstr_A(HIGH,LOW)+shf_low(KAXIS)*myflg(JAXIS),   &
                                kstr_A(LOW,LOW)+shf_low(KAXIS)*myflg(KAXIS):kstr_A(HIGH,LOW)+shf_low(KAXIS)*myflg(KAXIS)) = &
        facezData(FORC_FACE_VAR,istr_F(LOW,LOW)+shf_low(KAXIS)*myflg(IAXIS):istr_F(HIGH,LOW)+shf_low(KAXIS)*myflg(IAXIS),   &
                                jstr_F(LOW,LOW)+shf_low(KAXIS)*myflg(JAXIS):jstr_F(HIGH,LOW)+shf_low(KAXIS)*myflg(JAXIS),   &
                                kstr_F(LOW,LOW)+shf_low(KAXIS)*myflg(KAXIS):kstr_F(HIGH,LOW)+shf_low(KAXIS)*myflg(KAXIS))
#endif

        ! High Boundary:
        facexData(AUX_FACE_VAR, istr_A(LOW,HIGH):istr_A(HIGH,HIGH),   &
                                jstr_A(LOW,HIGH):jstr_A(HIGH,HIGH),   & 
                                kstr_A(LOW,HIGH):kstr_A(HIGH,HIGH)) = &
        facexData(FORC_FACE_VAR,istr_F(LOW,HIGH):istr_F(HIGH,HIGH),   &
                                jstr_F(LOW,HIGH):jstr_F(HIGH,HIGH),   &
                                kstr_F(LOW,HIGH):kstr_F(HIGH,HIGH))
       

        faceyData(AUX_FACE_VAR, istr_A(LOW,HIGH):istr_A(HIGH,HIGH),   &
                                jstr_A(LOW,HIGH):jstr_A(HIGH,HIGH),   &
                                kstr_A(LOW,HIGH):kstr_A(HIGH,HIGH)) = &
        faceyData(FORC_FACE_VAR,istr_F(LOW,HIGH):istr_F(HIGH,HIGH),   &
                                jstr_F(LOW,HIGH):jstr_F(HIGH,HIGH),   &
                                kstr_F(LOW,HIGH):kstr_F(HIGH,HIGH))

#if NDIM == MDIM         
        facezData(AUX_FACE_VAR, istr_A(LOW,HIGH):istr_A(HIGH,HIGH),     &
                                jstr_A(LOW,HIGH):jstr_A(HIGH,HIGH),     &
                                kstr_A(LOW,HIGH):kstr_A(HIGH,HIGH)) =   &
        facezData(FORC_FACE_VAR,istr_F(LOW,HIGH):istr_F(HIGH,HIGH),     &
                                jstr_F(LOW,HIGH):jstr_F(HIGH,HIGH),     &
                                kstr_F(LOW,HIGH):kstr_F(HIGH,HIGH))
#endif
 
        ! Release face data (velocities):
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

#if NDIM == MDIM
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     enddo

#ifdef FLASH_GRID_PARAMESH
     ! Now setup variables for guardcell fill of AUX_FACE_VAR only on direction
     ! myaxis.
     gcMask = .FALSE.
     gcMask(NUNK_VARS+AUX_FACE_VAR) = .TRUE.                 ! force x
     gcMask(NUNK_VARS+1*NFACE_VARS+AUX_FACE_VAR) = .TRUE.    ! force y
#if NDIM == 3
     gcMask(NUNK_VARS+2*NFACE_VARS+AUX_FACE_VAR) = .TRUE.    ! force z
#endif
     call Grid_fillGuardCells(CENTER_FACES,myaxis,minLayers=0,&
                              maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
#else /* #ifdef FLASH_GRID_PARAMESH */
    beginDataType=FACEX_DATATYPE
    if(NDIM==1)endDataType=FACEX_DATATYPE
    if(NDIM==2)endDataType=FACEY_DATATYPE
    if(NDIM==3)endDataType=FACEZ_DATATYPE
    do i = beginDataType,endDataType
       call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,gr_gridDataStruct(i))
       recvLeft(:,:)=1
       sendLeft(:,:)=1
       recvRight(:,:)=1
       sendRight(:,:)=1 !! do a default initialization of all starting points
       !! and then adjust individual ones as needed
       !! index of interior cell which will be first GC in block to right
       !! since the first index in the data strucutures in the variables,
       !! The "x" entry in the data structure corresponds to IAXIS+1
       j=myaxis

       sendRight(j,j+1) = blkLimits(HIGH,j)-blkLimits(LOW,j)+2-gr_offset(i,j)
        
       !index of first interior cell to be sent to be GC on block to left
       sendLeft(j, j+1) = blkLimits(LOW,j)+gr_offset(i,j)
        
       recvLeft(j, j+1) = blkLimits(HIGH,j)+1 !recv index of GC 
       call gr_shiftData(gr_axisComm(j), gr_exch(i,j), &
               sendRight(j,:), sendLeft(j,:), &
               recvRight(j,:),recvLeft(j,:),gr_gridDataStruct(i))
    end do
#endif /* #ifdef FLASH_GRID_PARAMESH */


  ! Load Resulting forces on Guard Cell arrays:
     do lb=1,count

        blockID = listOfBlocks(lb)

        ! Get face data (velocities):
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
        call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif  

        select case(myaxis)
        case(IAXIS)

        ! Force y:
        iaxis_forceyData(                      1, GRID_JLO:GRID_JHI+1, GRID_KLO:GRID_KHI,LOW,blockID) = &
               faceyData(AUX_FACE_VAR,  NGUARD-1, GRID_JLO:GRID_JHI+1, GRID_KLO:GRID_KHI)

        iaxis_forceyData(                      1, GRID_JLO:GRID_JHI+1, GRID_KLO:GRID_KHI,HIGH,blockID) = &
               faceyData(AUX_FACE_VAR,GRID_IHI+2, GRID_JLO:GRID_JHI+1, GRID_KLO:GRID_KHI)    
 
#if NDIM == MDIM
        ! Force z:
        iaxis_forcezData(                      1,   GRID_JLO:GRID_JHI, GRID_KLO:GRID_KHI+1,LOW,blockID) = &
               facezData(AUX_FACE_VAR,  NGUARD-1,   GRID_JLO:GRID_JHI, GRID_KLO:GRID_KHI+1)

        iaxis_forcezData(                      1,   GRID_JLO:GRID_JHI, GRID_KLO:GRID_KHI+1,HIGH,blockID) = &
               facezData(AUX_FACE_VAR,GRID_IHI+2,   GRID_JLO:GRID_JHI, GRID_KLO:GRID_KHI+1)        
#endif


        case(JAXIS)

        ! Force x:
        jaxis_forcexData(              GRID_ILO:GRID_IHI+1,        1, GRID_KLO:GRID_KHI,LOW,blockID) = &
               facexData(AUX_FACE_VAR, GRID_ILO:GRID_IHI+1, NGUARD-1, GRID_KLO:GRID_KHI)

        jaxis_forcexData(              GRID_ILO:GRID_IHI+1,          1, GRID_KLO:GRID_KHI,HIGH,blockID) = &
               facexData(AUX_FACE_VAR, GRID_ILO:GRID_IHI+1, GRID_JHI+2, GRID_KLO:GRID_KHI)

#if NDIM == MDIM
        ! Force z:
        jaxis_forcezData(              GRID_ILO:GRID_IHI,          1, GRID_KLO:GRID_KHI+1,LOW,blockID) = &
               facezData(AUX_FACE_VAR, GRID_ILO:GRID_IHI,   NGUARD-1, GRID_KLO:GRID_KHI+1)

        jaxis_forcezData(              GRID_ILO:GRID_IHI,          1, GRID_KLO:GRID_KHI+1,HIGH,blockID) = &        
               facezData(AUX_FACE_VAR, GRID_ILO:GRID_IHI, GRID_JHI+2, GRID_KLO:GRID_KHI+1)

        case(KAXIS)

        ! Force x:
        kaxis_forcexData(              GRID_ILO:GRID_IHI+1, GRID_JLO:GRID_JHI,        1,LOW,blockID) = &
               facexData(AUX_FACE_VAR, GRID_ILO:GRID_IHI+1, GRID_JLO:GRID_JHI, NGUARD-1)

        kaxis_forcexData(              GRID_ILO:GRID_IHI+1, GRID_JLO:GRID_JHI,          1,HIGH,blockID) = &
               facexData(AUX_FACE_VAR, GRID_ILO:GRID_IHI+1, GRID_JLO:GRID_JHI, GRID_KHI+2)

        ! Force y:
        kaxis_forceyData(              GRID_ILO:GRID_IHI, GRID_JLO:GRID_JHI+1,        1,LOW,blockID) = &
               faceyData(AUX_FACE_VAR, GRID_ILO:GRID_IHI, GRID_JLO:GRID_JHI+1, NGUARD-1)

        kaxis_forceyData(              GRID_ILO:GRID_IHI, GRID_JLO:GRID_JHI+1,          1,HIGH,blockID) = &
               faceyData(AUX_FACE_VAR, GRID_ILO:GRID_IHI, GRID_JLO:GRID_JHI+1, GRID_KHI+2)
#endif

        end select

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     enddo
  enddo

  do myaxis = IAXIS,KAXIS
     if((NDIM .ne. MDIM) .and. (myaxis .eq. KAXIS)) cycle

     ! Finally sum resulting guardcell data in AUX_FACE_VAR back to
     ! FORC_FACE_VAR.
     do lb=1,count
        blockID = listOfBlocks(lb)
        ! Get face data (velocities):
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
        call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

        select case (myaxis)
        case(IAXIS)

        ! LOW
        faceyData(FORC_FACE_VAR,GRID_ILO,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI) = &
        faceyData(FORC_FACE_VAR,GRID_ILO,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI) + &
                  iaxis_forceyData(1,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI,LOW,blockID)

#if NDIM == MDIM
        facezData(FORC_FACE_VAR,GRID_ILO,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1) = &
        facezData(FORC_FACE_VAR,GRID_ILO,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1) + &
                  iaxis_forcezData(1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1,LOW,blockID)
#endif

        ! HIGH
        faceyData(FORC_FACE_VAR,GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI) = &
        faceyData(FORC_FACE_VAR,GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI) + &
                  iaxis_forceyData(1,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI,HIGH,blockID)

#if NDIM == MDIM
        facezData(FORC_FACE_VAR,GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1) = &
        facezData(FORC_FACE_VAR,GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1) + &
                  iaxis_forcezData(1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1,HIGH,blockID)        
#endif

        case(JAXIS)

        ! LOW
        facexData(FORC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO,GRID_KLO:GRID_KHI) = &
        facexData(FORC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO,GRID_KLO:GRID_KHI) + &
                  jaxis_forcexData(GRID_ILO:GRID_IHI+1,1,GRID_KLO:GRID_KHI,LOW,blockID)

#if NDIM == MDIM
        facezData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO,GRID_KLO:GRID_KHI+1) = &
        facezData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO,GRID_KLO:GRID_KHI+1) + &
                  jaxis_forcezData(GRID_ILO:GRID_IHI,1,GRID_KLO:GRID_KHI+1,LOW,blockID)
#endif

        ! HIGH
        facexData(FORC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JHI,GRID_KLO:GRID_KHI) = &
        facexData(FORC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JHI,GRID_KLO:GRID_KHI) + &
                  jaxis_forcexData(GRID_ILO:GRID_IHI+1,1,GRID_KLO:GRID_KHI,HIGH,blockID)

#if NDIM == MDIM 
        facezData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JHI,GRID_KLO:GRID_KHI+1) = &
        facezData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JHI,GRID_KLO:GRID_KHI+1) + &
                  jaxis_forcezData(GRID_ILO:GRID_IHI,1,GRID_KLO:GRID_KHI+1,HIGH,blockID)        
#endif

        case(KAXIS)

        ! LOW
        facexData(FORC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO) = &
        facexData(FORC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO) + &
                  kaxis_forcexData(GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,1,LOW,blockID)

        faceyData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO) = &
        faceyData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO) + &
                  kaxis_forceyData(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,1,LOW,blockID)        

        ! HIGH
        facexData(FORC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KHI) = &
        facexData(FORC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KHI) + &
                  kaxis_forcexData(GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,1,HIGH,blockID)

        faceyData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KHI) = &
        faceyData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KHI) + &
                  kaxis_forceyData(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,1,HIGH,blockID)      

        end select

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

#if NDIM == MDIM
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     enddo

  enddo   


  ! Now do the collocated gc fills:
  ! First load themodified Forcing fields to the work arrays:
  do lb=1,count
     blockID = listOfBlocks(lb) 
     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)     
     forcexData(GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI,blockID) = &
     facexData(FORC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)
     forceyData(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI,blockID) = & 
     faceyData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI) 
#if NDIM == MDIM
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     forcezData(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1,blockID) = &
     facezData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1)
#endif
     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo

  !! Three directional guardcell fills:
  do myaxis = IAXIS,KAXIS
     if((NDIM .ne. MDIM) .and. (myaxis .eq. KAXIS)) cycle

     ! Indexes for data to be transfered:
     ! IAXIS
     istr_F(LOW,LOW:HIGH) = GRID_ILO;    istr_F(HIGH,LOW:HIGH)= GRID_IHI+1
     istr_A(LOW,LOW:HIGH) = GRID_ILO;    istr_A(HIGH,LOW:HIGH)= GRID_IHI+1

     ! JAXIS
     jstr_F(LOW,LOW:HIGH) = GRID_JLO;    jstr_F(HIGH,LOW:HIGH)= GRID_JHI+1
     jstr_A(LOW,LOW:HIGH) = GRID_JLO;    jstr_A(HIGH,LOW:HIGH)= GRID_JHI+1

     ! KAXIS
     kstr_F(LOW,LOW:HIGH) = GRID_KLO;    kstr_F(HIGH,LOW:HIGH)= GRID_KHI+1
     kstr_A(LOW,LOW:HIGH) = GRID_KLO;    kstr_A(HIGH,LOW:HIGH)= GRID_KHI+1

     select case (myaxis)
     case(IAXIS)
       ! Limits LOW and HIGH for for LOW,HIGH boundaries, i.e. istr_F(lim,bound) 
       ! IAXIS
       ! Low Boundary
       istr_F(LOW,LOW)  = NGUARD
       istr_F(HIGH,LOW) = NGUARD+1

       istr_A(LOW,LOW)  = istr_F(LOW,LOW)  + 2
       istr_A(HIGH,LOW) = istr_F(HIGH,LOW) + 2

       ! High Boundary
       istr_F(LOW,HIGH) = GRID_IHI + 1
       istr_F(HIGH,HIGH)= GRID_IHI + 2

       istr_A(LOW,HIGH) = istr_F(LOW,HIGH) - 2
       istr_A(HIGH,HIGH)= istr_F(HIGH,HIGH)- 2

       shf_low(IAXIS) = 0
       shf_low(JAXIS) =-1
       shf_low(KAXIS) =-1

       myflg(IAXIS)   = 1
       myflg(JAXIS)   = 0
       myflg(KAXIS)   = 0

     case(JAXIS)

       ! JAXIS
       ! Low Boundary
       jstr_F(LOW,LOW)  = NGUARD
       jstr_F(HIGH,LOW) = NGUARD+1

       jstr_A(LOW,LOW)  = jstr_F(LOW,LOW)  + 2
       jstr_A(HIGH,LOW) = jstr_F(HIGH,LOW) + 2

       ! High Boundary
       jstr_F(LOW,HIGH) = GRID_JHI + 1
       jstr_F(HIGH,HIGH)= GRID_JHI + 2

       jstr_A(LOW,HIGH) = jstr_F(LOW,HIGH) - 2
       jstr_A(HIGH,HIGH)= jstr_F(HIGH,HIGH)- 2
       
       shf_low(IAXIS) =-1
       shf_low(JAXIS) = 0
       shf_low(KAXIS) =-1

       myflg(IAXIS)   = 0
       myflg(JAXIS)   = 1
       myflg(KAXIS)   = 0

#if NDIM == MDIM
     case(KAXIS)

       ! KAXIS
       ! Low Boundary
       kstr_F(LOW,LOW)  = NGUARD
       kstr_F(HIGH,LOW) = NGUARD+1
       
       kstr_A(LOW,LOW)  = kstr_F(LOW,LOW)  + 2
       kstr_A(HIGH,LOW) = kstr_F(HIGH,LOW) + 2

       ! High Boundary
       kstr_F(LOW,HIGH) = GRID_KHI + 1
       kstr_F(HIGH,HIGH)= GRID_KHI + 2

       kstr_A(LOW,HIGH) = kstr_F(LOW,HIGH) - 2
       kstr_A(HIGH,HIGH)= kstr_F(HIGH,HIGH)- 2

       shf_low(IAXIS) =-1
       shf_low(JAXIS) =-1
       shf_low(KAXIS) = 0

       myflg(IAXIS)   = 0
       myflg(JAXIS)   = 0
       myflg(KAXIS)   = 1
#endif

     end select

     ! First set to zero the AUX_FACE variable, and load the FORC field in
     ! direction myaxis:
     do lb=1,count

        blockID = listOfBlocks(lb)

        ! Get face data (velocities):
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)

        facexData(AUX_FACE_VAR,:,:,:) = 0.
        faceyData(AUX_FACE_VAR,:,:,:) = 0.

#if NDIM == MDIM
        call Grid_getBlkPtr(blockID,facezData,FACEZ)
        facezData(AUX_FACE_VAR,:,:,:) = 0.
#endif

        facexData(FORC_FACE_VAR,:,:,:) = forcexData(:,:,:,blockID)
        faceyData(FORC_FACE_VAR,:,:,:) = forceyData(:,:,:,blockID)
        facezData(FORC_FACE_VAR,:,:,:) = forcezData(:,:,:,blockID)

        ! Low Boundary:
        facexData(AUX_FACE_VAR, istr_A(LOW,LOW)+shf_low(IAXIS)*myflg(IAXIS):istr_A(HIGH,LOW)+shf_low(IAXIS)*myflg(IAXIS),   &
                                jstr_A(LOW,LOW)+shf_low(IAXIS)*myflg(JAXIS):jstr_A(HIGH,LOW)+shf_low(IAXIS)*myflg(JAXIS),   &
                                kstr_A(LOW,LOW)+shf_low(IAXIS)*myflg(KAXIS):kstr_A(HIGH,LOW)+shf_low(IAXIS)*myflg(KAXIS)) = &
        facexData(FORC_FACE_VAR,istr_F(LOW,LOW)+shf_low(IAXIS)*myflg(IAXIS):istr_F(HIGH,LOW)+shf_low(IAXIS)*myflg(IAXIS),   &
                                jstr_F(LOW,LOW)+shf_low(IAXIS)*myflg(JAXIS):jstr_F(HIGH,LOW)+shf_low(IAXIS)*myflg(JAXIS),   &
                                kstr_F(LOW,LOW)+shf_low(IAXIS)*myflg(KAXIS):kstr_F(HIGH,LOW)+shf_low(IAXIS)*myflg(KAXIS)) 


        faceyData(AUX_FACE_VAR, istr_A(LOW,LOW)+shf_low(JAXIS)*myflg(IAXIS):istr_A(HIGH,LOW)+shf_low(JAXIS)*myflg(IAXIS),   &
                                jstr_A(LOW,LOW)+shf_low(JAXIS)*myflg(JAXIS):jstr_A(HIGH,LOW)+shf_low(JAXIS)*myflg(JAXIS),   &
                                kstr_A(LOW,LOW)+shf_low(JAXIS)*myflg(KAXIS):kstr_A(HIGH,LOW)+shf_low(JAXIS)*myflg(KAXIS)) = &
        faceyData(FORC_FACE_VAR,istr_F(LOW,LOW)+shf_low(JAXIS)*myflg(IAXIS):istr_F(HIGH,LOW)+shf_low(JAXIS)*myflg(IAXIS),   &
                                jstr_F(LOW,LOW)+shf_low(JAXIS)*myflg(JAXIS):jstr_F(HIGH,LOW)+shf_low(JAXIS)*myflg(JAXIS),   &
                                kstr_F(LOW,LOW)+shf_low(JAXIS)*myflg(KAXIS):kstr_F(HIGH,LOW)+shf_low(JAXIS)*myflg(KAXIS))


#if NDIM == MDIM         
        facezData(AUX_FACE_VAR, istr_A(LOW,LOW)+shf_low(KAXIS)*myflg(IAXIS):istr_A(HIGH,LOW)+shf_low(KAXIS)*myflg(IAXIS),   &
                                jstr_A(LOW,LOW)+shf_low(KAXIS)*myflg(JAXIS):jstr_A(HIGH,LOW)+shf_low(KAXIS)*myflg(JAXIS),   &
                                kstr_A(LOW,LOW)+shf_low(KAXIS)*myflg(KAXIS):kstr_A(HIGH,LOW)+shf_low(KAXIS)*myflg(KAXIS)) = &
        facezData(FORC_FACE_VAR,istr_F(LOW,LOW)+shf_low(KAXIS)*myflg(IAXIS):istr_F(HIGH,LOW)+shf_low(KAXIS)*myflg(IAXIS),   &
                                jstr_F(LOW,LOW)+shf_low(KAXIS)*myflg(JAXIS):jstr_F(HIGH,LOW)+shf_low(KAXIS)*myflg(JAXIS),   &
                                kstr_F(LOW,LOW)+shf_low(KAXIS)*myflg(KAXIS):kstr_F(HIGH,LOW)+shf_low(KAXIS)*myflg(KAXIS))
#endif

        ! High Boundary:
        facexData(AUX_FACE_VAR, istr_A(LOW,HIGH):istr_A(HIGH,HIGH),   &
                                jstr_A(LOW,HIGH):jstr_A(HIGH,HIGH),   & 
                                kstr_A(LOW,HIGH):kstr_A(HIGH,HIGH)) = &
        facexData(FORC_FACE_VAR,istr_F(LOW,HIGH):istr_F(HIGH,HIGH),   &
                                jstr_F(LOW,HIGH):jstr_F(HIGH,HIGH),   &
                                kstr_F(LOW,HIGH):kstr_F(HIGH,HIGH))
       

        faceyData(AUX_FACE_VAR, istr_A(LOW,HIGH):istr_A(HIGH,HIGH),   &
                                jstr_A(LOW,HIGH):jstr_A(HIGH,HIGH),   &
                                kstr_A(LOW,HIGH):kstr_A(HIGH,HIGH)) = &
        faceyData(FORC_FACE_VAR,istr_F(LOW,HIGH):istr_F(HIGH,HIGH),   &
                                jstr_F(LOW,HIGH):jstr_F(HIGH,HIGH),   &
                                kstr_F(LOW,HIGH):kstr_F(HIGH,HIGH))

#if NDIM == MDIM         
        facezData(AUX_FACE_VAR, istr_A(LOW,HIGH):istr_A(HIGH,HIGH),     &
                                jstr_A(LOW,HIGH):jstr_A(HIGH,HIGH),     &
                                kstr_A(LOW,HIGH):kstr_A(HIGH,HIGH)) =   &
        facezData(FORC_FACE_VAR,istr_F(LOW,HIGH):istr_F(HIGH,HIGH),     &
                                jstr_F(LOW,HIGH):jstr_F(HIGH,HIGH),     &
                                kstr_F(LOW,HIGH):kstr_F(HIGH,HIGH))
#endif
 
        ! Release face data (velocities):
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

#if NDIM == MDIM
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     enddo

#ifdef FLASH_GRID_PARAMESH
     ! Now setup variables for guardcell fill of AUX_FACE_VAR only on direction
     ! myaxis.
     gcMask = .FALSE.
     gcMask(NUNK_VARS+AUX_FACE_VAR) = .TRUE.                 ! force x
     gcMask(NUNK_VARS+1*NFACE_VARS+AUX_FACE_VAR) = .TRUE.    ! force y
#if NDIM == 3
     gcMask(NUNK_VARS+2*NFACE_VARS+AUX_FACE_VAR) = .TRUE.    ! force z
#endif
     call Grid_fillGuardCells(CENTER_FACES,myaxis,minLayers=0,&
                              maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
#else /*  #ifdef FLASH_GRID_PARAMESH */
    beginDataType=FACEX_DATATYPE
    if(NDIM==1)endDataType=FACEX_DATATYPE
    if(NDIM==2)endDataType=FACEY_DATATYPE
    if(NDIM==3)endDataType=FACEZ_DATATYPE
    do i = beginDataType,endDataType
       call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,gr_gridDataStruct(i))
       recvLeft(:,:)=1
       sendLeft(:,:)=1
       recvRight(:,:)=1
       sendRight(:,:)=1 !! do a default initialization of all starting points
       !! and then adjust individual ones as needed
       !! index of interior cell which will be first GC in block to right
       !! since the first index in the data strucutures in the variables,
       !! The "x" entry in the data structure corresponds to IAXIS+1
       j=myaxis

       sendRight(j,j+1) = blkLimits(HIGH,j)-blkLimits(LOW,j)+2-gr_offset(i,j)
        
       !index of first interior cell to be sent to be GC on block to left
       sendLeft(j, j+1) = blkLimits(LOW,j)+gr_offset(i,j)
        
       recvLeft(j, j+1) = blkLimits(HIGH,j)+1 !recv index of GC 
       call gr_shiftData(gr_axisComm(j), gr_exch(i,j), &
               sendRight(j,:), sendLeft(j,:), &
               recvRight(j,:),recvLeft(j,:),gr_gridDataStruct(i))
    end do
#endif /*  #ifdef FLASH_GRID_PARAMESH */

  ! Load Resulting forces on Guard Cell arrays:
     do lb=1,count

        blockID = listOfBlocks(lb)

        ! Get face data (velocities):
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
        call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif  

        select case(myaxis)
        case(IAXIS)

        ! Force x:        
        iaxis_forcexData(                             1:NGUARD, GRID_JLO:GRID_JHI, GRID_KLO:GRID_KHI,LOW,blockID)  = &
               facexData(AUX_FACE_VAR,                1:NGUARD, GRID_JLO:GRID_JHI, GRID_KLO:GRID_KHI)

        iaxis_forcexData(                             1:NGUARD, GRID_JLO:GRID_JHI, GRID_KLO:GRID_KHI,HIGH,blockID) = &
               facexData(AUX_FACE_VAR,GRID_IHI+2:GRID_IHI_GC+1, GRID_JLO:GRID_JHI, GRID_KLO:GRID_KHI)

        case(JAXIS)

        ! Force y:
        jaxis_forceyData(              GRID_ILO:GRID_IHI, 1:NGUARD, GRID_KLO:GRID_KHI,LOW,blockID) = &
               faceyData(AUX_FACE_VAR, GRID_ILO:GRID_IHI, 1:NGUARD, GRID_KLO:GRID_KHI)

        jaxis_forceyData(              GRID_ILO:GRID_IHI,                 1:NGUARD, GRID_KLO:GRID_KHI,HIGH,blockID) = &
               faceyData(AUX_FACE_VAR, GRID_ILO:GRID_IHI, GRID_JHI+2:GRID_JHI_GC+1, GRID_KLO:GRID_KHI)          

#if NDIM == MDIM

        case(KAXIS)

        ! Force z:
        kaxis_forcezData(              GRID_ILO:GRID_IHI, GRID_JLO:GRID_JHI, 1:NGUARD,LOW,blockID) = &
               facezData(AUX_FACE_VAR, GRID_ILO:GRID_IHI, GRID_JLO:GRID_JHI, 1:NGUARD)

        kaxis_forcezData(              GRID_ILO:GRID_IHI, GRID_JLO:GRID_JHI,                 1:NGUARD,HIGH,blockID) = &
               facezData(AUX_FACE_VAR, GRID_ILO:GRID_IHI, GRID_JLO:GRID_JHI, GRID_KHI+2:GRID_KHI_GC+1)

#endif

        end select

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     enddo
  enddo

  do myaxis = IAXIS,KAXIS
     if((NDIM .ne. MDIM) .and. (myaxis .eq. KAXIS)) cycle

     ! Finally sum resulting guardcell data in AUX_FACE_VAR back to
     ! FORC_FACE_VAR.
     do lb=1,count
        blockID = listOfBlocks(lb)
        ! Get face data (velocities):
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
        call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

        select case (myaxis)
        case(IAXIS)

        ! LOW
        facexData(FORC_FACE_VAR,GRID_ILO:GRID_ILO+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI) = &
        facexData(FORC_FACE_VAR,GRID_ILO:GRID_ILO+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI) + &
                  iaxis_forcexData(NGUARD-1:NGUARD,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI,LOW,blockID)

        ! HIGH
        facexData(FORC_FACE_VAR,GRID_IHI:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI) = &
        facexData(FORC_FACE_VAR,GRID_IHI:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI) + &
                  iaxis_forcexData(1:2,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI,HIGH,blockID)        

        case(JAXIS)

        ! LOW
        faceyData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JLO+1,GRID_KLO:GRID_KHI) = &
        faceyData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JLO+1,GRID_KLO:GRID_KHI) + &
                  jaxis_forceyData(GRID_ILO:GRID_IHI,NGUARD-1:NGUARD,GRID_KLO:GRID_KHI,LOW,blockID)

        ! HIGH
        faceyData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JHI:GRID_JHI+1,GRID_KLO:GRID_KHI) = &
        faceyData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JHI:GRID_JHI+1,GRID_KLO:GRID_KHI) + &
                  jaxis_forceyData(GRID_ILO:GRID_IHI,1:2,GRID_KLO:GRID_KHI,HIGH,blockID)

        case(KAXIS)

        ! LOW
#if NDIM == MDIM
        facezData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KLO+1) = &
        facezData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KLO+1) + &
                  kaxis_forcezData(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,NGUARD-1:NGUARD,LOW,blockID)
#endif

        ! HIGH
#if NDIM == MDIM
        facezData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KHI:GRID_KHI+1) = &
        facezData(FORC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KHI:GRID_KHI+1) + &
                  kaxis_forcezData(GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,1:2,HIGH,blockID)
#endif
        end select

        ! Release face data (velocities):
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     enddo
  enddo   
#endif /*  #ifdef INVERSE_GCELL_FILL */
  call Timers_stop("Inverse_GCfill")

#endif /* DO_MARKER_FORCING */

  ! Force inside Bodies:
  call Timers_start("ib_forceInsideBody")
  call ib_forceInsideBody()
  call Timers_stop("ib_forceInsideBody")


  !! Do the sum of forces to Velocities:
  !! ----------------------------------
  do lb=1,count

     blockID = listOfBlocks(lb)

     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 

     ! X velocities:
     facexData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1, &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &                       
     facexData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1, &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) + &  
     ib_dt*                                                                &
     facexData(FORC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1, &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))

     ! Y velocities:
     faceyData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &                       
     faceyData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) + &  
     ib_dt*                                                                &
     faceyData(FORC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))

#ifdef ITER_FORCING

     facexData(FORO_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1, &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &                       
     facexData(FORO_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1, &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) + &  
     facexData(FORC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1, &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))

     faceyData(FORO_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &                       
     faceyData(FORO_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) + &  
     faceyData(FORC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
#endif
#if NDIM == MDIM
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     ! Z velocities:
     facezData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1) = &                       
     facezData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1) + &  
     ib_dt*                                                                &
     facezData(FORC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)
#ifdef ITER_FORCING
     ! Z velocities:
     facezData(FORO_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1) = &                       
     facezData(FORO_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1) + &  
     facezData(FORC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)
#endif
#endif

     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo


! Save the Lagrangian forcing
#ifdef ITER_FORCING

  do b = 1, gr_sbNumBodies
    bodyInfo => gr_sbBodyInfo(b)
    if(gr_sbParticleCount(b) .gt. 0) then    
    bodyInfo % particles(FULO_PART_PROP,:) = & 
       bodyInfo % particles(FULO_PART_PROP,:) + &
       bodyInfo % particles(FUL_PART_PROP,:)

    bodyInfo % particles(FVLO_PART_PROP,:) = & 
       bodyInfo % particles(FVLO_PART_PROP,:) + & 
       bodyInfo % particles(FVL_PART_PROP,:) 
#if NDIM == MDIM
    bodyInfo % particles(FWLO_PART_PROP,:) = &
       bodyInfo % particles(FWLO_PART_PROP,:) + &
       bodyInfo % particles(FWL_PART_PROP,:) 
#endif
    endif
  enddo

#endif


#ifdef ITER_FORCING
! Fill guard cells
   ! APPLY BC AND FILL GUARDCELLS FOR INTERMEDIATE VELOCITIES:
   ! ----- -- --- ---- ---------- --- ------------ ----------
   gcMask = .FALSE.
   gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! ustar
   gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! vstar
#if NDIM == 3
   gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! wstar
#endif
   call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
        maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
 
   ! FIX FLUXES FOR USTAR: (Only for AMR grids)
   ! --- ------ --- -----
#ifdef FLASH_GRID_PARAMESH
   ! Fix fluxes at block boundaries
   call ins_fluxfix(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                    count,listofBlocks)
#endif

#endif

#ifdef CHK_BND_VLC
    !write(*,*) 'CHK_BND_VLC', CHK_BND_VLC
     ! Output interpolated velocity    
    call Driver_getNStep(NStep)
    if( mod(NStep,CHK_BND_VLC) == 0) then
    write(indNstep, '(I6.6)') NStep/CHK_BND_VLC
    call gr_intepVel('parts-uv-step-'//indNStep//'-iter',iForceIter)
    endif
#endif

#ifdef ITER_FORCING
  enddo ! iForceIter. Iiterative Forcing
#endif

  ! Recover variables
#ifdef ITER_FORCING

  do lb=1,count

     blockID = listOfBlocks(lb)

     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     facexData(FORC_FACE_VAR,:,:,:) = facexData(FORO_FACE_VAR,:,:,:)
     faceyData(FORC_FACE_VAR,:,:,:) = faceyData(FORO_FACE_VAR,:,:,:)

#if NDIM == MDIM
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     facezData(FORC_FACE_VAR,:,:,:) = facezData(FORO_FACE_VAR,:,:,:)
#endif

     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo

  do b = 1, gr_sbNumBodies
    bodyInfo => gr_sbBodyInfo(b)
    if(gr_sbParticleCount(b) .gt. 0) then    
    bodyInfo % particles(FUL_PART_PROP,:) = & 
       bodyInfo % particles(FULO_PART_PROP,:) 
    bodyInfo % particles(FVL_PART_PROP,:) = & 
       bodyInfo % particles(FVLO_PART_PROP,:) 
#if NDIM == MDIM
    bodyInfo % particles(FWL_PART_PROP,:) = &
       bodyInfo % particles(FWLO_PART_PROP,:) 
#endif
    endif
  enddo

#endif


#ifdef TEST_IB_FORCES

  ! ----------------------------------------------------------------------------
  ! Calculate the Lagrangian Forces on each body, only Stationary Bodies  

  ! Body Center: Use Origin of Eulerian Coordinate System
  Body_Cen(IAXIS) = 0. 
  Body_Cen(JAXIS) = 0. 
  Body_Cen(KAXIS) = 0.
  do b = 1, gr_sbNumBodies
     bodyInfo => gr_sbBodyInfo(b)

     Fxtot(b) = 0.
     Fytot(b) = 0.
     Momz(b)  = 0.

     Fxtoti(b) = 0.
     Fytoti(b) = 0.
     Momzi(b)  = 0.

     Fztot(b)  = 0. 
     Fztoti(b) = 0.
     Momx(b)   = 0.
     Momxi(b)  = 0.
     Momy(b)   = 0.
     Momyi(b)  = 0. 
  
     MomArmi_x=0.; MomArmi_y=0.; MomArmi_z=0.;  

     ! Here we could compute another Body_Cen..

     if (bodyInfo%myPE .eq. bodyInfo % bodymaster ) then
     ! Send Body center:
     call mpi_bcast(Body_Cen, NDIM, FLASH_REAL, bodyInfo%bodyMaster, gr_meshComm, ierr)
     else
     ! Receive body center:
     call mpi_bcast(Body_Cen, NDIM, FLASH_REAL, bodyInfo%bodyMaster, gr_meshComm, ierr)
     endif

     if (bodyInfo % myPE == bodyInfo % bodyMaster) then

        ! Find number of particles in the Master Processor:
        localPart = 0    
        do j=1,size(bodyInfo%particlesPerProc,DIM=2)
           if(bodyInfo%particlesPerProc(1,j) .eq. bodyInfo%bodyMaster) then
              localPart = bodyInfo%particlesPerProc(2,j)
           endif
        enddo

        do i = 1, localPart 
           if (int(bodyInfo % particles(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then

              MomArmi_x = bodyInfo%particles(POSX_PART_PROP,i)-Body_Cen(IAXIS)
              MomArmi_y = bodyInfo%particles(POSY_PART_PROP,i)-Body_Cen(JAXIS)

              Voli =  bodyInfo%particles(AREA_PART_PROP,i)*bodyInfo%particles(HL_PART_PROP,i)

              Fxtoti(b) =  Fxtoti(b) + bodyInfo%particles(FUL_PART_PROP,i)* Voli
              Fytoti(b) =  Fytoti(b) + bodyInfo%particles(FVL_PART_PROP,i)* Voli

              Momzi (b) =  Momzi(b)  - bodyInfo%particles(FUL_PART_PROP,i)*MomArmi_y*Voli +  &
                                       bodyInfo%particles(FVL_PART_PROP,i)*MomArmi_x*Voli

#if NDIM == MDIM
              MomArmi_z = bodyInfo%particles(POSZ_PART_PROP,i)-Body_Cen(KAXIS)
              Fztoti(b) = Fztoti(b) + bodyInfo%particles(FWL_PART_PROP,i)* Voli
              Momxi (b) =  Momxi(b) - bodyInfo%particles(FVL_PART_PROP,i)*MomArmi_z*Voli +  &
                                      bodyInfo%particles(FWL_PART_PROP,i)*MomArmi_y*Voli              
              Momyi (b) =  Momyi(b) - bodyInfo%particles(FWL_PART_PROP,i)*MomArmi_x*Voli +  &
                                      bodyInfo%particles(FUL_PART_PROP,i)*MomArmi_z*Voli           
#endif 
           endif
        enddo

     else !MyPe not BodyMaster

        gettingFrom = gr_sbParticleCount(b)
        if (gettingFrom > 0) then
          recvCount = gettingFrom
           do p = 1, recvCount

              MomArmi_x = bodyInfo%particles(POSX_PART_PROP,p)-Body_Cen(IAXIS)
              MomArmi_y = bodyInfo%particles(POSY_PART_PROP,p)-Body_Cen(JAXIS)

              Voli =  bodyInfo%particles(AREA_PART_PROP,p)*bodyInfo%particles(HL_PART_PROP,p)

              Fxtoti(b) =  Fxtoti(b) + bodyInfo%particles(FUL_PART_PROP,p)* Voli
              Fytoti(b) =  Fytoti(b) + bodyInfo%particles(FVL_PART_PROP,p)* Voli
              Momzi (b) =  Momzi(b)  - bodyInfo%particles(FUL_PART_PROP,p)*MomArmi_y*Voli +  &
                                       bodyInfo%particles(FVL_PART_PROP,p)*MomArmi_x*Voli

#if NDIM == MDIM
              MomArmi_z = bodyInfo%particles(POSZ_PART_PROP,p)-Body_Cen(KAXIS)
              Fztoti(b) = Fztoti(b) + bodyInfo%particles(FWL_PART_PROP,p)* Voli
              Momxi (b) =  Momxi(b) - bodyInfo%particles(FVL_PART_PROP,p)*MomArmi_z*Voli +  &
                                      bodyInfo%particles(FWL_PART_PROP,p)*MomArmi_y*Voli
              Momyi (b) =  Momyi(b) - bodyInfo%particles(FWL_PART_PROP,p)*MomArmi_x*Voli +  &
                                      bodyInfo%particles(FUL_PART_PROP,p)*MomArmi_z*Voli

#endif 
           enddo
        endif
     endif

     ! Sum force values among processors for each body, only stationary Bodies
     ! Fx:
     call MPI_Allreduce(Fxtoti(b),Fxtot(b), CONSTANT_ONE, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)     

     ! Fy:
     call MPI_Allreduce(Fytoti(b),Fytot(b), CONSTANT_ONE, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)   

     ! Mz:
     call MPI_Allreduce(Momzi(b),Momz(b), CONSTANT_ONE, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)   

#if NDIM == MDIM
     ! Fz:
     call MPI_Allreduce(Fztoti(b),Fztot(b), CONSTANT_ONE, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)   

     ! Mx:
     call MPI_Allreduce(Momxi(b),Momx(b), CONSTANT_ONE, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)   

     ! My:
     call MPI_Allreduce(Momyi(b),Momy(b), CONSTANT_ONE, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)   
#endif

     if (gr_meshMe .eq. MASTER_PE)  then
     write(*,*) ' '
     write(*,'(A35,3g18.10,A10,I8)') ' Total Lagr Fx , Fy & Fz=',FxTot(b), FyTot(b), FzTot(b),'on body ',b
     write(*,'(A35,3g18.10,A10,I8)') ' Total Lagr Mx , My & Mz=', Momx(b),  Momy(b),  Momz(b),'on body ',b
     endif

  enddo


  ! Calculate the Total Eulerian Forces for all bodies:
  ! loop over the leaf blocks
  Fx = 0.; Fy=0.; Fz=0.; MomxL=0.; MomyL=0.; MomzL=0.;
  do lb = 1,count

     blockID = listofblocks(lb)
     
     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == MDIM
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif 
     
     ! Get dx,dy,dz
     call Grid_getDeltas(blockID,del)
     
     ! The Block center coord
     call Grid_getBlkCenterCoords(blockId,coord)
     
     ! The Block Physical size
     call Grid_getBlkPhysicalSize(blockID,bsize)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 
     
     ! Local cell size:
     dx = del(IAXIS)
     dy = del(JAXIS)
     dz = 1.
#if NDIM == MDIM
     dz = del(KAXIS)
#endif
     dxdydz=dx*dy*dz
     
     ! Calculating the total forces 

     ! X - Direction
     Fx = SUM (facexData(FORC_FACE_VAR,                 &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)    &
          ))*dxdydz + Fx 
     
     ! Y - Direction
     Fy = SUM (faceyData( FORC_FACE_VAR,                &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)    &
          ))*dxdydz + Fy      

#if NDIM == MDIM
     ! Z - Direction
     Fz = SUM (facezData( FORC_FACE_VAR,                &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)    &
          ))*dxdydz + Fz    
#endif
     
     ! Calc the moment about the Z (2D) (+ve is anti clock-wise)
     do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)

        zcell = coord(KAXIS) - 0.5*bsize(KAXIS) + &
                real(k-NGUARD-1)*dz + 0.5*dz
 
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
   
           ycell = coord(JAXIS) - 0.5*bsize(JAXIS) + &
                   real(j-NGUARD-1)*dy + 0.5*dy

           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

              xcell = coord(IAXIS) - 0.5*bsize(IAXIS) + &
                      real(i-NGUARD-1)*dx + 0.5*dx

              MomzL = -dxdydz*facexData(FORC_FACE_VAR,i,j,k)*(ycell-0.) + &
                       dxdydz*faceyData(FORC_FACE_VAR,i,j,k)*(xcell-0.) + &
                       MomzL

#if NDIM == MDIM
              MomxL = -dxdydz*faceyData(FORC_FACE_VAR,i,j,k)*(zcell-0.) + &
                       dxdydz*facezData(FORC_FACE_VAR,i,j,k)*(ycell-0.) + &
                       MomxL 

              MomyL = -dxdydz*facezData(FORC_FACE_VAR,i,j,k)*(xcell-0.) + &
                       dxdydz*facexData(FORC_FACE_VAR,i,j,k)*(zcell-0.) + &
                       MomyL 
#endif

           enddo
        enddo
     enddo

     ! Release face data (Forces):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     
#if NDIM == MDIM
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     
     
  enddo
  
  ! Fx:
  FxProc = Fx
  call MPI_Allreduce(FxProc, FxTotL, 1, FLASH_REAL,&
       MPI_SUM, MPI_COMM_WORLD, ierr)
  
  ! Fy:
  FyProc = Fy  
  call MPI_Allreduce(FyProc, FyTotL, 1, FLASH_REAL,&
       MPI_SUM, MPI_COMM_WORLD, ierr)
  
  ! Mz:
  Momzproc = MomzL
  call MPI_Allreduce(Momzproc,MomzL, 1, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)
 
#if NDIM == MDIM

  ! Fz:
  FzProc = Fz
  call MPI_Allreduce(FzProc, FzTotL, 1, FLASH_REAL,&
       MPI_SUM, MPI_COMM_WORLD, ierr)
  
  ! Mx:
  MomxProc = MomxL  
  call MPI_Allreduce(MomxProc, MomxL, 1, FLASH_REAL,&
       MPI_SUM, MPI_COMM_WORLD, ierr)
  
  ! My:
  MomyProc = MomyL
  call MPI_Allreduce(MomyProc,MomyL, 1, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)

#endif

  if (gr_meshMe .eq. MASTER_PE)  then
     write(*,*) ' '
#if NDIM == MDIM
     write(*,'(A35,3g18.10)') ' Total Euler Fx , Fy & Fz=',FxTotL, FyTotL, FzTotL
     write(*,'(A35,3g18.10)') ' Total Euler Mx , My & Mz=',MomxL,  MomyL, MomzL
#else
     write(*,'(A35,3g18.10)') ' Total Euler Fx , Fy =',FxTotL, FyTotL
     write(*,'(A35,3g18.10)') ' Total Euler Mz=',MomzL
#endif

  endif

#endif /* TEST_IB_FORCES */

!  call Timers_stop("update_forces")
end subroutine gr_sbUpdateForces

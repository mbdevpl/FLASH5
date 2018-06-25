
!!             integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount)
!!             real(IN)    :: dt)
! This routine computes the drag, lift and moment about Z for every immersed body a 
! given time step. 


SUBROUTINE ib_CalcForce(blockCount, blockList, dt)

  use Driver_data, only : dr_simTime

  use Grid_data, ONLY: gr_meshMe

  use ImBound_data, ONLY: ib_nbd,ib_dsxu2,ib_dsyv2,ib_FuL2,ib_indFuL,ib_FvL2,ib_indFvL, &
                          ib_xb,ib_yb,ib_sb,ib_alphax,ib_alphay,ib_ABODY,ib_AELEM,ib_Fd, &
                          ib_xbus,ib_ybus,freq_t,Ro,ao        

  use Grid_interface, ONLY : Grid_getDeltas,           &
                             Grid_getBlkCenterCoords,  &
                             Grid_getBlkPtr,           &
                             Grid_releaseBlkPtr,       &
                             Grid_getBlkIndexLimits,   &
                             Grid_getBlkPhysicalSize

  use Driver_interface, ONLY : Driver_abortFlash

  use ib_interface, only : ib_forces2D

  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
  real, INTENT(IN) :: dt
  !! -----------------------------------------------------


  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData

  ! Local variables
  integer :: lb, blockID    ! block counter
  integer :: ierr  ! error flag for the MPI_allreduce call
  integer :: ibd   ! Counter over the number of immersed bodies
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real    :: dx , dy, dxdy ! local grid spacing in x and y direction
  real    :: del(MDIM), coord(MDIM), bsize(MDIM)
  real    :: Fx, Fy, FxProc, FyProc, Momz, Momzproc 
  real    :: FxTot, FyTot, Clift, Cdrag, Cmoment
  real    :: ib_BodyCen(MDIM,ib_nbd)
  real    :: xcell,ycell
  integer :: i,j,k

  real :: dsbm,dsbp,dsb,hl,volL
  real :: FxL,FyL,MomzL

  integer :: ind1,ind2

  logical, save :: firstflg = .true.


  integer :: ind1e, ind2e
  real :: tita
  real Qx,Qy,Qtita
  real Lseg,t(MDIM), kfact, kfact2
  real :: tau1tau2_1,tau1tau2_2,pr1pr2_1,pr1pr2_2
  real :: t1mt2_1,t2pt1_1,t1mt2_2,t2pt1_2 
  real :: dQtitasdrp1x,dQtitasdrp1y,dQtitasdrp2x,dQtitasdrp2y
  real :: Qxs,Qys,Qtitas


  ! Center of the immersed boundary i
  ibd = 1
  ib_BodyCen(IAXIS,ibd)= 0.
  ib_BodyCen(JAXIS,ibd)= 0.


  ! Eulerian grid computation of forces:

  ! loop over the leaf blocks
  Fx = 0.; Fy=0.; Momz=0.;
  do lb = 1,blockCount

     blockID = blockList(lb)
     
     ! Get face data (velocities):
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif 
     
     ! Get dx,dy
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
     dxdy=dx*dy
     
     ! Calculating the total forces 

     ! X - Direction
     Fx = SUM (facexData(FORC_FACE_VAR,                 &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)    &
          ))*dxdy + Fx 
     
     ! Y - Direction
     Fy = SUM (faceyData( FORC_FACE_VAR,                &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)    &
          ))*dxdy + Fy      
     
     ! Calc the moment about the Z (2D) (+ve is anti clock-wise)
     do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)

!!$           yface = coord(JAXIS) - 0.5*bsize(JAXIS) + &
!!$                   real(j-NGUARD-1)*dy
   
           ycell = coord(JAXIS) - 0.5*bsize(JAXIS) + &
                   real(j-NGUARD-1)*dy + 0.5*dy

           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

!!$              xface = coord(IAXIS) - 0.5*bsize(IAXIS) + &
!!$                      real(i-NGUARD-1)*dx

              xcell = coord(IAXIS) - 0.5*bsize(IAXIS) + &
                      real(i-NGUARD-1)*dx + 0.5*dx

              Momz = -dxdy*facexData(FORC_FACE_VAR,i,j,k)*(ycell-ib_BodyCen(JAXIS,ibd)) + &
                      dxdy*faceyData(FORC_FACE_VAR,i,j,k)*(xcell-ib_BodyCen(IAXIS,ibd)) + &
                      Momz

           enddo
        enddo
     enddo

!!$     MomZ= -Fx*( coord(JAXIS)-ib_BodyCen(JAXIS,ibd))  +  & 
!!$            Fy*( coord(IAXIS)-ib_BodyCen(IAXIS,ibd))) +  &
!!$            MomZ


     ! Release face data (Forces):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     
     
  enddo
  
  FxProc = Fx
!  FxTot  = FxProc 
  call MPI_Allreduce(FxProc, FxTot, 1, FLASH_REAL,&
       MPI_SUM, MPI_COMM_WORLD, ierr)
  
  FyProc = Fy
!  FyTot  = FyProc  
  call MPI_Allreduce(FyProc, FyTot, 1, FLASH_REAL,&
       MPI_SUM, MPI_COMM_WORLD, ierr)

  Momzproc = Momz
  call MPI_Allreduce(Momzproc,Momz, 1, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)
 
  ! Calculating the drag and lift, and Moment coefficients
  Cdrag  = 2.* FxTot  ! Reference length, velocity and density == 1
  Clift  = 2.* FyTot
  Cmoment= 2.* MomZ

  if (gr_meshMe .eq. MASTER_PE)  then
     write(*,*) ' '
     write(*,'(A35,3g18.12)') ' Total Euler Fx , Fy & Mz=',FxTot, FyTot, Momz
     write(*,'(A35,3g18.12)') ' Euler Drag , Lift & Moment ceoff=',Cdrag, Clift,Cmoment
  endif
  
  if (firstflg) then

     if (gr_meshMe .eq. MASTER_PE)  then
     open(unit=113,file='./IOData/Cyl_forcesEuler.res',form='formatted',status='replace')
     write(113,'(4g16.8)')dr_simTime,FxTot,FyTot,Momz
     close(113)
     endif

     !firstflg = .false.

  else

     if (gr_meshMe .eq. MASTER_PE)  then
     open(unit=113,file='./IOData/Cyl_forcesEuler.res',form='formatted',status='old',position='append')
     write(113,'(4g16.8)')dr_simTime,FxTot,FyTot,Momz
     close(113)
     endif

  endif  
  

  ! Lagrangian Grid Force computation:
  FxL=0.; FyL=0.; MomzL=0.;
  do ibd=1,ib_nbd

     ! Find bounding box of body
     ind1 = ib_ABODY(ibd) % lb + 1
     ind2 = ib_ABODY(ibd) % lb + ib_ABODY(ibd) % mb

     do i = ind1,ind2

        
        if(i .eq. ind2) then
        dsbp = 0.
        else
        dsbp = ib_sb(i+1)-ib_sb(i)
        endif

        if(i .eq. ind1) then
        dsbm = 0.
        else
        dsbm = ib_sb(i)-ib_sb(i-1) 
        endif

        dsb = .5*(dsbm+dsbp)

        !write(*,*) 'dsb=',dsb


        if (ib_indFuL(i) .eq. 0.) call Driver_abortFlash(" Lagrangian CalcForce: Marker U with no stencil !!")
        if (ib_indFvL(i) .eq. 0.) call Driver_abortFlash(" Lagrangian CalcForce: Marker V with no stencil !!")


        hl  = (ib_dsxu2(i)/ib_indFuL(i)+ib_dsyv2(i)/ib_indFvL(i))/(ib_alphax+ib_alphay)

        volL = dsb*hl

        FxL = FxL + volL*ib_FuL2(i)/ib_indFuL(i)

        FyL = FyL + volL*ib_FvL2(i)/ib_indFvL(i)

        MomzL = MomzL - volL*ib_FuL2(i)/ib_indFuL(i)*(ib_yb(i)-ib_BodyCen(JAXIS,ibd)) + &
                        volL*ib_FvL2(i)/ib_indFvL(i)*(ib_xb(i)-ib_BodyCen(IAXIS,ibd))   
   


     enddo
  enddo


  if (gr_meshMe .eq. MASTER_PE)  then
     write(*,*) ' '
     write(*,'(A35,3g18.12)') ' Total Lagrange Fx , Fy & Mz=',FxL, FyL, MomzL
     write(*,'(A35,3g18.12)') ' Lagrange Drag , Lift & Moment ceoff=',2.*FxL, 2.*FyL,2.*MomzL
  endif


  if (firstflg) then

     if (gr_meshMe .eq. MASTER_PE)  then
     open(unit=113,file='./IOData/Cyl_forcesLagr.res',form='formatted',status='replace')
     write(113,'(4g16.8)')dr_simTime,FxL,FyL,MomzL
     close(113)
     endif

!     firstflg = .false.

  else

     if (gr_meshMe .eq. MASTER_PE)  then
     open(unit=113,file='./IOData/Cyl_forcesLagr.res',form='formatted',status='old',position='append')
     write(113,'(4g16.8)')dr_simTime,FxL,FyL,MomzL
     close(113)
     endif
  endif  



  ! Compute forces integrating surface forces:
  ! Distributed Forces:
  call ib_forces2D(blockCount, blockList, dt)


  ! Now do the integration:
  Qx=0.; Qy=0.; Qtita=0.;
  kfact2 = 1./6.
  tita = 0.
  do ibd=1,ib_nbd

     ! Use body elements:
     ind1e = ib_ABODY(ibd) % elb + 1
     ind2e = ib_ABODY(ibd) % elb + ib_ABODY(ibd) % emb

     do i = ind1e,ind2e

        ind1=ib_AELEM(2,i)
        ind2=ib_AELEM(3,i)

        !Lseg = sqrt((ib_xb(ind2)-ib_xb(ind1))**2+(ib_yb(ind2)-ib_yb(ind1))**2);
        Lseg = ib_sb(i+1)-ib_sb(i)   

 
        t = Lseg**(-1)*(/ ib_xb(ind2)-ib_xb(ind1), ib_yb(ind2)-ib_yb(ind1), 0./);
    
        Qxs = 0.5*Lseg*(t(1)*(ib_Fd(1,ind1)+ib_Fd(1,ind2)) - &
                        t(2)*(ib_Fd(2,ind1)+ib_Fd(2,ind2)));
        Qys = 0.5*Lseg*(t(2)*(ib_Fd(1,ind1)+ib_Fd(1,ind2)) + &
                        t(1)*(ib_Fd(2,ind1)+ib_Fd(2,ind2)));
    
        tau1tau2_1 = 2.*ib_Fd(1,ind1) +    ib_Fd(1,ind2);
        tau1tau2_2 =    ib_Fd(1,ind1) + 2.*ib_Fd(1,ind2);
        
        pr1pr2_1 = 2.*ib_Fd(2,ind1) +    ib_Fd(2,ind2);
        pr1pr2_2 =    ib_Fd(2,ind1) + 2.*ib_Fd(2,ind2);
    
        kfact = -kfact2*Lseg;
    
        t1mt2_1 = t(1)*tau1tau2_1-t(2)*pr1pr2_1;
        t2pt1_1 = t(2)*tau1tau2_1+t(1)*pr1pr2_1;
        t1mt2_2 = t(1)*tau1tau2_2-t(2)*pr1pr2_2;
        t2pt1_2 = t(2)*tau1tau2_2+t(1)*pr1pr2_2;
    
        dQtitasdrp1x=kfact*(sin(tita)*t1mt2_1-cos(tita)*t2pt1_1);
        dQtitasdrp1y=kfact*(sin(tita)*t2pt1_1+cos(tita)*t1mt2_1);
        dQtitasdrp2x=kfact*(sin(tita)*t1mt2_2-cos(tita)*t2pt1_2);
        dQtitasdrp2y=kfact*(sin(tita)*t2pt1_2+cos(tita)*t1mt2_2);
    
    
        Qtitas = dQtitasdrp1x*ib_xbus(ind1) + dQtitasdrp1y*ib_ybus(ind1) + dQtitasdrp2x*ib_xbus(ind2) + dQtitasdrp2y*ib_ybus(ind2);
        ! HERE the x and y coordinates of points 1 and 2 should be
        ! the ones of the deformed body in local coordinates. In this
        ! Part they are the same of the undeformed (rigid).
    
        Qx = Qx + Qxs;
        Qy = Qy + Qys;
        Qtita  = Qtita  + Qtitas;
     enddo
  enddo


  if (gr_meshMe .eq. MASTER_PE)  then
     write(*,*) ' '
     write(*,'(A35,3g14.6)') ' Total Integrated Fx , Fy & Mz=',Qx, Qy, Qtita
     write(*,'(A35,3g14.6)') ' Integrated Drag , Lift & Moment ceoff=',2.*Qx, 2.*Qy,2.*Qtita
  endif


  if (firstflg) then

     if (gr_meshMe .eq. MASTER_PE)  then
     open(unit=113,file='./IOData/Cyl_forcesInteg.res',form='formatted',status='replace')
     write(113,'(4g16.8)')dr_simTime,Qx,Qy,Qtita
     close(113)
     endif

     firstflg = .false.

  else
     if (gr_meshMe .eq. MASTER_PE)  then
     open(unit=113,file='./IOData/Cyl_forcesInteg.res',form='formatted',status='old',position='append')
     write(113,'(4g16.8)')dr_simTime,Qx,Qy,Qtita
     close(113)
     endif
  endif  


end SUBROUTINE ib_CalcForce

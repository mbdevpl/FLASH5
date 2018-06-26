
!!             integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount)
!!             real(IN)    :: dt)
  ! This routine computes the immersed boundary forcing for a 
  ! given time step.  A diffuse interface method is employed.


SUBROUTINE ib_forcing(blockCount, blockList, dt)

  use Grid_data, ONLY: gr_meshMe

  use ImBound_data

  use ib_interface, ONLY : ib_stencils,ib_interpLpoints,ib_extrapEpoints

  use Grid_interface, ONLY : Grid_getDeltas,          &
                             Grid_getBlkCenterCoords, &
                             Grid_getBlkPhysicalSize, &
                             Grid_getBlkPtr,          &
                             Grid_releaseBlkPtr

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
  real, INTENT(IN) :: dt
  !! -----------------------------------------------------

  integer, parameter :: ng = NGUARD
!!$  integer, parameter :: nxb = NXB
!!$  integer, parameter :: nyb = NYB
!!$  integer, parameter :: nzb = NZB
  
  integer, parameter :: nxc = NXB + NGUARD + 1
  integer, parameter :: nyc = NYB + NGUARD + 1
#if NDIM == 3
  integer, parameter :: nzc = NZB + NGUARD + 1
#else
  integer, parameter :: nzc = 1
#endif

  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData

  ! counters 
  integer :: i,j,k,lb,blockID
  ! Number of points found in the block for u, v and w grids
  integer :: npu,npv,npw 
  integer, dimension(ib_nmaxa) :: lpindexu,lpindexv

  real :: del(MDIM),coord(MDIM),bsize(MDIM)
  real :: invdt
 
  integer :: ierr

  ! Uncomment to set position and velocities to 0
  !spcoord = 0.0
  !ubdd = 0.0
  !vbdd = 0.0
  !ubd0 = 0.0
  !vbd0 = 0.0

  invdt = 1./dt

  ! Indexes of stencils
  ib_ielemu = 0
  ib_jelemu = 0
  ib_ielemv = 0
  ib_jelemv = 0

#if NDIM == 3
  ib_kelemu = 0
  ib_kelemv = 0
  ib_ielemw = 0
  ib_jelemw = 0
  ib_kelemw = 0
#endif

  ib_dsxu=0.; ib_dsyu=0.
  ib_dsxv=0.; ib_dsyv=0.
  ib_dsxu2=0.; ib_dsyv2=0.
  ib_indFuL(1:ib_nnoda) = 0.
  ib_indFvL(1:ib_nnoda) = 0.
  ib_FuL2(1:ib_nnoda) = 0.
  ib_FvL2(1:ib_nnoda) = 0.
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

     ! Get block center coords:
     call Grid_getBlkCenterCoords(blockId,coord)

     ! Get block size:
     call Grid_getBlkPhysicalSize(blockID,bsize)


#if NDIM == 2
     ! Obtain Stencils for interpolation to Lagrangian points:
     ! U velocities
     call ib_stencils(ng,nxc+1,nyc,ib_nbd,ib_xb,ib_yb,ib_sb,ib_dsxu,ib_dsyu,            &
                  ib_ielemu,ib_jelemu,flaguo(:,:,:,blockID),flagui(:,:,:,blockID),      &
                  ib_nxL,ib_nyL,blockID,IAXIS,npu,lpindexu,                             &
                  GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC,del,coord,bsize)
                                     
     ! V velocities
     call ib_stencils(ng,nxc,nyc+1,ib_nbd,ib_xb,ib_yb,ib_sb,ib_dsxv,ib_dsyv,            &
                  ib_ielemv,ib_jelemv,flagvo(:,:,:,blockID),flagvi(:,:,:,blockID),      &
                  ib_nxL,ib_nyL,blockID,JAXIS,npv,lpindexv,                             &
                  GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC,del,coord,bsize)
                                     
#elif NDIM == 3
     if (gr_meshMe .eq. MASTER_PE) write(*,*) 'Stencils not 3D yet!'
     call Driver_abortFlash("ib_forcing: Not 3D ready yet!")
#endif               

     if(npu .gt. 0) then

        !write(*,*) 'Into interp U .., npu=',npu


        ! Interpolation of the values of velocity to Lagrangian points:
#if NDIM == 2
        call ib_interpLpoints(ib_nbd,ib_xb,ib_yb,ib_sb,ib_dsxu,ib_dsyu,                 &
                              ib_ielemu,ib_jelemu,ib_phileu,ib_UL,                      &
                     facexData(VELC_FACE_VAR,:,:,:),blockID,IAXIS,npu,lpindexu,         &
                     GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC,ng,nxc+1,nyc,del,coord,bsize)
                                        
#endif

        ib_FuL(1:ib_nnoda) = invdt*(ib_ubd(1:ib_nnoda) - ib_UL(1:ib_nnoda))

        ib_indFuL(lpindexu(1:npu)) = ib_indFuL(lpindexu(1:npu)) + 1.

        ib_dsxu2(lpindexu(1:npu)) = ib_dsxu2(lpindexu(1:npu)) + ib_dsxu(lpindexu(1:npu))       
              
        ib_FuL2(lpindexu(1:npu)) = ib_FuL2(lpindexu(1:npu)) + ib_FuL(lpindexu(1:npu))

        


        ! Extrapolation of Forces to the Eulerian Points (Check what 
        ! happens to the positions accelerations and velocities) and 
        ! sum to ustar and wstar:
#if NDIM == 2
        call ib_extrapEpoints(ib_nbd,ib_xb,ib_yb,ib_sb,                   &
                     ib_ielemu,ib_jelemu,ib_phileu,ib_FuL,                &
                     facexData(FORC_FACE_VAR,:,:,:),blockID,npu,lpindexu, &
                     GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC,del,coord,bsize) 
#endif                                  

        ! Sum to ustar: u** = u* + dt*fu
        facexData(VELC_FACE_VAR,:,:,:) = facexData(VELC_FACE_VAR,:,:,:) + &
                                         dt*facexData(FORC_FACE_VAR,:,:,:)
              
     endif 


     if(npv .gt. 0) then

#if NDIM == 2
        call ib_interpLpoints(ib_nbd,ib_xb,ib_yb,ib_sb,ib_dsxv,ib_dsyv,                 &
                     ib_ielemv,ib_jelemv,ib_philev,ib_VL,                               &
                     faceyData(VELC_FACE_VAR,:,:,:),blockID,JAXIS,npv,lpindexv,         &
                     GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC,ng,nxc,nyc+1,del,coord,bsize)
                                        
#endif

        ib_FvL(1:ib_nnoda) = invdt*(ib_vbd(1:ib_nnoda) - ib_VL(1:ib_nnoda))

        ib_indFvL(lpindexv(1:npv)) = ib_indFvL(lpindexv(1:npv)) + 1.        

        ib_dsyv2(lpindexv(1:npv)) = ib_dsyv2(lpindexv(1:npv)) + ib_dsyv(lpindexv(1:npv))
              
        ib_FvL2(lpindexv(1:npv)) = ib_FvL2(lpindexv(1:npv)) + ib_FvL(lpindexv(1:npv))


#if NDIM == 2
        call ib_extrapEpoints(ib_nbd,ib_xb,ib_yb,ib_sb,                   &
                     ib_ielemv,ib_jelemv,ib_philev,ib_FvL,                &
                     faceyData(FORC_FACE_VAR,:,:,:),blockID,npv,lpindexv, &
                     GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC,del,coord,bsize) 
#endif                                  

        ! Sum to vstar: v** = v* + dt*fv
        faceyData(VELC_FACE_VAR,:,:,:) = faceyData(VELC_FACE_VAR,:,:,:) + &
                                         dt*faceyData(FORC_FACE_VAR,:,:,:)


     endif

     ! Set the value of BFLAGS(1,lb) and BFLAGS(2,lb)
     ! To npu and npv to individualize which blocks have
     ! some Lagrangean points on inside.
     ib_bflags(IAXIS,blockID) = npu
     ib_bflags(JAXIS,blockID) = npv

               
     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

    
  enddo


  ! Sum all the index arrays and forces from all processors to avoid duplication of force in some Lagrangian markers.
  ib_FuL(1:ib_nnoda) = ib_indFuL(1:ib_nnoda)
  call MPI_Allreduce(ib_FuL(1:ib_nnoda),ib_indFuL(1:ib_nnoda), ib_nnoda, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)

  ib_FuL(1:ib_nnoda) = ib_FuL2(1:ib_nnoda)
  call MPI_Allreduce(ib_FuL(1:ib_nnoda),ib_FuL2(1:ib_nnoda), ib_nnoda, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)
  
  ib_FvL(1:ib_nnoda) = ib_indFvL(1:ib_nnoda)
  call MPI_Allreduce(ib_FvL(1:ib_nnoda),ib_indFvL(1:ib_nnoda), ib_nnoda, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)

  ib_FvL(1:ib_nnoda) = ib_FvL2(1:ib_nnoda)
  call MPI_Allreduce(ib_FvL(1:ib_nnoda),ib_FvL2(1:ib_nnoda), ib_nnoda, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)


  ib_FuL(1:ib_nnoda) = ib_dsxu2(1:ib_nnoda)
  call MPI_Allreduce(ib_FuL(1:ib_nnoda),ib_dsxu2(1:ib_nnoda), ib_nnoda, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)


!  ib_FuL(1:ib_nnoda) = ib_dsyu(1:ib_nnoda)
!  call MPI_Allreduce(ib_FuL(1:ib_nnoda),ib_dsyu(1:ib_nnoda), ib_nnoda, FLASH_REAL, &
!       MPI_SUM, MPI_COMM_WORLD, ierr)


!  ib_FvL(1:ib_nnoda) = ib_dsxv(1:ib_nnoda)
!  call MPI_Allreduce(ib_FvL(1:ib_nnoda),ib_dsxv(1:ib_nnoda), ib_nnoda, FLASH_REAL, &
!       MPI_SUM, MPI_COMM_WORLD, ierr)

  ib_FvL(1:ib_nnoda) = ib_dsyv2(1:ib_nnoda)
  call MPI_Allreduce(ib_FvL(1:ib_nnoda),ib_dsyv2(1:ib_nnoda), ib_nnoda, FLASH_REAL, &
       MPI_SUM, MPI_COMM_WORLD, ierr)






END SUBROUTINE ib_forcing


! #################################################################################
! Subroutine ib_forces
! ---------------------------------------------------------------------------------
subroutine ib_forces2D(blockCount, blockList, dt)

  ! This routine computes the immersed boundary forcing for a 
  ! given time step.  A diffuse interface method is employed.

  use Grid_interface, ONLY : Grid_getDeltas,          &
                             Grid_getBlkPtr,          &
                             Grid_releaseBlkPtr

  use ImBound_data

  use ib_interface, ONLY : ib_calcdistforce2D

  use IncompNS_data, ONLY : ins_invRe

  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
  real, INTENT(IN) :: dt
  !! -----------------------------------------------------
      
  integer, parameter :: nxi = NGUARD + NXB
  integer, parameter :: nyj = NGUARD + NYB      
  integer, parameter :: nxc = NGUARD + NXB + 1
  integer, parameter :: nyc = NGUARD + NYB + 1

  real :: nu

  integer :: i,j,ierr
  integer :: lb, blockID
  real :: dx,dy
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC+1) :: vort

  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,solnData

  real :: del(MDIM)

  ib_Fd(CONSTANT_ONE,:) = 0. ! Tangential Force on Markers 
  ib_Fd(CONSTANT_TWO,:) = 0. ! Pressure on Markers

  nu = ins_invRe

  !write(*,*) 'Into Calcdistforce ..'

  do lb = 1,blockCount
         
     blockID = blockList(lb)
     
     if ((ib_bflags(IAXIS,blockID) .ne. 0) .or. &
         (ib_bflags(JAXIS,blockID) .ne. 0)) then

        ! Get face data (velocities):
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)

        ! Get dx,dy
        call Grid_getDeltas(blockID,del)
                 
        ! Initialize vorticity:
        vort(:,:) = 0.
         
        ! Calculate the vorticity of the block:
        dx = del(IAXIS) 
        dy = del(JAXIS) 
               
        ! Works only for two or more layers of guardcells
        vort(CONSTANT_TWO:GRID_IHI_GC,CONSTANT_TWO:GRID_JHI_GC) = &
         ( faceyData(VELC_FACE_VAR,CONSTANT_TWO:GRID_IHI_GC,CONSTANT_TWO:GRID_JHI_GC,CONSTANT_ONE) -           &
           faceyData(VELC_FACE_VAR,CONSTANT_TWO-1:GRID_IHI_GC-1,CONSTANT_TWO:GRID_JHI_GC,CONSTANT_ONE) )/dx -  &
         ( facexData(VELC_FACE_VAR,CONSTANT_TWO:GRID_IHI_GC,CONSTANT_TWO:GRID_JHI_GC,CONSTANT_ONE) -           &
           facexData(VELC_FACE_VAR,CONSTANT_TWO:GRID_IHI_GC,CONSTANT_TWO-1:GRID_JHI_GC-1,CONSTANT_ONE) )/dy

        vort(CONSTANT_ONE,CONSTANT_TWO:GRID_JHI_GC)  = vort(CONSTANT_TWO,CONSTANT_TWO:GRID_JHI_GC)
        vort(GRID_IHI_GC+1,CONSTANT_TWO:GRID_JHI_GC) = vort(GRID_IHI_GC,CONSTANT_TWO:GRID_JHI_GC)
        vort(CONSTANT_TWO:GRID_IHI_GC,CONSTANT_ONE)  = vort(CONSTANT_TWO:GRID_IHI_GC,CONSTANT_TWO)
        vort(CONSTANT_TWO:GRID_IHI_GC,GRID_JHI_GC+1) = vort(CONSTANT_TWO:GRID_IHI_GC,GRID_JHI_GC)

!        vort(NGUARD:nxc+1,NGUARD:nyc+1) =  &
!         ( faceyData(VELC_FACE_VAR,NGUARD:nxc+1,NGUARD:nyc+1,CONSTANT_ONE)         & 
!         - faceyData(VELC_FACE_VAR,NGUARD-1:nxc,NGUARD:nyc+1,CONSTANT_ONE) )/dx -  &
!         ( facexData(VELC_FACE_VAR,NGUARD:nxc+1,NGUARD:nyc+1,CONSTANT_ONE) -       &
!           facexData(VELC_FACE_VAR,NGUARD:nxc+1,NGUARD-1:nyc,CONSTANT_ONE) )/dy 

        ! Obtain Distributed Tangential stresses in
        ! marker points:
        call ib_calcdistforce2D(nu,NGUARD,nxc,nyc,ib_nbd,ib_nnoda,                   &
                                ib_xb,ib_yb,ib_sb,ib_nxL,ib_nyL,ib_xbe,ib_ybe,       &
                                ib_stencil,ib_Fd(CONSTANT_ONE,:),vort,CONSTANT_ZERO, &
                                ib_ubdd,ib_vbdd,blockID,nxc+NGUARD,nyc+NGUARD)   


        ! Obtain Pressures in marker points:
        !                                 nxi,nyj
        call ib_calcdistforce2D(nu,NGUARD,nxi,nyj,ib_nbd,ib_nnoda,                                                &
                                ib_xb,ib_yb,ib_sb,ib_nxL,ib_nyL,ib_xbe,ib_ybe,                                    &
                                ib_stencil,ib_Fd(CONSTANT_TWO,:),                                                 &
                                solnData(PRES_VAR,GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,CONSTANT_ONE),  &
                                CONSTANT_ONE,ib_ubdd,ib_vbdd,blockID,nxi+NGUARD,nyj+NGUARD)



     ! Release face data (velocities):
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

     endif
            
  enddo

  ! Allreduce sum:
  ib_Fdaux(CONSTANT_ONE:CONSTANT_TWO,1:ib_nnoda) = ib_Fd(CONSTANT_ONE:CONSTANT_TWO,1:ib_nnoda)

  call MPI_Allreduce(ib_Fdaux(CONSTANT_ONE:CONSTANT_TWO,1:ib_nnoda),                      &
                     ib_Fd(CONSTANT_ONE:CONSTANT_TWO,1:ib_nnoda), CONSTANT_TWO*ib_nnoda,  &
                     FLASH_REAL,MPI_SUM, MPI_COMM_WORLD, ierr)


!  ! Write Files:
!  open(unit=55,file='./IOData/Pressure.dat',status='unknown',form='formatted')
!  do i = 1,ib_nnoda
!    write(55,'(2g24.16)')ib_sb(i),ib_Fd(2,i)
!  enddo
!  close(55)

!  open(unit=55,file='./IOData/Viscous.dat',status='unknown',form='formatted')
!  do i = 1,ib_nnoda
!    write(55,'(2g24.16)')ib_sb(i),ib_Fd(1,i)
!  enddo 
!  close(55)      

  return
   

END subroutine ib_forces2D



! Subroutine calcdistforce:
!
! --------------------------------------------------------------

subroutine ib_calcdistforce2D(nu,ng,nx,ny,nbd,nnoda,xb,yb,sb, &
                              nxp,nyp,xbe,ybe,                &
                              stencil,zL,VarO,presflag,       &
                              ubdd,vbdd,lb,nx1,ny1)

  use Grid_interface, ONLY : Grid_getDeltas,          &
                             Grid_getBlkCenterCoords, &
                             Grid_getBlkPhysicalSize
 
  use  ImBound_data, ONLY   : ib_alphax,ib_alphay,ib_ABODY, ib_nmaxa,ib_interp,ib_npol

  use ib_interface, ONLY : ib_buildABLan

  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  ! --- ARGUMENTS -----------------------------------------
  integer, INTENT(IN) :: ng,lb
  integer, INTENT(IN) :: nx,ny,nbd,nnoda,stencil,presflag,nx1,ny1
  real  , INTENT(IN)  :: nu,xb(ib_nmaxa),yb(ib_nmaxa),sb(ib_nmaxa)
  real  , INTENT(IN)  ::  nxp(ib_nmaxa),nyp(ib_nmaxa)
  real  , INTENT(INOUT)  ::  xbe(ib_nmaxa),ybe(ib_nmaxa)
  real  , INTENT(INOUT)  ::  zL(ib_nmaxa)
  real  , INTENT(IN) :: ubdd(ib_nmaxa),vbdd(ib_nmaxa)
  real  ,INTENT(IN)   :: VarO(nx1,ny1)
  ! -------------------------------------------------------


  ! Local Variables:
  integer i,j,ii,ibd,ind1,ind2
  integer ipos, iposp1, jpos, jposp1
  real ::  dxloc, dyloc, factor
  real  :: coord(MDIM),bsize(MDIM), Del(MDIM)

  integer :: npoints, indx1, indx2, indy1, indy2
  real ::  cublim(4),auxx,auxy,distx,disty,xp,yp
  integer :: indbx1(nbd),indbx2(nbd),indby1(nbd),indby2(nbd)
  integer :: icpoint,jcpoint,flagio
  real ::  alphax,alphay


  real ::  dsx,dsy,h
  integer :: ielem(stencil),jelem(stencil)

  real ::  d,zp,dpdn

  real :: A(ib_npol,ib_npol), B(ib_npol,stencil)
  real :: p(ib_npol), phi(stencil), gamma(ib_npol), indx(ib_npol)

      
  integer :: bndx1,bndx2,bndy1,bndy2,count,buildflag
  real :: xlower,xupper,ylower,yupper,dx,dy,dxaux,dyaux,xi,yj
  real :: xmin,xmax,ymin,ymax,eps

  real :: x(stencil),y(stencil)
 
  real :: dVx,dVy
  
  !real :: aux(ib_nmaxa)
  !real :: VarOngng,VarOnxng,VarOnxny,VarOngny


  ! Get dx,dy
  call Grid_getDeltas(lb,del)
  call Grid_getBlkCenterCoords(lb,coord)
  call Grid_getBlkPhysicalSize(lb,bsize)

  dx=Del(IAXIS)
  dy=Del(JAXIS)

  dsx = ib_alphax*dx
  dsy = ib_alphay*dy

  h = 0.5*(dsx+dsy)  

  eps = 1E-2*MIN(dx,dy)

  ! Find block boundaries
  xlower = coord(IAXIS) - bsize(IAXIS)/2.0 
  xupper = coord(IAXIS) + bsize(IAXIS)/2.0 
  ylower = coord(JAXIS) - bsize(JAXIS)/2.0 
  yupper = coord(JAXIS) + bsize(JAXIS)/2.0 
      


  if(presflag .eq. CONSTANT_ZERO) then        ! Vorticity
     dxaux = 0.0
     dyaux = 0.0
     
     bndx1 = CONSTANT_TWO !ng + 1
     bndx2 = GRID_IHI_GC  !nx        
     bndy1 = CONSTANT_TWO !ng + 1
     bndy2 = GRID_JHI_GC  !ny        
  elseif(presflag .eq. CONSTANT_ONE) then    ! Pressure
     dxaux = 0.5*dx
     dyaux = 0.5*dy
     
     bndx1 = CONSTANT_TWO  !ng 
     bndx2 = GRID_IHI_GC-1 !nx + 1      
     bndy1 = CONSTANT_TWO  !ng 
     bndy2 = GRID_JHI_GC-1 !ny + 1      
  endif

  count = 0

  ! Obtain positions of external pressure points:
  do ibd=1,nbd  

     ! Find bounding box of body
     ind1 = ib_ABODY(ibd) % lb + 1
     ind2 = ib_ABODY(ibd) % lb + ib_ABODY(ibd) % mb

     xmin = MINVAL( xb(ind1:ind2) )
     xmax = MAXVAL( xb(ind1:ind2) )
     ymin = MINVAL( yb(ind1:ind2) )
     ymax = MAXVAL( yb(ind1:ind2) )    



     ! Loop through Eulerian points box:
     do ii = ind1 , ind2

        xp = xb(ii);
        yp = yb(ii);

        if( xp .ge. xlower .and. xp .lt. (xupper) .and. &
            yp .ge. ylower .and. yp .lt. (yupper) ) then
               

           xbe(ii) = xp + nxp(ii)*h
           ybe(ii) = yp + nyp(ii)*h


           ! Obtain stencil for External Pressure point:
           xp = xbe(ii)
           yp = ybe(ii)

           ! Find closest point:
           distx = 10.e10;
           disty = 10.e10;
           icpoint = 0
           jcpoint = 0
           do i = bndx1,bndx2
              xi = coord(IAXIS) - 0.5*bsize(IAXIS) + &
                   real(i - ng - 1)*dx + dxaux
              auxx = abs(xp - xi);
              if (auxx+eps .le. distx) then
                 icpoint = i;
                 distx = auxx;
              endif
           enddo
    
           do j = bndy1,bndy2
              yj = coord(JAXIS) - 0.5*bsize(JAXIS) + &
                   real(j - ng - 1)*dy + dyaux
              auxy = abs(yp - yj);
              if (auxy+eps .le. disty) then
                 jcpoint = j;
                 disty = auxy;
              endif
           enddo


           ! MODIFY CLOSEST POINT REGARDING THE EXTERNAL
           ! POINT POSITION RESPECT TO THE BLOCK CORNERS:
           !if (.false.) then !((presflag .eq. CONSTANT_ZERO) .and.                &
              !((icpoint .eq. ng+1) .or. (icpoint .eq. nx)) .and.  &
              !((jcpoint .eq. ng+1) .or. (jcpoint .eq. ny))) then

              !zp = VarO(icpoint,jcpoint);
               
           !elseif ((presflag .eq. CONSTANT_ONE) .and.             &
              !((icpoint .eq. ng) .or. (icpoint .eq. nx+1)) .and.  &
              !((jcpoint .eq. ng) .or. (jcpoint .eq. ny+1))) then

              !zp = VarO(icpoint,jcpoint);
              
           !else

              ! Build stencil structure:
              ! Center North South East West:
              !if (stencil .eq. 5) then
              ielem = (/ icpoint,icpoint,icpoint,icpoint+1,icpoint-1 /)
              jelem = (/ jcpoint,jcpoint+1,jcpoint-1,jcpoint,jcpoint /)
              !elseif (stencil .eq. 9) then
              !ielem = (/ icpoint,icpoint,icpoint,icpoint+1,icpoint-1,icpoint-1,icpoint-1,icpoint+1,icpoint+1 /)
              !jelem = (/ jcpoint,jcpoint+1,jcpoint-1,jcpoint,jcpoint,jcpoint-1,jcpoint+1,jcpoint-1,jcpoint+1 /) 
              !endif

              x = coord(IAXIS) - 0.5*bsize(IAXIS) + &
                  real(ielem(1:stencil) - ng - 1)*dx + dxaux
              y = coord(JAXIS) - 0.5*bsize(JAXIS) + &
                  real(jelem(1:stencil) - ng - 1)*dy + dyaux


              ! Build A and B matrices:
              call ib_buildABLan(stencil,ib_npol,1.*dsx,1.*dsy,xp,yp, &
                                 x,y,ib_interp,A,B,buildflag); 


              if (buildflag == CONSTANT_ONE) then
                 write(*,*) 'ZONE =', lb
                 write(*,*) 'dsx,dsy =',dsx,dsy
                 write(*,*) 'xp,yp =',xp,yp
                 write(*,*) 'icpoint,jcpoint =',icpoint,jcpoint
                 do i = 1,stencil
                    write(*,*) ielem(i),jelem(i),x(i),y(i)
                 enddo
              endif
              

              ! Obtain gamma coefficients:
              ! Solve for systems coefficients:
              if (ib_interp == CONSTANT_ONE) then
                 p(1) = 1.; p(2) = xp; p(3) = yp;     
              endif

              gamma = p
            
              call ludcmp(A,ib_npol,ib_npol,indx,d)
              call lubksb(A,ib_npol,ib_npol,indx,gamma)

              ! Obtain Shape functions:
              phi =0.
              do i = 1 , stencil
                 do j = 1, ib_npol
                    phi(i) = phi(i) + gamma(j)*B(j,i)
                 enddo
              enddo
              
            
              ! Value of the function in xp and yp:
              zp = 0.;
              do i = 1 , stencil      
                 zp = zp + phi(i)*VarO(ielem(i),jelem(i));   
              enddo

           !endif

           if (presflag .eq. CONSTANT_ONE) then         ! Pressure
              dpdn = -(ubdd(ii)*nxp(ii) + vbdd(ii)*nyp(ii));
              zL(ii) = zp - dpdn*h;
           else                              ! Tangent stress
              zL(ii) = nu*zp;
           endif

        endif

     enddo
  enddo
  
  return

End Subroutine ib_calcdistforce2D


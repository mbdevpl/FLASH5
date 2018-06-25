
! This routine fixes inconsistent fluxes across block interfaces 
! of different refinement levels as defined by the normal 
! velocity component.
!
! Adapted to Flash 3, Marcos Vanella, September 2007.

!-------------------------------------------------------------------
!- kpd - These routines are only used for ParaMesh/AMR runs, as they
!           fix inconsitent fluxes across non-conformal boundaries.
!
!      - Subroutine ins_fluxfix is used for predictor velocities
!      - Subroutine ins_fluxfix_p is used for pressure GRADIENT 
!           after the pressure Poisson is solved
!-------------------------------------------------------------------

  SUBROUTINE ins_fluxfix(ng,nxc,nyc,nzc,nxi,nyj,nzk,blockCount,&
                         blockList)

  
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes
  use IncompNS_data, ONLY : ins_meshMe

  implicit none
  
#include "constants.h"
#include "Flash.h"  
#include "IncompNS.h"

  integer, intent(IN) :: ng,nxc,nyc,nzc,nxi,nyj,nzk, &
                         blockCount     
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList


  integer, dimension(MDIM) :: datasize
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData

  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_u
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_v
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_w

  integer :: sx,sy,sz,ex,ey,ez
  integer :: lb,blockID,level
  real :: del(MDIM)
  
  flxint_u = 0.
  flxint_v = 0.
  flxint_w = 0.

  level = -1

  ! Dump fluxes into flux variables:   
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take cell depth to be 1:
     del(DIR_Z) = 1. 
#endif

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   

     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(1)-NGUARD
     ey = dataSize(2)-NGUARD*K2D
     ez = dataSize(3)-NGUARD*K3D

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif
     ! copy interface velocities to flux variables
     ! X direction:
     !- kpd - U-velocities that lie on left and right boundaries
     flxint_u(VELC_FLUX,sx,sy:ey,sz:ez) = &
     facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez)*del(DIR_Y)*del(DIR_Z)
     flxint_u(VELC_FLUX,ex+1,sy:ey,sz:ez) = &
     facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez)*del(DIR_Y)*del(DIR_Z)

     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, IAXIS, &
                           flxint_u, dataSize)             ! dy*dz


     ! Y direction:
     flxint_v(VELC_FLUX,sx:ex,sy,sz:ez) = &
     faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez)*del(DIR_X)*del(DIR_Z)
     flxint_v(VELC_FLUX,sx:ex,ey+1,sz:ez) = &
     faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez)*del(DIR_X)*del(DIR_Z)

     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, JAXIS, &
                           flxint_v, dataSize)             ! dx*dz     


#if NDIM == 3
     ! Z direction:
     flxint_w(VELC_FLUX,sx:ex,sy:ey,sz) = &
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz)*del(DIR_X)*del(DIR_Y)
     flxint_w(VELC_FLUX,sx:ex,sy:ey,ez+1) = &
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1)*del(DIR_X)*del(DIR_Y)

     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, KAXIS, &
                           flxint_w, dataSize)             ! dx*dy 
#endif


     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo

!  write(*,*) 'Mype=',Ins_MeshMe, 'Before Grid Conserve fluxes'

  ! fix fluxes at refinement interfaces
  call Grid_conserveFluxes( ALLDIR, level)

!  write(*,*) 'Mype=',Ins_MeshMe, 'After Grid Conserve fluxes'

 
  ! Get fixed fluxes:
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take cell depth to be 1:
     del(DIR_Z) = 1. 
#endif

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   


     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(1)-NGUARD
     ey = dataSize(2)-NGUARD*K2D
     ez = dataSize(3)-NGUARD*K3D


     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     ! update interface velocities with flux variables
     ! X direction:
     call Grid_getFluxData(blockID, IAXIS, &
     flxint_u, dataSize)            

     facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez) = &
       flxint_u(VELC_FLUX,sx,sy:ey,sz:ez)/(del(DIR_Y)*del(DIR_Z))  ! dy*dz
     facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez) = &
       flxint_u(VELC_FLUX,ex+1,sy:ey,sz:ez)/(del(DIR_Y)*del(DIR_Z))


     ! Y direction:
     call Grid_getFluxData(blockID, JAXIS, &
     flxint_v, dataSize)                 

     faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez) = &
       flxint_v(VELC_FLUX,sx:ex,sy,sz:ez)/(del(DIR_X)*del(DIR_Z))  ! dx*dz
     faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez) = &
       flxint_v(VELC_FLUX,sx:ex,ey+1,sz:ez)/(del(DIR_X)*del(DIR_Z))

#if NDIM == 3
     ! Z direction:
     call Grid_getFluxData(blockID, KAXIS, &
     flxint_w, dataSize)              

     facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz) = &
       flxint_w(VELC_FLUX,sx:ex,sy:ey,sz)/(del(DIR_X)*del(DIR_Y))  ! dx*dy
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1) = &
       flxint_w(VELC_FLUX,sx:ex,sy:ey,ez+1)/(del(DIR_X)*del(DIR_Y))
#endif

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo


  return

END SUBROUTINE ins_fluxfix

!#############################################################################
!#############################################################################
!#############################################################################

  SUBROUTINE ins_fluxfixRho1(ng,nxc,nyc,nzc,nxi,nyj,nzk,blockCount,&
                         blockList)

  
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes
  use IncompNS_data, ONLY : ins_meshMe

  implicit none
  
#include "constants.h"
#include "Flash.h"  
#include "IncompNS.h"

  integer, intent(IN) :: ng,nxc,nyc,nzc,nxi,nyj,nzk, &
                         blockCount     
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList


  integer, dimension(MDIM) :: datasize
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData

  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_u
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_v
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_w

  integer :: sx,sy,sz,ex,ey,ez
  integer :: lb,blockID,level
  real :: del(MDIM)
  
  flxint_u = 0.
  flxint_v = 0.
  flxint_w = 0.

  level = -1

  ! Dump fluxes into flux variables:   
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take cell depth to be 1:
     del(DIR_Z) = 1. 
#endif

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   

     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(1)-NGUARD
     ey = dataSize(2)-NGUARD*K2D
     ez = dataSize(3)-NGUARD*K3D

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif
     ! copy interface velocities to flux variables
     ! X direction:
     !- kpd - U-velocities that lie on left and right boundaries
     flxint_u(RH1F_FLUX,sx,sy:ey,sz:ez) = &
     facexData(RH1F_FACE_VAR,sx,sy:ey,sz:ez)*del(DIR_Y)*del(DIR_Z)
     flxint_u(RH1F_FLUX,ex+1,sy:ey,sz:ez) = &
     facexData(RH1F_FACE_VAR,ex+1,sy:ey,sz:ez)*del(DIR_Y)*del(DIR_Z)

     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, IAXIS, &
                           flxint_u, dataSize)             ! dy*dz


     ! Y direction:
     flxint_v(RH1F_FLUX,sx:ex,sy,sz:ez) = &
     faceyData(RH1F_FACE_VAR,sx:ex,sy,sz:ez)*del(DIR_X)*del(DIR_Z)
     flxint_v(RH1F_FLUX,sx:ex,ey+1,sz:ez) = &
     faceyData(RH1F_FACE_VAR,sx:ex,ey+1,sz:ez)*del(DIR_X)*del(DIR_Z)

     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, JAXIS, &
                           flxint_v, dataSize)             ! dx*dz     


#if NDIM == 3
     ! Z direction:
     flxint_w(RH1F_FLUX,sx:ex,sy:ey,sz) = &
     facezData(RH1F_FACE_VAR,sx:ex,sy:ey,sz)*del(DIR_X)*del(DIR_Y)
     flxint_w(RH1F_FLUX,sx:ex,sy:ey,ez+1) = &
     facezData(RH1F_FACE_VAR,sx:ex,sy:ey,ez+1)*del(DIR_X)*del(DIR_Y)

     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, KAXIS, &
                           flxint_w, dataSize)             ! dx*dy 
#endif


     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo

!  write(*,*) 'Mype=',Ins_MeshMe, 'Before Grid Conserve fluxes'

  ! fix fluxes at refinement interfaces
  call Grid_conserveFluxes( ALLDIR, level)

!  write(*,*) 'Mype=',Ins_MeshMe, 'After Grid Conserve fluxes'

 
  ! Get fixed fluxes:
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take cell depth to be 1:
     del(DIR_Z) = 1. 
#endif

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   


     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(1)-NGUARD
     ey = dataSize(2)-NGUARD*K2D
     ez = dataSize(3)-NGUARD*K3D


     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     ! update interface velocities with flux variables
     ! X direction:
     call Grid_getFluxData(blockID, IAXIS, &
     flxint_u, dataSize)            

     facexData(RH1F_FACE_VAR,sx,sy:ey,sz:ez) = &
       flxint_u(RH1F_FLUX,sx,sy:ey,sz:ez)/(del(DIR_Y)*del(DIR_Z))  ! dy*dz
     facexData(RH1F_FACE_VAR,ex+1,sy:ey,sz:ez) = &
       flxint_u(RH1F_FLUX,ex+1,sy:ey,sz:ez)/(del(DIR_Y)*del(DIR_Z))


     ! Y direction:
     call Grid_getFluxData(blockID, JAXIS, &
     flxint_v, dataSize)                 

     faceyData(RH1F_FACE_VAR,sx:ex,sy,sz:ez) = &
       flxint_v(RH1F_FLUX,sx:ex,sy,sz:ez)/(del(DIR_X)*del(DIR_Z))  ! dx*dz
     faceyData(RH1F_FACE_VAR,sx:ex,ey+1,sz:ez) = &
       flxint_v(RH1F_FLUX,sx:ex,ey+1,sz:ez)/(del(DIR_X)*del(DIR_Z))

#if NDIM == 3
     ! Z direction:
     call Grid_getFluxData(blockID, KAXIS, &
     flxint_w, dataSize)              

     facezData(RH1F_FACE_VAR,sx:ex,sy:ey,sz) = &
       flxint_w(RH1F_FLUX,sx:ex,sy:ey,sz)/(del(DIR_X)*del(DIR_Y))  ! dx*dy
     facezData(RH1F_FACE_VAR,sx:ex,sy:ey,ez+1) = &
       flxint_w(RH1F_FLUX,sx:ex,sy:ey,ez+1)/(del(DIR_X)*del(DIR_Y))
#endif

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo


  return

END SUBROUTINE ins_fluxfixRho1

!#############################################################################
!#############################################################################
!#############################################################################
 
  SUBROUTINE ins_fluxfixRho2(ng,nxc,nyc,nzc,nxi,nyj,nzk,blockCount,&
                         blockList)

  
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes
  use IncompNS_data, ONLY : ins_meshMe

  implicit none
  
#include "constants.h"
#include "Flash.h"  
#include "IncompNS.h"

  integer, intent(IN) :: ng,nxc,nyc,nzc,nxi,nyj,nzk, &
                         blockCount     
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList


  integer, dimension(MDIM) :: datasize
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData

  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_u
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_v
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_w

  integer :: sx,sy,sz,ex,ey,ez
  integer :: lb,blockID,level
  real :: del(MDIM)
  
  flxint_u = 0.
  flxint_v = 0.
  flxint_w = 0.

  level = -1

  ! Dump fluxes into flux variables:   
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take cell depth to be 1:
     del(DIR_Z) = 1. 
#endif

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   

     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(1)-NGUARD
     ey = dataSize(2)-NGUARD*K2D
     ez = dataSize(3)-NGUARD*K3D

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif
     ! copy interface velocities to flux variables
     ! X direction:
     !- kpd - U-velocities that lie on left and right boundaries
     flxint_u(RH2F_FLUX,sx,sy:ey,sz:ez) = &
     facexData(RH2F_FACE_VAR,sx,sy:ey,sz:ez)*del(DIR_Y)*del(DIR_Z)
     flxint_u(RH2F_FLUX,ex+1,sy:ey,sz:ez) = &
     facexData(RH2F_FACE_VAR,ex+1,sy:ey,sz:ez)*del(DIR_Y)*del(DIR_Z)

     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, IAXIS, &
                           flxint_u, dataSize)             ! dy*dz


     ! Y direction:
     flxint_v(RH2F_FLUX,sx:ex,sy,sz:ez) = &
     faceyData(RH2F_FACE_VAR,sx:ex,sy,sz:ez)*del(DIR_X)*del(DIR_Z)
     flxint_v(RH2F_FLUX,sx:ex,ey+1,sz:ez) = &
     faceyData(RH2F_FACE_VAR,sx:ex,ey+1,sz:ez)*del(DIR_X)*del(DIR_Z)

     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, JAXIS, &
                           flxint_v, dataSize)             ! dx*dz     


#if NDIM == 3
     ! Z direction:
     flxint_w(RH2F_FLUX,sx:ex,sy:ey,sz) = &
     facezData(RH2F_FACE_VAR,sx:ex,sy:ey,sz)*del(DIR_X)*del(DIR_Y)
     flxint_w(RH2F_FLUX,sx:ex,sy:ey,ez+1) = &
     facezData(RH2F_FACE_VAR,sx:ex,sy:ey,ez+1)*del(DIR_X)*del(DIR_Y)

     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, KAXIS, &
                           flxint_w, dataSize)             ! dx*dy 
#endif


     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo

!  write(*,*) 'Mype=',Ins_MeshMe, 'Before Grid Conserve fluxes'

  ! fix fluxes at refinement interfaces
  call Grid_conserveFluxes( ALLDIR, level)

!  write(*,*) 'Mype=',Ins_MeshMe, 'After Grid Conserve fluxes'

 
  ! Get fixed fluxes:
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take cell depth to be 1:
     del(DIR_Z) = 1. 
#endif

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   


     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(1)-NGUARD
     ey = dataSize(2)-NGUARD*K2D
     ez = dataSize(3)-NGUARD*K3D


     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     ! update interface velocities with flux variables
     ! X direction:
     call Grid_getFluxData(blockID, IAXIS, &
     flxint_u, dataSize)            

     facexData(RH2F_FACE_VAR,sx,sy:ey,sz:ez) = &
       flxint_u(RH2F_FLUX,sx,sy:ey,sz:ez)/(del(DIR_Y)*del(DIR_Z))  ! dy*dz
     facexData(RH2F_FACE_VAR,ex+1,sy:ey,sz:ez) = &
       flxint_u(RH2F_FLUX,ex+1,sy:ey,sz:ez)/(del(DIR_Y)*del(DIR_Z))


     ! Y direction:
     call Grid_getFluxData(blockID, JAXIS, &
     flxint_v, dataSize)                 

     faceyData(RH2F_FACE_VAR,sx:ex,sy,sz:ez) = &
       flxint_v(RH2F_FLUX,sx:ex,sy,sz:ez)/(del(DIR_X)*del(DIR_Z))  ! dx*dz
     faceyData(RH2F_FACE_VAR,sx:ex,ey+1,sz:ez) = &
       flxint_v(RH2F_FLUX,sx:ex,ey+1,sz:ez)/(del(DIR_X)*del(DIR_Z))

#if NDIM == 3
     ! Z direction:
     call Grid_getFluxData(blockID, KAXIS, &
     flxint_w, dataSize)              

     facezData(RH2F_FACE_VAR,sx:ex,sy:ey,sz) = &
       flxint_w(RH2F_FLUX,sx:ex,sy:ey,sz)/(del(DIR_X)*del(DIR_Y))  ! dx*dy
     facezData(RH2F_FACE_VAR,sx:ex,sy:ey,ez+1) = &
       flxint_w(RH2F_FLUX,sx:ex,sy:ey,ez+1)/(del(DIR_X)*del(DIR_Y))
#endif

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo


  return

END SUBROUTINE ins_fluxfixRho2

!#############################################################################
!#############################################################################
!#############################################################################


!#############################################################################
!#############################################################################
!#############################################################################
! This routine fixes inconsistent fluxes across block interfaces 
! of different refinement levels as defined by the pressure
! gradient.
!
!- kpd - flxint_u is later used for the corrector velocity step
!#############################################################################


  SUBROUTINE ins_fluxfix_p(ng,nxc,nyc,nzc,nxi,nyj,nzk,pvar,&
                         blockCount,blockList)

  use Grid_interface, ONLY : Grid_getDeltas,         &                             
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes

  implicit none
  
#include "constants.h"
#include "Flash.h"  
#include "IncompNS.h"

  integer, intent(in) :: ng,nxc,nyc,nzc,nxi,nyj,nzk,pvar,&
                         blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList


  integer, dimension(MDIM) :: datasize
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  
  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData

  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_u
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_v
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_w

  integer :: sx,sy,sz,ex,ey,ez
  integer :: lb,blockID,level,i,j,k
  real :: del(MDIM),Mdens  

  ! Dump fluxes into flux variables:   
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take cell depth to be 1:
     del(DIR_Z) = 1. 
#endif

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   

     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(1)-NGUARD
     ey = dataSize(2)-NGUARD*K2D
     ez = dataSize(3)-NGUARD*K3D

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM ==3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

!--------------------------------------------
! Copy interface velocities to flux variables
!--------------------------------------------

  !-------------
  ! X direction:
  !-------------
     !- kpd - "i" is fixed at NGUARDCELLS+1 (dz =1 in 2-d)
     do k=sz,ez
     do j=sy,ey
     do i=sx,sx

     !- kpd - This is the inverse of the mixture density on the cell face...
     Mdens = ( facexData(RH1F_FACE_VAR,i,j,k) + facexData(RH2F_FACE_VAR,i,j,k) ) ! Inverse density in face.

     flxint_u(VELC_FLUX,i,j,k) =                                        &
                                 Mdens*(solnData(pvar,i  ,j,k)-         &
                                        solnData(pvar,i-1,j,k))         &
                                      *del(DIR_Y)*del(DIR_Z)/del(DIR_X)

     !print*,"FLUX IN",lb,i,j,solnData(pvar,i  ,j,k),solnData(pvar,i-1,j,k),flxint_u(VELC_FLUX,i,j,1)

     enddo
     enddo
     enddo

     !- kpd - "i" is fixed at NGUARDCELLS+NXB+1 (dz =1 in 2-d)
     do k=sz,ez
     do j=sy,ey
     do i=ex+1,ex+1

     !- kpd - This is the inverse of the mixture density on the cell face...
     Mdens = ( facexData(RH1F_FACE_VAR,i,j,k) + facexData(RH2F_FACE_VAR,i,j,k) ) ! Inverse density in face.

     flxint_u(VELC_FLUX,i,j,k) =                                      &
                                Mdens*(solnData(pvar,i  ,j,k)-        &
                                       solnData(pvar,i-1,j,k))        &
                                     *del(DIR_Y)*del(DIR_Z)/del(DIR_X)

     enddo
     enddo
     enddo

     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, IAXIS, &
                           flxint_u, dataSize)             ! dy*dz

  !-------------
  ! Y direction:
  !-------------
     do k=sz,ez
     do j=sy,sy
     do i=sx,ex

     !- kpd - This is the inverse of the mixture density on the cell face...
     Mdens = ( faceyData(RH1F_FACE_VAR,i,j,k) + faceyData(RH2F_FACE_VAR,i,j,k) ) ! Inverse density in face.

     flxint_v(VELC_FLUX,i,j,k) =                                      &
                                Mdens*(solnData(pvar,i,j,k)-          &
                                       solnData(pvar,i,j-1,k))        &
                                     *del(DIR_X)*del(DIR_Z)/del(DIR_Y)

     enddo
     enddo
     enddo

     do k=sz,ez
     do j=ey+1,ey+1
     do i=sx,ex

     !- kpd - This is the inverse of the mixture density on the cell face...
     Mdens = ( faceyData(RH1F_FACE_VAR,i,j,k) + faceyData(RH2F_FACE_VAR,i,j,k) ) ! Inverse density in face.

     flxint_v(VELC_FLUX,i,j,k) =                                      &
                                Mdens*(solnData(pvar,i,j,k)-          &
                                       solnData(pvar,i,j-1,k))        &
                                     *del(DIR_X)*del(DIR_Z)/del(DIR_Y)

     enddo
     enddo
     enddo 


     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, JAXIS, &
                           flxint_v, dataSize)             ! dx*dz     

#if NDIM == 3


     ! Z direction:
     do k=sz,sz
     do j=sy,ey
     do i=sx,ex

     !- kpd - This is the inverse of the mixture density on the cell face...
     Mdens = ( facezData(RH1F_FACE_VAR,i,j,k) + facezData(RH2F_FACE_VAR,i,j,k) ) ! Inverse density in face.

     flxint_w(VELC_FLUX,i,j,k) =    &
     Mdens*(solnData(pvar,i,j,k)- &
            solnData(pvar,i,j,k-1))*del(DIR_X)*del(DIR_Y)/del(DIR_Z)

     enddo
     enddo
     enddo

     do k=ez+1,ez+1
     do j=sy,ey
     do i=sx,ex

     !- kpd - This is the inverse of the mixture density on the cell face...
     Mdens = ( facezData(RH1F_FACE_VAR,i,j,k) + facezData(RH2F_FACE_VAR,i,j,k) ) ! Inverse density in face.

     flxint_w(VELC_FLUX,i,j,k) =    &
     Mdens*(solnData(pvar,i,j,k)- &
            solnData(pvar,i,j,k-1))*del(DIR_X)*del(DIR_Y)/del(DIR_Z)

     enddo      
     enddo      
     enddo      

     !- kpd - Copy the fluxes into the paraMesh data structure
     call Grid_putFluxData(blockID, KAXIS, &
                           flxint_w, dataSize)             ! dx*dy 
#endif

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

     !- kpd - 
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

  enddo


  ! fix fluxes at refinement interfaces
  !- kpd - At fine-coarse interface, set flux = to sum of fine block fluxes
  call Grid_conserveFluxes( ALLDIR, level)

  return

END SUBROUTINE ins_fluxfix_p


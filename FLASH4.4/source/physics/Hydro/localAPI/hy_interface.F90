!!****ih* source/physics/Hydro/HydroMain/unsplit/hy_interface
!!
!! NAME
!!  hy_interface
!!
!! SYNOPSIS
!!  use hy_interface
!!
!! DESCRIPTION
!!  This is an interface specific for the two unsplit solvers:
!!   1. MHD_StaggeredMesh
!!   2. Hydro_Unsplit
!!  This interface defines its public interfaces.
!!
!!***

!!REORDER(4): [xyz]flux, faceFlux, fl[xyz], Sp, U

Module hy_interface

  implicit none


#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  interface
     subroutine hy_avgState(sweepDir,VL,VR,Vavg)
       implicit none
       integer, intent(IN) :: sweepDir
       real, dimension(HY_VARINUM3), intent(IN)  :: VL,VR
       real, dimension(HY_VARINUM2), intent(OUT) :: Vavg
     end subroutine hy_avgState
  end interface



  interface
     subroutine hy_getRiemannState(block,U,blkLimits,loGC,hiGC,dt,del, &
                                       ogravX,ogravY,ogravZ,&
                                       scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,&
                                       hy_SpcR,&
                                       hy_SpcL,hy_SpcSig,normalFieldUpdate)

       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t), intent(IN)   :: block
       integer, intent(IN),dimension(LOW:HIGH,MDIM):: blkLimits !, blkLimitsGC
       integer, intent(IN), dimension(MDIM):: loGC, hiGC
       real,    intent(IN)   :: dt
       real,    intent(IN),dimension(MDIM) :: del
!!$       real, dimension(blkLimitsGC(HIGH,IAXIS),  &
!!$                       blkLimitsGC(HIGH,JAXIS),  &
!!$                       blkLimitsGC(HIGH,KAXIS)), &
!!$                       intent(IN) :: ogravX,ogravY,ogravZ
  real, dimension(loGC(IAXIS):,loGC(JAXIS):,loGC(KAXIS):), intent(IN) :: ogravX,ogravY,ogravZ
       real,pointer,dimension(:,:,:,:)::U
       real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr
       real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
       logical, intent(IN), optional :: normalFieldUpdate
     end subroutine hy_getRiemannState
 end interface



  interface
     subroutine hy_entropyFix(lambda,lambdaL,lambdaR)
       implicit none
       real, dimension(HY_WAVENUM), intent(INOUT) :: lambda
       real, dimension(HY_WAVENUM), intent(IN)    :: lambdaL,lambdaR
     end subroutine hy_entropyFix
  end interface

  interface
     subroutine hy_advance(simTime, dt, dtOld)
       implicit none
       real, intent(IN) ::  simTime, dt, dtOld
     end subroutine hy_advance
  end interface

  interface
     subroutine hy_getFaceFlux ( block,blkLimits,blkLimitsGC,datasize,del,&
                                     loFl, xflux,yflux,zflux,&
                                     scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,scrch_Ptr,&
                                     hy_SpcR,hy_SpcL,hy_SpcSig,lastCall)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t), intent(IN)  :: block
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
       integer, dimension(MDIM), intent(IN)         :: datasize
       real,    dimension(MDIM), intent(IN)         :: del
       integer, dimension(MDIM+1),intent(IN)        :: loFl
       real, intent(OUT) :: XFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZED TO PREVENT INDEX REORDERING!
       real, intent(OUT) :: YFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZATION INTENTIONAL!
       real, intent(OUT) :: ZFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZATION INTENTIONAL!
       real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr
       real, pointer, dimension(:,:,:,:) :: scrch_Ptr
       real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
       logical, optional, intent(IN) :: lastCall
     end subroutine hy_getFaceFlux
  end interface



  interface
     Subroutine hy_dataReconstOneStep(block,U,loGC,hiGC,order,ix,iy,iz, &
                                          dt,del,ogravX,ogravY,ogravZ,&
                                          DivU,FlatCoeff,   &
                                          TransX_updateOnly,&
                                          TransY_updateOnly,&
                                          TransZ_updateOnly,&
                                          Wp,Wn,sig,&
                                          lambda,leftEig,rghtEig,&
                                          cellCfl,&
                                          hy_SpcR,hy_SpcL,hy_SpcSig)
       use block_metadata,   ONLY : block_metadata_t

       implicit none
       type(block_metadata_t),intent(IN) :: block
       real,dimension(:,:,:,:),pointer :: U
       integer, intent(IN), dimension(MDIM):: loGC, hiGC
       integer,intent(IN) :: order,ix,iy,iz
       real,   intent(IN) :: dt
       real,   intent(IN), dimension(MDIM) :: del
       real, dimension(loGC(IAXIS):hiGC(IAXIS),loGC(JAXIS):hiGC(JAXIS),loGC(KAXIS):hiGC(KAXIS)), &
            intent(IN), target :: ogravX,ogravY,ogravZ
       real, dimension(NDIM,loGC(IAXIS):hiGC(IAXIS),loGC(JAXIS):hiGC(JAXIS),loGC(KAXIS):hiGC(KAXIS)), &
            intent(IN) :: FlatCoeff
       real, dimension(loGC(IAXIS):hiGC(IAXIS),loGC(JAXIS):hiGC(JAXIS),loGC(KAXIS):hiGC(KAXIS)), &
            intent(IN) :: DivU
       logical, intent(IN) ::  TransX_updateOnly, TransY_updateOnly, TransZ_updateOnly
       real, dimension(HY_VARINUMMAX,           NDIM),intent(OUT) :: Wp, Wn
       real, dimension(HY_VARINUMMAX,           NDIM),intent(OUT) :: sig
       real, dimension(              HY_WAVENUM,NDIM),intent(OUT) :: lambda
       real, dimension(HY_VARINUM,   HY_WAVENUM,NDIM),intent(OUT) :: leftEig
       real, dimension(HY_VARINUM,   HY_WAVENUM,NDIM),intent(OUT) :: rghtEig
       real, optional, intent(INOUT) :: cellCfl ! We may pass back a lowered CFL factor here - KW
       real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
     end subroutine hy_dataReconstOnestep
  end interface



  interface
     Subroutine hy_DataReconstructNormalDir_MH&
          (dir,dt,delta,Data1D,DataGrav1D,&
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: dir
       real,   intent(IN) :: dt,delta
       real, pointer, dimension(:,:) :: Data1D
       real, pointer, dimension(:)   :: DataGrav1D
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     End Subroutine Hy_DataReconstructNormalDir_MH
  end interface




  interface
     Subroutine hy_DataReconstructNormalDir_PPM&
          (dir,dt,delta,Data1D,DataGrav1D,&
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: dir
       real,   intent(IN) :: dt,delta
       real, pointer, dimension(:,:) :: Data1D
       real, pointer, dimension(:)   :: DataGrav1D
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     End Subroutine Hy_DataReconstructNormalDir_PPM
  end interface



  interface
     Subroutine hy_DataReconstructNormalDir_WENO&
          (wenoMethod,dir,dt,delta,Data1D,DataGrav1D,&
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: wenoMethod,dir
       real,   intent(IN) :: dt,delta
       real, pointer, dimension(:,:) :: Data1D
       real, pointer, dimension(:)   :: DataGrav1D
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     End Subroutine Hy_DataReconstructNormalDir_WENO
  end interface



  interface
     Subroutine hy_DataReconstructNormalDir_GP&
          (dir,dt,delta,DataMultiD,DataGravMultiD,&
          x1,x2,x3,radius, &
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: dir
       real,   intent(IN) :: dt,delta
#if NDIM < 3
       ! Although we define this array for 1D here for the sake of compilations,
       ! GP does not support 1D calculations.
       real, pointer, dimension(:,:,:)   :: DataMultiD,DataGravMultiD
#elif NDIM == 3
       real, pointer, dimension(:,:,:,:) :: DataMultiD,DataGravMultiD
#endif
       real, pointer, dimension(:) :: x1
       real, pointer, dimension(:) :: x2
       real, pointer, dimension(:) :: x3
       integer, intent(IN) :: radius
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     end Subroutine hy_DataReconstructNormalDir_GP
    end interface



  interface
     subroutine hy_Roe(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm,Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_Roe
  end interface



  interface
     subroutine hy_HLL(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1),  intent(OUT) :: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_HLL
  end interface



  interface
     subroutine hy_HLLC(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_HLLC
  end interface



  interface
     subroutine hy_Marquina(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_Marquina
  end interface


  interface
     subroutine hy_MarquinaModified(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_MarquinaModified
  end interface


  interface
     subroutine hy_LLF(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_LLF
  end interface



  interface
     subroutine hy_setMinTimeStep(block,i,j,k,delta,speed)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t), intent(IN)   :: block
       integer, INTENT(in) :: i,j,k
       real, INTENT(in) :: delta,speed
     end subroutine hy_setMinTimeStep
  end interface




!!$  interface
!!$     subroutine hy_upwindTransverseFlux(order,vector,lambda,leig,reig,sig)
!!$       implicit none
!!$       integer,intent(IN) :: order
!!$       real,intent(IN),dimension(HY_VARINUMMAX,-2:2)   :: vector
!!$       real,intent(IN),dimension(HY_WAVENUM) :: lambda
!!$       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: leig
!!$       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: reig
!!$       real,intent(OUT), dimension(HY_VARINUMMAX) :: sig
!!$     end subroutine hy_upwindTransverseFlux
!!$  end interface


  interface
     subroutine hy_upwindTransverseFlux&
          (dir,order,vm1,vc0,vp1,lambda,leig,reig,sigSize,sig,speciesScalar)
       implicit none
       integer,intent(IN) :: dir,order
       real,pointer,dimension(:)  :: vm1,vc0,vp1
       real,intent(IN),dimension(HY_WAVENUM) :: lambda
       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: leig
       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: reig
!!$       real,pointer,dimension(:) :: lambda
!!$       real,pointer,dimension(:,:) :: leig,reig
       integer,intent(IN) :: sigSize
       real,intent(OUT),dimension(sigSize) :: sig
       logical, intent(IN),optional :: speciesScalar
     end subroutine hy_upwindTransverseFlux
  end interface



  interface
     subroutine hy_TVDslope(dir,VL,V0,VR,lambda,leig,delbar)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX),intent(IN)  :: VL,V0,VR
       real, dimension(HY_WAVENUM),   intent(IN)  :: lambda
       real, dimension(HY_VARINUM,HY_WAVENUM), intent(IN) :: leig
       real, dimension(HY_VARINUMMAX),intent(OUT) :: delbar
     end subroutine hy_TVDslope
  end interface


  interface
     subroutine hy_TVDslopeUpwind(dir,VLL,VL,V0,VR,VRR,lambdaL,lambda,lambdaR,leig,delbar)
       implicit none
       integer, intent(IN) :: dir
       real,dimension(HY_VARINUMMAX),intent(IN)  :: VLL,VL,V0,VR,VRR
       real,dimension(HY_WAVENUM),  intent(IN)  :: lambdaL,lambda,lambdaR
       real,dimension(HY_VARINUM,HY_WAVENUM),intent(IN) :: leig
       real,dimension(HY_VARINUMMAX),intent(OUT) :: delbar
     end subroutine hy_TVDslopeUpwind
  end interface



  interface
     subroutine hy_unsplit( blockDesc,Uin,blkLimitsGC,&
                      Uout,blkLimits,&
                      del,dt, dtOld )
       use block_metadata,   ONLY : block_metadata_t
       implicit none

       type(block_metadata_t), intent(IN) :: blockDesc
       real,    INTENT(IN) :: dt, dtOld
       integer,dimension(LOW:HIGH,MDIM),INTENT(IN) :: blkLimits,blkLimitsGC
       real, dimension(:,:,:,:),pointer :: Uin
       real,dimension(:,:,:,:), pointer :: Uout
       real,dimension(MDIM),INTENT(IN) :: del
     end subroutine hy_unsplit
  end interface



  interface
     subroutine hy_unsplitUpdate(blockDesc,Uin,Uout,rangeSwitch,dt,del,dataSize,blkLimits,&    
                                     blGC,loFl,xflux,yflux,zflux,gravX,gravY,gravZ,&
                                     scrch_Ptr)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t), intent(IN)   :: blockDesc
       real,dimension(:,:,:,:),pointer :: Uin, Uout
       integer,intent(IN) :: rangeSwitch
       real, intent(IN)   :: dt
       real, intent(IN)   :: del(MDIM)
       integer,dimension(MDIM),intent(IN) :: dataSize
       integer,intent(IN) :: blkLimits(LOW:HIGH,MDIM)
       integer,intent(IN) :: blGC(LOW:HIGH,MDIM)
       integer, dimension(MDIM+1),intent(IN)        :: loFl
       real, intent(IN) :: XFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZED TO PREVENT INDEX REORDERING!
       real, intent(IN) :: YFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZATION INTENTIONAL!
       real, intent(IN) :: ZFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZATION INTENTIONAL!
    real, dimension(blGC(LOW,IAXIS):blGC(HIGH,IAXIS),blGC(LOW,JAXIS):blGC(HIGH,JAXIS),blGC(LOW,KAXIS):blGC(HIGH,KAXIS)),& 
         intent(IN) :: gravX,gravY,gravZ
       real, pointer, dimension(:,:,:,:) :: scrch_Ptr
     end subroutine hy_unsplitUpdate
  end interface





    interface
       subroutine hy_unsplitUpdateMultiTemp&
            (block,rangeSwitch,blkLimits,datasize,dt,del,xflux,yflux,zflux, scrch_Ptr)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t), intent(IN)   :: block
       integer,intent(IN) :: rangeSwitch
       integer,intent(IN) :: blkLimits(LOW:HIGH,MDIM)
       integer,intent(IN) :: datasize(MDIM)      
       real, intent(IN) :: dt
       real, intent(IN) :: del(MDIM)
#ifdef FIXEDBLOCKSIZE
         real, intent(in) :: xflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
         real, intent(in) :: yflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
         real, intent(in) :: zflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
#else
         real, intent(in) :: xflux(NFLUXES,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))  
         real, intent(in) :: yflux(NFLUXES,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))  
         real, intent(in) :: zflux(NFLUXES,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))
#endif
         real, pointer, dimension(:,:,:,:) :: scrch_Ptr
       end subroutine hy_unsplitUpdateMultiTemp
    end interface

    interface
       subroutine hy_multiTempAfter(blockCount, blockList, dt)
         implicit none
         integer, intent(in) :: blockCount
         integer, intent(in) :: blockList(blockCount)
         real,    intent(in) :: dt
       end subroutine hy_multiTempAfter
    end interface

    interface
       subroutine hy_ragelike(uextra, soln, dens_old, &
            duele_adv, duion_adv, durad_adv, &
            xc, yc, zc, &
            uele_new, uion_new, urad_new, &
            stepFactor)
         implicit none
         real, intent(in)  :: uextra
         real, intent(in)  :: soln(NUNK_VARS)
         real, intent(in)  :: dens_old
         real, intent(in)  :: duele_adv
         real, intent(in)  :: duion_adv
         real, intent(in)  :: durad_adv
         real, intent(in)  :: xc 
         real, intent(in)  :: yc 
         real, intent(in)  :: zc
         real, intent(out) :: uele_new
         real, intent(out) :: uion_new
         real, intent(out) :: urad_new
         real, intent(in),OPTIONAL  :: stepFactor
       end subroutine hy_ragelike
    end interface

  interface
     subroutine hy_updateSpeciesMassScalar&
          (order,densNew,Sp,U,FL,FR,GL,GR,HL,HR,dx,dy,dz,dt,SpNew)
  implicit none
  integer, intent(IN) :: order
  real, intent(IN) :: densNew
#if NDIM == 1
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,1,1)   :: Sp
  real, intent(IN), dimension(6,-3:3,1,1) :: U
#elif NDIM == 2
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,-3:3,1)   :: Sp
  real, intent(IN), dimension(6,-3:3,-3:3,1) :: U
#elif NDIM == 3
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,-3:3,-3:3)   :: Sp
  real, intent(IN), dimension(6,-3:3,-3:3,-3:3)   :: U
#endif
  real, intent(IN)  :: FL,FR,GL,GR,HL,HR,dx,dy,dz,dt
  real, intent(OUT), dimension(NSPECIES+NMASS_SCALARS) :: SpNew
     end subroutine hy_updateSpeciesMassScalar
  end interface


  interface
     subroutine hy_energyFix(block,U,blkLimits,dt,dtOld,del,eosMode)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       
       type(block_metadata_t), intent(IN) :: block
       real, pointer, dimension(:,:,:,:) :: U
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
       real, intent(IN) :: dt,dtOld
       real, dimension(MDIM), intent(IN) :: del
       integer, intent(IN) :: eosMode
     end subroutine hy_energyFix
  end interface



  interface 
     subroutine hy_unitConvert(Uin,blkLimitsGC,convertDir)
       implicit none
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
       integer, intent(IN) :: convertDir
       real, dimension(:,:,:,:) :: Uin
     end subroutine hy_unitConvert
  end interface
  


  interface
     subroutine hy_addViscousFluxes&
          (block,blkLimitsGC,ix,iy,iz,Flux,mu,sweepDir)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t), intent(IN)   :: block
       integer, INTENT(IN) :: ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE 
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, & 
                  GRID_JLO_GC:GRID_JHI_GC, & 
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: mu
#else 
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), & 
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), & 
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),& 
                  intent(IN) :: mu
#endif 
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_addViscousFluxes
    end interface



  interface
     subroutine hy_addThermalFluxes&
          (block,blkLimitsGC,ix,iy,iz,Flux,kappa,sweepDir)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t), intent(IN)   :: block
       integer, INTENT(IN) :: ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE 
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, & 
                  GRID_JLO_GC:GRID_JHI_GC, & 
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: kappa
#else 
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), & 
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), & 
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),& 
                  intent(IN) :: kappa
#endif 
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_addThermalFluxes
    end interface



    interface
       subroutine hy_eigenParameters&
            (V,dir,U_normal,C_fast,C_alfn,C_slow,A_f,A_s,B_beta,C_hyp)
         implicit none
         real, dimension(HY_VARINUM2), intent(IN)  :: V
         integer, intent(IN)  :: dir
         real, intent(OUT) :: U_normal,C_fast
         real, intent(OUT), optional :: C_alfn,C_slow,A_f,A_s
         real, dimension(MDIM), intent(OUT),optional :: B_beta
         real, intent(OUT), optional :: C_hyp
       end subroutine hy_eigenParameters
    end interface



    interface
       subroutine hy_eigenValue(EigValue,U_normal,C_fast,C_alfn,C_slow,C_hyp)
         implicit none
         real,dimension(HY_WAVENUM), intent(OUT) :: EigValue
         real,intent(IN) :: U_normal,C_fast
         real,intent(IN), optional :: C_alfn,C_slow,C_hyp
       end subroutine hy_eigenValue
    end interface



    interface
       subroutine hy_eigenVector&
            (LeftEigvec,RightEigvec,V,dir,cons,C_fast,C_alfn,C_slow,A_f,A_s,B_beta,C_hyp)
         implicit none
         real, dimension(HY_VARINUM,HY_WAVENUM), intent(OUT) :: LeftEigvec
         real, dimension(HY_VARINUM,HY_WAVENUM), intent(OUT) :: RightEigvec
         real, dimension(HY_VARINUM2), intent(IN) :: V
         integer, intent(IN) :: dir
         logical, intent(IN) :: cons
         real,    intent(IN) :: C_fast
         real,    intent(IN), optional :: C_alfn,C_slow,A_f,A_s
         real, dimension(MDIM), intent(IN), optional  :: B_beta
         real,    intent(IN), optional :: C_hyp
       end subroutine hy_eigenVector
    end interface


    interface
       subroutine hy_prepareNewGravityAccel(gcMaskLogged)
         logical,intent(in) :: gcMaskLogged
       end subroutine hy_prepareNewGravityAccel
    end interface


    interface
       subroutine hy_putGravity&
            (blockDesc,blGC,Uin,dataSize,dt,dtOld,gravX,gravY,gravZ,potentialIndex,&
             lastCall)
         use block_metadata, ONLY : block_metadata_t
         implicit none
         type(block_metadata_t), intent(IN)   :: blockDesc
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blGC
         real,dimension(:,:,:,:),pointer :: Uin
         integer, dimension(MDIM), intent(IN) :: dataSize
         real,    intent(IN) :: dt, dtOld
         real,dimension(blGC(LOW,IAXIS):blGC(HIGH,IAXIS), blGC(LOW,JAXIS):blGC(HIGH,JAXIS), blGC(LOW,KAXIS):blGC(HIGH,KAXIS)),&
              intent(OUT) :: gravX,gravY,gravZ
         integer, intent(IN), OPTIONAL :: potentialIndex
         logical, intent(IN), OPTIONAL :: lastCall
       end subroutine hy_putGravity
    end interface



    interface
       Subroutine hy_addGravity&
            (blockDesc,blkLimits,loGC,hiGC,dt,gravX,gravY,gravZ)
         use block_metadata, ONLY : block_metadata_t
         implicit none
         type(block_metadata_t), intent(IN)   :: blockDesc
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
         integer, intent(IN), dimension(MDIM):: loGC, hiGC
         real,    intent(IN) :: dt
         real, dimension(loGC(IAXIS):hiGC(IAXIS),loGC(JAXIS):hiGC(JAXIS),loGC(KAXIS):hiGC(KAXIS)), intent(IN) :: &
              gravX,gravY,gravZ
       end Subroutine hy_addGravity
    end interface



    interface
       subroutine hy_shockDetectBlk(Uin,blkLimitsGC,Uout,blkLimits,del )
       implicit none
       integer,dimension(LOW:HIGH,MDIM),INTENT(IN) :: blkLimits,blkLimitsGC
       real, dimension(:,:,:,:),INTENT(IN) :: Uin
       real,dimension(:,:,:,:),INTENT(INOUT) :: Uout
       real,dimension(MDIM),INTENT(IN) :: del
     end subroutine hy_shockDetectBlk
  end interface

  interface
     subroutine hy_shockDetect()
       implicit none
     end subroutine hy_shockDetect
  end interface

    interface
       subroutine hy_prim2con(V,CU)
         implicit none
         real ,dimension(HY_VARINUM2), intent(IN) :: V
         real ,dimension(HY_VARINUM),  intent(OUT) :: CU
       end subroutine hy_prim2con
    end interface



    interface
       subroutine hy_con2prim(CU,game,V)
         implicit none
         real ,dimension(HY_VARINUM), intent(IN)  :: CU
         real, intent(IN) :: game
         real ,dimension(HY_VARINUM), intent(OUT) :: V
       end subroutine hy_con2prim
    end interface



    interface
       subroutine hy_prim2flx(dir,V,F)
         implicit none
         integer, intent(IN)  :: dir
         real, dimension(*),            intent(IN) :: V
         real, dimension(HY_VARINUM1),  intent(OUT) :: F
       end subroutine hy_prim2flx
    end interface



    interface
       Subroutine hy_checkRHjumpCond(dir,dens,velocity,pres,gamc,Wp,Wn,SWp,SWn)
         implicit none
         integer, intent(IN)  :: dir
         real, dimension(MDIM), intent(IN) :: velocity
         real, intent(IN) :: dens, pres, gamc
         real, dimension(HY_DENS:HY_EINT), intent(IN) :: Wp,Wn
         logical, intent(OUT) :: SWp,SWn
       end Subroutine hy_checkRHjumpCond
    end interface

    interface
     subroutine hy_gravityStep(simTime, dt, dtOld)
       implicit none
       real, intent(IN) ::  simTime, dt, dtOld
     end subroutine Hy_gravityStep
     
     subroutine hy_gravityStepBlk(blockDesc, blkLimitsGC, Uin, blkLimits, Uout, del,timeEndAdv,dt,dtOld)
       use block_metadata, ONLY : block_metadata_t
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       
       real, pointer, dimension(:,:,:,:) :: Uout
       real, pointer, dimension(:,:,:,:) :: Uin
       
       real,dimension(MDIM),intent(IN) :: del
       integer,dimension(LOW:HIGH,MDIM),intent(INoUt) ::blkLimits,blkLimitsGC 
       type(block_metadata_t), intent(IN) :: blockDesc
     end subroutine hy_gravityStepBlk
  end interface

#ifdef FLASH_UGLM_MHD
    interface
       Subroutine hy_updateSourceGLM(block,dt,del,blkLimits,blkLimitsGC,C_hyp,C_par)
       use block_metadata, ONLY : block_metadata_t
         implicit none
         type(block_metadata_t), intent(IN)   :: block
         real,    intent(IN) :: dt
         real,    dimension(MDIM),   intent(IN) :: del
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC
         real,    intent(IN) :: C_hyp,C_par
       end Subroutine hy_updateSourceGLM
    end interface
#endif


  ! Dongwook thinks that the following two interfaces should be taken care of for
  ! general implementation later.
!!$  interface
!!$     Subroutine hy_counterGP_init(radiusGP, counterGP, blkLimitsGC)
!!$       implicit none
!!$       real,    intent(IN)  :: radiusGP
!!$       integer, intent(OUT) :: counterGP
!!$       integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimitsGC
!!$     end subroutine hy_counterGP_init
!!$  end interface
!!$
  interface
     Subroutine hy_initGP(RinvGP, WpGP, WmGP, blkLimitsGC)
       implicit none
       real, intent(INOUT), dimension(:,:) :: RinvGP
       real, intent(INOUT), dimension(:,:) :: WpGP
       real, intent(INOUT), dimension(:,:) :: WmGP
       integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimitsGC
     end Subroutine hy_initGP
  end interface


!! FOR UNSPLIT STAGGERED MESH MHD SOLVER -------------------------------------------
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  interface
     subroutine hy_addResistiveFluxes&
          (block,blkLimitsGC,ix,iy,iz,Flux,eta,sweepDir)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t), intent(IN)   :: block
       integer, INTENT(IN) :: block,ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE 
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, & 
                  GRID_JLO_GC:GRID_JHI_GC, & 
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: eta
#else 
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), & 
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), & 
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),& 
                  intent(IN) :: eta
#endif 
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_addResistiveFluxes
    end interface



  interface
     subroutine hy_addOhmicHeating&
          (block,blkLimits,ix,iy,iz,Qohm,eta)
       use block_metadata, ONLY : block_metadata_t
     implicit none
     type(block_metadata_t), intent(IN)   :: block
     integer, INTENT(IN) :: block,ix,iy,iz
     integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits
     real, intent(INOUT) :: Qohm
     real, intent(IN)    :: eta  
     end subroutine hy_addOhmicHeating
  end interface


  interface
     subroutine hy_staggeredDivb(block,dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t), intent(IN)   :: block
       integer, intent(IN) :: block
       real,    intent(IN) :: dt
       real,    dimension(MDIM),   intent(IN) :: del
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC
       logical, intent(IN) :: halfTimeAdvance
     end subroutine hy_staggeredDivb
  end interface



  interface
     subroutine hy_HLLD(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_HLLD
  end interface



  interface
     subroutine hy_getElectricFields( block,blkLimits,blkLimitsGC,del,flx,fly,flz)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t), intent(IN)   :: block
       integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
       real,    intent(IN), dimension(MDIM)  :: del  
#ifdef FIXEDBLOCKSIZE
       real, DIMENSION(NFLUXES,               &
                       GRID_ILO_GC:GRID_IHI_GC,  &
                       GRID_JLO_GC:GRID_JHI_GC,  &
                       GRID_KLO_GC:GRID_KHI_GC), &
                       intent(IN) :: flx,fly,flz
#else
       real, DIMENSION(NFLUXES,               &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), &
                       intent(IN) :: flx,fly,flz
#endif
     end subroutine hy_getElectricFields
  end interface



  interface
     subroutine hy_getCurrents( &
        block, rangeSwitch, blkLimits,datasize, del, Jp, Jm, mode_switch,&
        scrch_Ptr,&
        ix,iy,iz)
       use block_metadata, ONLY : block_metadata_t
        implicit none
        ! Arguments:
         type(block_metadata_t), intent(IN)   :: block
        integer,intent(in) :: block, rangeSwitch
        integer,intent(in) :: blkLimits  (LOW:HIGH,MDIM)
        integer,intent(in) :: datasize(MDIM), mode_switch
        real,   intent(in) :: del(MDIM)
#ifdef FIXEDBLOCKSIZE
        real, intent(inout) :: Jp(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
        real, intent(inout) :: Jm(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
#else
        real, intent(inout) :: Jp(3,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))  
        real, intent(inout) :: Jm(3,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))
#endif
        real, pointer, dimension(:,:,:,:) :: scrch_Ptr
        integer, intent(in), OPTIONAL :: ix,iy,iz 
     end subroutine hy_getCurrents
  end interface



  interface
     subroutine hy_getFluxDeriv( ix,iy,iz,blkLimitsGC,&
                                     fluxType,DerivDir,   &
                                     faceFlux,            &
                                     Flux1Deriv,          &
                                     Flux2Deriv           )
       implicit none
       integer, intent(IN) :: ix,iy,iz
       integer, intent(IN) :: fluxType,DerivDir
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real,    dimension(NFLUXES,                  &
                          GRID_ILO_GC:GRID_IHI_GC,  &
                          GRID_JLO_GC:GRID_JHI_GC,  &
                          GRID_KLO_GC:GRID_KHI_GC), &
                          intent(IN) :: faceFlux
#else
       real,    dimension(NFLUXES,                  &
                          blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), &
                          intent(IN) :: faceFlux
#endif
       real,    intent(OUT):: Flux1Deriv,Flux2Deriv
     end subroutine hy_getFluxDeriv
  end interface



  interface
     Subroutine hy_addBiermannBatteryTerms(block,blkLimitsGC,ix,iy,iz,Flux,sweepDir)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       integer, INTENT(IN) :: ix,iy,iz
       type(block_metadata_t), intent(IN)   :: block
       
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
       integer, INTENT(IN) :: sweepDir
     end Subroutine hy_addBiermannBatteryTerms
  end interface
  
  interface
     subroutine hy_biermannSource(blockCount, blockList, dt)
       implicit none
       integer, INTENT(IN) :: blockCount  
       integer, INTENT(IN), dimension(blockCount) :: blockList
       real,    INTENT(IN) :: dt
     end subroutine hy_biermannSource
  end interface
  interface
     subroutine hy_putGravity(blkLimitsGC,Uin,dataSize,dt,dtOld,gravX,gravY,gravZ)
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
       integer, dimension(MDIM), intent(IN) :: dataSize
       real,    intent(IN) :: dt, dtOld
       real,dimension(:,:,:,:),pointer :: Uin
     end subroutine hy_putGravity
  end interface
#endif 
!endif #ifdef FLASH_USM_MHD

  interface
     subroutine hy_computeFluxes(blockDesc, blkLimitsGC, Uin, blkLimits, Uout, del,timeEndAdv, dt, dtOld,  &
          sweepOrder )
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t) :: blockDesc
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
       real, pointer, dimension(:,:,:,:) :: Uout,Uin
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       integer, INTENT(IN) :: sweepOrder
       real,dimension(MDIM),intent(IN) :: del
     end subroutine hy_computeFluxes

     subroutine hy_advanceBlk(blockDesc, blkLimitsGC, Uin, blkLimits, Uout, del,timeEndAdv, dt, dtOld,  &
          sweepOrder )
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t) :: blockDesc
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
       real, pointer, dimension(:,:,:,:) :: Uout,Uin
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       integer, INTENT(IN) :: sweepOrder
       real,dimension(MDIM),intent(IN) :: del
     end subroutine hy_advanceBlk
  end interface

    interface
     subroutine hy_updateSolution(blockDesc, blkLimitsGC, Uin, blkLimits, Uout, del,timeEndAdv, dt, dtOld,  &
          sweepOrder )
       use block_metadata, ONLY : block_metadata_t
       implicit none
       type(block_metadata_t) :: blockDesc
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
       real, pointer, dimension(:,:,:,:) :: Uout,Uin
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       integer, INTENT(IN) :: sweepOrder
       real,dimension(MDIM),intent(IN) :: del
     end subroutine hy_updateSolution
  end interface

End Module hy_interface

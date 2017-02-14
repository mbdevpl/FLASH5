!!****ih* source/physics/Hydro/HydroMain/localAPI/hy_simpleInterface
!!
!! NAME
!!  hy_simpleInterface
!!
!! SYNOPSIS
!!  use hy_simpleInterface
!!
!! DESCRIPTION
!!  This is a local interface specific for simple unsplit solvers:
!!
!!***

!!REORDER(4): [xyz]flux, faceFlux, fl[xyz], Sp, U

Module hy_simpleInterface

  implicit none


#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"

#ifdef FLASH_HYDRO_UNSPLIT
#include "UHD.h"

!#ifdef FLASH_USM_MHD
!#define HY_VARINUM 8
!#define HY_WAVENUM 7
!#else
!#define HY_VARINUM 5
!#define HY_WAVENUM 5
!#endif
!#define HY_VARINUM2 HY_VARINUM+2
!#define HY_VARINUM3 HY_VARINUM+3
!#define HY_VARINUMMAX HY_VARINUM+4


  interface
     subroutine hy_uhd_avgState(sweepDir,VL,VR,Vavg)
       implicit none
       integer, intent(IN) :: sweepDir
       real, dimension(HY_VARINUM3), intent(IN)  :: VL,VR
       real, dimension(HY_VARINUM2), intent(OUT) :: Vavg
     end subroutine hy_uhd_avgState
  end interface



  interface
     subroutine hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del, &
                                       ogravX,ogravY,ogravZ,&
                                       hgravX,hgravY,hgravZ,&
                                       normalFieldUpdate)

       implicit none
       integer, intent(IN)   :: blockID
       integer, intent(IN),dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
       real,    intent(IN)   :: dt
       real,    intent(IN),dimension(MDIM) :: del
#ifdef FIXEDBLOCKSIZE
       real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) &
            :: ogravX,ogravY,ogravZ,hgravX,hgravY,hgravZ
#else
       real, dimension(blkLimitsGC(HIGH,IAXIS),  &
                       blkLimitsGC(HIGH,JAXIS),  &
                       blkLimitsGC(HIGH,KAXIS)), &
                       intent(IN) :: ogravX,ogravY,ogravZ,hgravX,hgravY,hgravZ
#endif
       logical, intent(IN), optional :: normalFieldUpdate
     end subroutine hy_uhd_getRiemannState
 end interface



  interface
     subroutine hy_uhd_entropyFix(lambda,lambdaL,lambdaR)
       implicit none
       real, dimension(HY_WAVENUM), intent(INOUT) :: lambda
       real, dimension(HY_WAVENUM), intent(IN)    :: lambdaL,lambdaR
     end subroutine hy_uhd_entropyFix
  end interface



  interface
     subroutine hy_uhd_getFaceFlux ( blockID,blkLimits,blkLimitsGC, datasize, del,&
                                     xflux, yflux,zflux,lastCall)
       implicit none
       integer, intent(IN)  :: blockID
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
       integer, dimension(MDIM), intent(IN)         :: datasize
       real,    dimension(MDIM), intent(IN)         :: del
#ifdef FIXEDBLOCKSIZE
       real, intent(OUT) :: xflux (NFLUXES, &
                       GRID_ILO_GC:GRID_IHI_GC, &
                       GRID_JLO_GC:GRID_JHI_GC, &
                       GRID_KLO_GC:GRID_KHI_GC)
       real, intent(OUT) :: yflux (NFLUXES, &
                       GRID_ILO_GC:GRID_IHI_GC, &
                       GRID_JLO_GC:GRID_JHI_GC, &
                       GRID_KLO_GC:GRID_KHI_GC)
       real, intent(OUT) :: zflux (NFLUXES, &
                       GRID_ILO_GC:GRID_IHI_GC, &
                       GRID_JLO_GC:GRID_JHI_GC, &
                       GRID_KLO_GC:GRID_KHI_GC)                
#else
       real, intent(OUT) :: xflux(NFLUXES, &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
       real, intent(OUT) :: yflux(NFLUXES, &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
       real, intent(OUT) :: zflux(NFLUXES, &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
#endif
       logical, optional, intent(IN) :: lastCall
     end subroutine hy_uhd_getFaceFlux
  end interface



  interface
     Subroutine hy_uhd_dataReconstOneStep(blockID,blkLimitsGC,ix,iy,iz, &
                                          dt,del,ogravX,ogravY,ogravZ,DivU,soundSpeed,V0, &
                                          Vxp,  Vxn,  Vyp,  Vyn,  Vzp,  Vzn,  &
                                          Vxpp, Vxnn, Vypp, Vynn, Vzpp, Vznn, &
                                          Vxppp,Vxnnn,Vyppp,Vynnn,Vzppp,Vznnn,&
                                          FlatCoeff,&
                                          TransX_updateOnly,&
                                          TransY_updateOnly,&
                                          TransZ_updateOnly,&
                                          Wxp, Wxn, Wyp, Wyn, Wzp, Wzn, &
                                          sig,lambda,leig,reig )

       implicit none
       integer,intent(IN) :: blockID
       integer, intent(IN),dimension(LOW:HIGH,MDIM):: blkLimitsGC
       integer,intent(IN) :: ix,iy,iz
       real,   intent(IN) :: dt
       real,   intent(IN), dimension(MDIM) :: del
#ifdef FIXEDBLOCKSIZE
       real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: ogravX,ogravY,ogravZ
       real, dimension(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: FlatCoeff
       real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: DivU,soundSpeed
#else
       real, dimension(blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
            intent(IN) :: ogravX,ogravY,ogravZ
       real, dimension(NDIM,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
            intent(IN) :: FlatCoeff
       real, dimension(blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
            intent(IN) :: DivU,soundSpeed
#endif
       logical, intent(IN) ::  TransX_updateOnly, TransY_updateOnly, TransZ_updateOnly
       real,   intent(INOUT),  dimension(HY_VARINUMMAX) :: V0, Vxp,  Vxn,  Vyp,  Vyn,  Vzp,  Vzn, &
                                                           Vxpp, Vxnn, Vypp, Vynn, Vzpp, Vznn,&
                                                           Vxppp,Vxnnn,Vyppp,Vynnn,Vzppp,Vznnn
       real,   intent(OUT), dimension(HY_VARINUMMAX) ::     Wxp, Wxn, Wyp, Wyn, Wzp, Wzn
       real, intent(OUT), dimension(HY_VARINUMMAX,NDIM) :: sig
       real, intent(OUT), dimension(HY_WAVENUM,NDIM) :: lambda
       real, intent(OUT), dimension(HY_WAVENUM,HY_VARINUM,NDIM) :: leig
       real, intent(OUT), dimension(HY_VARINUM,HY_WAVENUM,NDIM) :: reig

     end subroutine hy_uhd_dataReconstOnestep
  end interface



  interface
     subroutine hy_uhd_Roe(dir,Vm,Vp,Fstar,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm,Vp
       real, dimension(HY_VARINUM), intent(OUT):: Fstar
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_Roe
  end interface



  interface
     subroutine hy_uhd_HLL(dir,Vm,Vp,Fstar,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM),  intent(OUT) :: Fstar
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_HLL
  end interface



  interface
     subroutine hy_uhd_HLLC(dir,Vm,Vp,Fstar,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM), intent(OUT):: Fstar
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_HLLC
  end interface



  interface
     subroutine hy_uhd_Marquina(dir,Vm,Vp,Fstar,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM), intent(OUT):: Fstar
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_Marquina
  end interface



  interface
     subroutine hy_uhd_LLF(dir,Vm,Vp,Fstar,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM), intent(OUT):: Fstar
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_LLF
  end interface



  interface
     subroutine hy_uhd_TVDslope(dir,VLL,VL,V0,VR,VRR,lambdaL,lambda,lambdaR,leig,delbar)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX),intent(IN)  :: VLL,VL,V0,VR,VRR
       real, dimension(HY_WAVENUM),  intent(IN)  :: lambdaL,lambda,lambdaR
       real, dimension(HY_WAVENUM,HY_VARINUM), intent(IN) :: leig
       real, dimension(HY_VARINUMMAX),intent(OUT) :: delbar
     end subroutine hy_uhd_TVDslope
  end interface


  interface
     subroutine hy_uhd_TVDslopeUpwind(dir,VLL,VL,V0,VR,VRR,lambdaL,lambda,lambdaR,leig,delbar)
       implicit none
       integer, intent(IN) :: dir
       real,dimension(HY_VARINUMMAX),intent(IN)  :: VLL,VL,V0,VR,VRR
       real,dimension(HY_WAVENUM),  intent(IN)  :: lambdaL,lambda,lambdaR
       real,dimension(HY_WAVENUM,HY_VARINUM),intent(IN) :: leig
       real,dimension(HY_VARINUMMAX),intent(OUT) :: delbar
     end subroutine hy_uhd_TVDslopeUpwind
  end interface



  interface
     subroutine hy_hllUnsplit( blockCount, blockList, dt, dtOld )
       implicit none
       integer, INTENT(IN) :: blockCount  
       integer, INTENT(IN), dimension(blockCount) :: blockList
       real,    INTENT(IN) :: dt, dtOld
     end subroutine hy_hllUnsplit
  end interface

  interface
     subroutine hy_llfUnsplit( blockCount, blockList, dt, dtOld )
       implicit none
       integer, INTENT(IN) :: blockCount  
       integer, INTENT(IN), dimension(blockCount) :: blockList
       real,    INTENT(IN) :: dt, dtOld
     end subroutine hy_llfUnsplit
  end interface



  interface
     subroutine hy_hllUnsplitUpdate(blockID,dt,dtOld,del,dataSize,blkLimits,&
                                     blkLimitsGC,xflux,yflux,zflux,gravX,gravY,gravZ,tilingCtx)
       use GridTilingModule, ONLY: Grid_tilingContext_t
       implicit none
       integer,intent(IN) :: blockID
       real, intent(IN)   :: dt,dtOld
       real, intent(IN)   :: del(MDIM)
       integer,dimension(MDIM),intent(IN) :: dataSize
       integer,intent(IN) :: blkLimits(LOW:HIGH,MDIM)
       integer,intent(IN) :: blkLimitsGC(LOW:HIGH,MDIM)       
#ifdef FIXEDBLOCKSIZE
  real, intent(in) ::xflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
  real, intent(in) ::yflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
  real, intent(in) ::zflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
  real, dimension(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: gravX,gravY,gravZ

#else
  real, intent(in) :: xflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS))  
  real, intent(in) :: yflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS))  
  real, intent(in) :: zflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS))
  real, dimension(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)),& 
                  intent(IN) :: gravX,gravY,gravZ
#endif
       type(Grid_tilingContext_t),intent(INOUT) :: tilingCtx
     end subroutine hy_hllUnsplitUpdate
  end interface





    interface
       subroutine hy_hllUnsplitUpdateMultiTemp&
            (blockID,range_switch,blkLimits,datasize,dt,del,xflux,yflux,zflux)
         implicit none
         integer,intent(IN) :: blockID, range_switch
         integer,intent(IN) :: blkLimits  (LOW:HIGH,MDIM)
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
       end subroutine hy_hllUnsplitUpdateMultiTemp
    end interface




  interface
     subroutine hy_uhd_updateSpeciesMassScalar&
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
     end subroutine hy_uhd_updateSpeciesMassScalar
  end interface


  interface
     subroutine hy_uhd_energyFix(blockID,blkLimits,dt,del,eosMode)
       implicit none
       integer, intent(IN) :: blockID
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
       real, intent(IN) :: dt
       real, dimension(MDIM), intent(IN) :: del
       integer, intent(IN) :: eosMode
     end subroutine hy_uhd_energyFix
  end interface



  interface 
     subroutine hy_hllUnitConvert(blockID,convertDir)
       implicit none
       integer, intent(IN) :: blockID, convertDir
     end subroutine hy_hllUnitConvert
  end interface
  


  interface
     subroutine hy_uhd_addViscousFluxes&
          (blockID,blkLimitsGC,ix,iy,iz,Flux,mu,sweepDir)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
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
     end subroutine hy_uhd_addViscousFluxes
    end interface



  interface
     subroutine hy_uhd_addThermalFluxes&
          (blockID,blkLimitsGC,ix,iy,iz,Flux,kappa,sweepDir)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
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
     end subroutine hy_uhd_addThermalFluxes
    end interface



    interface
       subroutine hy_uhd_eigenParameters(V,dir,cons,U_normal,C_fast,C_alfn,C_slow,A_f,A_s,B_beta)
         implicit none
         real, dimension(HY_VARINUM2), intent(IN)  :: V
         integer, intent(IN)  :: dir
         logical, intent(IN)  :: cons
         real, intent(OUT) :: U_normal,C_fast
         real, intent(OUT), optional :: C_alfn,C_slow,A_f,A_s
         real, dimension(MDIM), intent(OUT),optional :: B_beta
       end subroutine hy_uhd_eigenParameters
    end interface



    interface
       subroutine hy_uhd_eigenValue(EigValue,U_normal,C_fast,C_alfn,C_slow)
         implicit none
         real,dimension(HY_WAVENUM), intent(OUT) :: EigValue
         real,intent(IN) :: U_normal,C_fast
         real,intent(IN), optional :: C_alfn,C_slow
       end subroutine hy_uhd_eigenValue
    end interface



    interface
       subroutine hy_uhd_eigenVector&
            (LeftEigvec,RightEigvec,V,dir,cons,C_fast,C_alfn,C_slow,A_f,A_s,B_beta)
         implicit none
         real, dimension(HY_WAVENUM,HY_VARINUM), intent(OUT) :: LeftEigvec
         real, dimension(HY_VARINUM,HY_WAVENUM), intent(OUT) :: RightEigvec
         real, dimension(HY_VARINUM2), intent(IN) :: V
         integer, intent(IN) :: dir
         logical, intent(IN) :: cons
         real, intent(IN) :: C_fast
         real, intent(IN),optional :: C_alfn,C_slow,A_f,A_s
         real, dimension(MDIM), intent(IN),optional  :: B_beta
       end subroutine hy_uhd_eigenVector
    end interface



    interface
       subroutine hy_uhd_putGravityUnsplit&
            (blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
         implicit none
         integer, intent(IN) :: blockID
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
         integer, dimension(MDIM), intent(IN) :: dataSize
         real,    intent(IN) :: dt, dtOld
#ifdef FIXEDBLOCKSIZE
         real, dimension(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(INOUT) :: &
              gravX,gravY,gravZ
#else
         real, dimension(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(INOUT) :: &
              gravX,gravY,gravZ
#endif
       end subroutine hy_uhd_putGravityUnsplit
    end interface



    interface
       Subroutine hy_uhd_addGravityUnsplit&
            (blockID,blkLimitsGC,dataSize,dt,gravX,gravY,gravZ)
         implicit none
         integer, intent(IN) :: blockID
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
         integer, dimension(MDIM), intent(IN) :: dataSize
         real,    intent(IN) :: dt
#ifdef FIXEDBLOCKSIZE
         real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: &
              gravX,gravY,gravZ
#else
         real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(IN) :: &
              gravX,gravY,gravZ
#endif
       end Subroutine hy_uhd_addGravityUnsplit
    end interface



    interface
       subroutine hy_uhd_shockDetect(blockID)
         implicit none
         integer, INTENT(IN) :: blockID
       end subroutine hy_uhd_shockDetect
    end interface



    interface
       subroutine hy_uhd_prim2con(V,CU)
         implicit none
         real ,dimension(HY_VARINUM2), intent(IN) :: V
         real ,dimension(HY_VARINUM),  intent(OUT) :: CU
       end subroutine hy_uhd_prim2con
    end interface



    interface
       subroutine hy_uhd_con2prim(CU,game,V)
         implicit none
         real ,dimension(HY_VARINUM), intent(IN)  :: CU
         real, intent(IN) :: game
         real ,dimension(HY_VARINUM), intent(OUT) :: V
       end subroutine hy_uhd_con2prim
    end interface



    interface
       subroutine hy_uhd_prim2flx(dir,V,F)
         implicit none
         integer, intent(IN)  :: dir
         real, dimension(HY_VARINUM2), intent(IN) :: V
         real, dimension(HY_VARINUM),  intent(OUT) :: F
       end subroutine hy_uhd_prim2flx
    end interface




    interface
       subroutine hy_uhd_checkRHjumpCond(dir,idx,idy,idz,V0,Vxr,Vxl,Vyr,Vyl,Vzr,Vzl,Wr,Wl,SWr,SWl)
         implicit none
         integer, intent(IN) :: dir
         real, intent(IN) :: idx,idy,idz
         real, dimension(HY_VARINUMMAX), intent(IN)    :: V0,Vxr,Vxl,Vyr,Vyl,Vzr,Vzl
         real, dimension(HY_VARINUMMAX), intent(INOUT) :: Wr,Wl
         logical, intent(OUT) :: SWr,SWl
       end subroutine hy_uhd_checkRHjumpCond
    end interface



!! FOR UNSPLIT STAGGERED MESH MHD SOLVER -------------------------------------------
#ifdef FLASH_USM_MHD
  interface
     subroutine hy_uhd_addResistiveFluxes&
          (blockID,blkLimitsGC,ix,iy,iz,Flux,eta,sweepDir)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
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
     end subroutine hy_uhd_addResistiveFluxes
    end interface



  interface
     subroutine hy_uhd_staggeredDivb(blockID,dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)
       implicit none
       integer, intent(IN) :: blockID
       real,    intent(IN) :: dt
       real,    dimension(MDIM),   intent(IN) :: del
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC
       logical, intent(IN) :: halfTimeAdvance
     end subroutine hy_uhd_staggeredDivb
  end interface



  interface
     subroutine hy_uhd_HLLD(dir,Vm,Vp,Fstar,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM), intent(OUT):: Fstar
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_HLLD
  end interface



  interface
     subroutine hy_uhd_getElectricFields( blockID,blkLimits,blkLimitsGC,del,flx,fly,flz)
       implicit none
       integer, intent(IN)  :: blockID
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
     end subroutine hy_uhd_getElectricFields
  end interface


  interface
     subroutine hy_uhd_getCurrents( &
        blockID, range_switch, blkLimits,datasize, del, Jp, Jm, mode_switch, &
        inDataPtr)
        implicit none
        ! Arguments:
        integer,intent(in) :: blockID, range_switch
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
        real,optional,POINTER_INTENT_IN :: inDataPtr(:,:,:,:)
     end subroutine hy_uhd_getCurrents
  end interface



  interface
     subroutine hy_uhd_getFluxDeriv( ix,iy,iz,blkLimitsGC,&
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
     end subroutine hy_uhd_getFluxDeriv
  end interface



  interface
     Subroutine hy_uhd_addBiermannBatteryTerms(blockID,blkLimitsGC,ix,iy,iz,Flux,sweepDir)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
       integer, INTENT(IN) :: sweepDir
     end Subroutine hy_uhd_addBiermannBatteryTerms
  end interface
  
  interface
     subroutine hy_uhd_biermannSource(blockCount, blockList, dt)
       implicit none
       integer, INTENT(IN) :: blockCount  
       integer, INTENT(IN), dimension(blockCount) :: blockList
       real,    INTENT(IN) :: dt
     end subroutine hy_uhd_biermannSource
  end interface


#endif 
!endif #ifdef FLASH_UNSPLIT_MHD
#endif 
!endif #ifdef FLASH_HYDRO_UNSPLIT

End Module hy_simpleInterface

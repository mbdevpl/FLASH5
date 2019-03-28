!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_Unsplit/hy_llfUnsplit
!!
!! NAME
!!
!!  hy_llfUnsplit
!!
!! SYNOPSIS
!!
!!  call hy_llfUnsplit( integer (IN) :: blockCount,
!!                      integer (IN) :: blockList(blockCount),
!!                      real    (IN) :: dt,
!!                      real    (IN) :: dtOld  )
!!
!! DESCRIPTION
!!
!!  Performs Hydro update in a directionally unsplit fashion over a set
!!  of blocks.
!!  dt gives the timestep over which to advance, and timeEndAdv gives the
!!  simulation time at the end of the update.
!! 
!!  This routine performs a guardcell fill and for each block: 
!!   - applies an Eos call to the guard cells (!DEV: only if/where necessary?)
!!   - computes fluxes
!!   - update all the cell values from the fluxes.
!!   - and finally, we apply an Eos call to the block (interiors).
!!
!!
!!
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList  -  array holding local IDs of blocks on which to advance
!!  dt         -  timestep
!!  dtOld      -  old timestep (not used here)
!!
!! NOTES
!!
!!  This is a simple demo version. Some of the numerous limitiations, compared
!!  to the more serious Hydro implementations available in FLASH:
!!
!!  * Flux correction is not implemented.
!!  * No reconstruction (and thus no limiting) of variables; only HLL as "Riemann solver".
!!  * No support for advecting mass fractions / abundances or other mass scalar
!!    variables.
!!  * No support for non-Cartesian geometries. If this works for any of them,
!!    it should be considered an accident.
!!  * No support for MHD, or anything related to magnetic fields.
!!  * No support for gravity or other body forces or source terms.
!!  * No support for flux-based diffusive terms.
!!  * No artificial viscosity term.
!!  * No support for eintSwitch .NE. 0.0.
!!
!! HISTORY
!!
!!  June  2013  - created KW, outer structure derived from hy_uhd_unsplit.F90 (Dongwook)
!!***

!!REORDER(4): U, fl[xyz]

#ifdef DEBUG_ALL
#define DEBUG_UHD
#endif
#define DEBUG_GRID_GCMASK

Subroutine hy_llfUnsplit ( tileLimits, Uin, plo, Uout, del, dt )

  use Grid_interface, ONLY : Grid_genGetBlkPtr,         &
                             Grid_genReleaseBlkPtr

  use GridTilingModule, ONLY: Grid_tilingContext_t, &
                              Grid_startLoopTiling, &
                              Grid_addToLoopTiling, &
                              Grid_endLoopTiling,   &
                              Grid_getBlkVarPtrs,   &
                              Grid_releaseBlkVarPtrs, &
                              Grid_getTileVarPtrs,   &
                              Grid_releaseTileVarPtrs

#include "Flash.h"

  use Hydro_data, ONLY : hy_fluxCorrect,      &
                         hy_gref,             &
                         hy_useGravity,       &
                         hy_order,            &
                         hy_gcMaskSize,       &
                         hy_gcMask,           &
                         hy_unsplitEosMode,   &
                         hy_eosModeAfter,     &
                         hy_useGravHalfUpdate,&
                         hy_useGravPotUpdate, &
                         hy_gravConsv,        &
                         hy_updateHydroFluxes,&
                         hy_geometry,         &
                         hy_fluxCorVars,      &
                         hy_threadTileList
#ifdef FLASH_USM_MHD
  use Hydro_data, ONLY : hy_E_upwind
#endif

  use Driver_interface, ONLY : Driver_abortFlash

  use Eos_interface, ONLY : Eos_wrapped

  use Logfile_interface, ONLY : Logfile_stampVarMask

!!$  use Timers_interface, ONLY : Timers_start, Timers_stop


  implicit none

#include "constants.h"
#include "Eos.h"
#include "UHD.h"


  !! ---- Argument List ----------------------------------
  integer, intent(IN)  :: tileLimits(LOW:HIGH, 1:MDIM)
  integer, intent(IN)  :: plo(*)
  real,    intent(IN)  :: UIN(plo(1):,plo(2):,plo(3):,plo(4):)  !CAPITALIZATION INTENTIONAL!
  real,    intent(OUT) :: UOUT(plo(1):,plo(2):,plo(3):,plo(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(IN)  :: del(1:MDIM)
  real,    intent(IN)  :: dt
  !! -----------------------------------------------------

!!$  integer, dimension(MDIM) :: datasize
  integer, dimension(LOW:HIGH,MDIM) :: tileLimits
  integer :: ib, i,j,k,blockID
  integer :: is,js,ks
  integer :: ix,iy,iz
  integer :: iL,iR, jL, jR, kL, kR
  real, dimension(MDIM) :: del
  real :: dtdx, dtdy, dtdz, vn, invNewDens
  real :: sMax
  real :: c, a2Avg, cAvg
  real :: vL, vR
  logical :: gcMask(hy_gcMaskSize)
  integer, dimension(2,MDIM) :: eosRange
  integer :: t, tileID

  real, pointer, dimension(:,:,:,:)   :: faceX, faceY, faceZ, auxC

  real, pointer, dimension(:,:,:,:) :: Uin,Uout

  type(Grid_tilingContext_t),pointer :: tilingCtx
  integer :: tileCount
  integer :: tileList(1024)

#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif

  integer,parameter,dimension(HY_VARINUM4) :: &
       outVarList=(/DENS_VAR,&
       VELX_VAR,VELY_VAR,VELZ_VAR,&
       PRES_VAR,&
       GAMC_VAR,GAME_VAR,EINT_VAR,&
       ENER_VAR/)


  !! End of data declaration ***********************************************

#ifdef FLASH_GRID_PARAMESH2
  call Driver_abortFlash("The unsplit Hydro solver only works with PARAMESH 3 or 4!")
#endif


#ifdef FLASH_GRID_PARAMESH3OR4
  if (hy_fluxCorrect) then
     call Driver_abortFlash("hy_llfUnsplit: flux correction is not implemented!")
  end if
#endif

  if (hy_useGravity) then
     call Driver_abortFlash("hy_llfUnsplit: support for gravity not implemented!")
  end if

  if (.NOT.hy_updateHydroFluxes) then
     return
  end if

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(hy_gcMask, .FALSE., '[hy_llfUnsplit]', 'gcNeed')
  end if
#endif
  !! Guardcell filling routine
!!$  call Grid_fillGuardCells(CENTER,ALLDIR,&
!!$       maskSize=hy_gcMaskSize, mask=hy_gcMask,makeMaskConsistent=.true.,doLogMask=.NOT.gcMaskLogged)


  !! ***************************************************************************
  !! There is only one overall loop in this simplified advancement             *
  !! ***************************************************************************
  !! Loop over the blocks


!!$  call Grid_startLoopTiling(CENTER,tilingCtx,&
!!$                  useInPtr=.TRUE.,  useInVarPtrs=.FALSE.,    &
!!$                  useOutPtr=.TRUE., useOutVarPtrs=.FALSE., outVarList=outVarList, &
!!$                  nAuxVars=1, nAuxGuard=1, blockList=blockList(1:blockCount))
!!$  call Grid_addToLoopTiling(tilingCtx,SCRATCH_FACEX,&
!!$                  useInPtr=.FALSE.,  useInVarPtrs=.FALSE.,    &
!!$                  useOutPtr=.FALSE., useOutVarPtrs=.FALSE., nAuxVars=5, nAuxGuard=0)
!!$  if (NDIM > 1) then
!!$     call Grid_addToLoopTiling(tilingCtx,SCRATCH_FACEY,&
!!$                  useInPtr=.FALSE.,  useInVarPtrs=.FALSE.,    &
!!$                  useOutPtr=.FALSE., useOutVarPtrs=.FALSE., nAuxVars=5, nAuxGuard=0)
!!$  end if
!!$  if (NDIM > 2) then
!!$     call Grid_addToLoopTiling(tilingCtx,SCRATCH_FACEZ,&
!!$                  useInPtr=.FALSE.,  useInVarPtrs=.FALSE.,    &
!!$                  useOutPtr=.FALSE., useOutVarPtrs=.FALSE., nAuxVars=5, nAuxGuard=0)
!!$  end if

!!$  do ib=1,blockCount
!!$
!!$     blockID = blockList(ib)
!!$
!!$     call Grid_getDeltas(blockID,del)

     dtdx = dt / del(IAXIS)
     if (NDIM > 1) dtdy = dt / del(JAXIS)
     if (NDIM > 2) dtdz = dt / del(KAXIS)

!!$     call Grid_getListOfTiles(blockID, tileList,tileCount)


!!$     !$omp parallel if (hy_threadTileList) &
!!$     !$omp default(none) &
!!$     !$omp firstprivate(dtdx,dtdy,dtdz,blockID) &
!!$     !$omp private(i,del,tileLimits,tileID,&
!!$     !$omp faceX,faceY,faceZ,Uin,Uout,auxC,&
!!$     !$omp c,sMax,a2Avg,cAvg,vn,is,il,ir,vl,vr,js,jl,jr,ks,kl,kr,invNewDens) &
!!$     !$omp shared(tilingCtx,blockCount,blockList,tileCount,tileList,&
!!$     !$omp hy_unsplitEosMode,hy_useGravity,hy_gref,hy_fluxCorrect,&
!!$     !$omp hy_updateHydroFluxes,hy_eosModeAfter,hy_useGravHalfUpdate,&
!!$     !$omp hy_useGravPotUpdate,hy_geometry,hy_fluxCorVars)

     !$omp do schedule(static)
     do t=1,tileCount
        tileID = tileList(t)

     !! NO call to Eos for guardcell regions - simplified! --------------------
     !! End of Eos call for guardcell regions ----------------------------------

        ! Note: Not handling gravity.



        call Grid_getTileVarPtrs(tileID,gridDataStruct=CENTER, &
             inPtr=Uin, nInGuard=(/1,1,1/), &
             outPtr=Uout,&
             outLimits=tileLimits, &
             dataInit=GRID_DATAINIT_NEGINFINITY,&
             auxPtr=auxC,tilingContext=tilingCtx)


        call Grid_getBlkVarPtrs(tileID,gridDataStruct=SCRATCH_FACEX,auxPtr=faceX,tilingContext=tilingCtx)
        if (NDIM > 1) then 
           call Grid_getBlkVarPtrs(tileID,gridDataStruct=SCRATCH_FACEY,auxPtr=faceY,tilingContext=tilingCtx)
        end if
        if (NDIM > 2) then 
           call Grid_getBlkVarPtrs(tileID,gridDataStruct=SCRATCH_FACEZ,auxPtr=faceZ,tilingContext=tilingCtx)
        end if


     !! ************************************************************************
     !! Calculate Riemann (interface) states

     !  No equivalent really to  call hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del)

        do k = tileLimits(LOW,KAXIS)-K3D,tileLimits(HIGH,KAXIS)+K3D
           do j = tileLimits(LOW,JAXIS)-K2D,tileLimits(HIGH,JAXIS)+K2D
              do i = tileLimits(LOW,IAXIS)-1,tileLimits(HIGH,IAXIS)+1
                 c = sqrt(Uin(GAMC_VAR,i,j,k)*Uin(PRES_VAR,i,j,k)/Uin(DENS_VAR,i,j,k))
                 auxC(1, i,j,k) = c
              end do
           end do
        end do

     !! ************************************************************************
     !! Calculate Godunov fluxes

     !  instead of  call hy_uhd_getFaceFlux(blockID,blkLimits,blkLimitsGC,datasize,del, ...)

        do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
           do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
              do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)+1
                 a2Avg = 0.5* (Uin(GAMC_VAR,i-1,j,k)+Uin(GAMC_VAR,i,j,k)) * &
                      (Uin(PRES_VAR,i-1,j,k)+Uin(PRES_VAR,i,j,k)) / &
                      (Uin(DENS_VAR,i-1,j,k)+Uin(DENS_VAR,i,j,k))
                 cAvg = sqrt(a2Avg)
                 iL = i-1; iR=i
                 vL = Uin(VELX_VAR,iL,j,k);  vR = Uin(VELX_VAR,iR,j,k)
                 vn = (vL+vR) * 0.5
                 sMax = max(abs(vL) + auxC(1, i-1,j,k), &
                      abs(vn) + cAvg, &
                      abs(vR) + auxC(1, i  ,j,k))

                 faceX(HY_DENS_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,iL,j,k) + vR * Uin(DENS_VAR,iR,j,k) &
                      - sMax*(Uin(DENS_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k))) * 0.5
                 faceX(HY_XMOM_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,iL,j,k)*Uin(VELX_VAR,iL,j,k) &
                      + vR * Uin(DENS_VAR,iR,j,k)*Uin(VELX_VAR,iR,j,k) &
                      - sMax*(Uin(DENS_VAR,iR,j,k)*Uin(VELX_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(VELX_VAR,iL,j,k)) )* 0.5
                 faceX(HY_XMOM_FLUX,i,j,k) = faceX(HY_XMOM_FLUX,i,j,k) &
                      + (Uin(PRES_VAR,iL,j,k) + Uin(PRES_VAR,iR,j,k)) * 0.5
                 faceX(HY_YMOM_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,iL,j,k)*Uin(VELY_VAR,iL,j,k) &
                      + vR * Uin(DENS_VAR,iR,j,k)*Uin(VELY_VAR,iR,j,k) &
                      - sMax*(Uin(DENS_VAR,iR,j,k)*Uin(VELY_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(VELY_VAR,iL,j,k)) )* 0.5
                 faceX(HY_ZMOM_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,iL,j,k)*Uin(VELZ_VAR,iL,j,k) &
                      + vR * Uin(DENS_VAR,iR,j,k)*Uin(VELZ_VAR,iR,j,k) &
                      - sMax*(Uin(DENS_VAR,iR,j,k)*Uin(VELZ_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(VELZ_VAR,iL,j,k)) )* 0.5
                 faceX(HY_ENER_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,iL,j,k)*Uin(ENER_VAR,iL,j,k) &
                      + vR * Uin(DENS_VAR,iR,j,k)*Uin(ENER_VAR,iR,j,k) &
                      - sMax*(Uin(DENS_VAR,iR,j,k)*Uin(ENER_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(ENER_VAR,iL,j,k)))* 0.5
                 faceX(HY_ENER_FLUX,i,j,k) = faceX(HY_ENER_FLUX,i,j,k) &
                      + (vL * Uin(PRES_VAR,iL,j,k) + vR * Uin(PRES_VAR,iR,j,k)) * 0.5

                 faceX(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) = faceX(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) * dtdx
              end do
           end do
        end do
#if NDIM > 1
        do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
           do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)+1
              do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
                 a2Avg = 0.5* (Uin(GAMC_VAR,i,j-1,k)+Uin(GAMC_VAR,i,j,k)) * &
                      (Uin(PRES_VAR,i,j-1,k)+Uin(PRES_VAR,i,j,k)) / &
                      (Uin(DENS_VAR,i,j-1,k)+Uin(DENS_VAR,i,j,k))
                 cAvg = sqrt(a2Avg)
                 jL = j-1; jR=j
                 vL = Uin(VELY_VAR,i,jL,k);  vR = Uin(VELY_VAR,i,jR,k)
                 vn = (vL+vR) * 0.5
                 sMax = max(abs(vL) + auxC(1, i,j-1,k), &
                      abs(vn) + cAvg, &
                      abs(vR) + auxC(1, i,j  ,k))

                 faceY(HY_DENS_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,i,jL,k) + vR * Uin(DENS_VAR,i,jR,k) &
                      - sMax*(Uin(DENS_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k))) * 0.5
                 faceY(HY_XMOM_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,i,jL,k)*Uin(VELX_VAR,i,jL,k) &
                      + vR * Uin(DENS_VAR,i,jR,k)*Uin(VELX_VAR,i,jR,k) &
                      - sMax*(Uin(DENS_VAR,i,jR,k)*Uin(VELX_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(VELX_VAR,i,jL,k)) )* 0.5
                 faceY(HY_YMOM_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,i,jL,k)*Uin(VELY_VAR,i,jL,k) &
                      + vR * Uin(DENS_VAR,i,jR,k)*Uin(VELY_VAR,i,jR,k) &
                      - sMax*(Uin(DENS_VAR,i,jR,k)*Uin(VELY_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(VELY_VAR,i,jL,k)) )* 0.5
                 faceY(HY_YMOM_FLUX,i,j,k) = faceY(HY_YMOM_FLUX,i,j,k) &
                      + (Uin(PRES_VAR,i,jL,k) + Uin(PRES_VAR,i,jR,k)) * 0.5
                 faceY(HY_ZMOM_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,i,jL,k)*Uin(VELZ_VAR,i,jL,k) &
                      + vR * Uin(DENS_VAR,i,jR,k)*Uin(VELZ_VAR,i,jR,k) &
                      - sMax*(Uin(DENS_VAR,i,jR,k)*Uin(VELZ_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(VELZ_VAR,i,jL,k)) )* 0.5
                 faceY(HY_ENER_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,i,jL,k)*Uin(ENER_VAR,i,jL,k) &
                      + vR * Uin(DENS_VAR,i,jR,k)*Uin(ENER_VAR,i,jR,k) &
                      - sMax*(Uin(DENS_VAR,i,jR,k)*Uin(ENER_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(ENER_VAR,i,jL,k)))* 0.5
                 faceY(HY_ENER_FLUX,i,j,k) = faceY(HY_ENER_FLUX,i,j,k) &
                      + (vL * Uin(PRES_VAR,i,jL,k) + vR * Uin(PRES_VAR,i,jR,k)) * 0.5

                 faceY(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) = faceY(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) * dtdy
              end do
           end do
        end do
#endif
#if NDIM > 2
        do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)+1
           do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
              do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
                 a2Avg = 0.5* (Uin(GAMC_VAR,i,j,k-1)+Uin(GAMC_VAR,i,j,k)) * &
                      (Uin(PRES_VAR,i,j,k-1)+Uin(PRES_VAR,i,j,k)) / &
                      (Uin(DENS_VAR,i,j,k-1)+Uin(DENS_VAR,i,j,k))
                 cAvg = sqrt(a2Avg)
                 kL = k-1; kR=i
                 vL = Uin(VELZ_VAR,i,j,kL);  vR = Uin(VELZ_VAR,i,j,kR)
                 vn = (vL+vR) * 0.5
                 sMax = max(abs(vL) + auxC(1, i,j,k-1), &
                      abs(vn) + cAvg, &
                      abs(vR) + auxC(1, i,j,k  ))

                 faceZ(HY_DENS_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,i,j,kL) + vR * Uin(DENS_VAR,i,j,kR) &
                      - sMax*(Uin(DENS_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL))) * 0.5
                 faceZ(HY_XMOM_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,i,j,kL)*Uin(VELX_VAR,i,j,kL) &
                      + vR * Uin(DENS_VAR,i,j,kR)*Uin(VELX_VAR,i,j,kR) &
                      - sMax*(Uin(DENS_VAR,i,j,kR)*Uin(VELX_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(VELX_VAR,i,j,kL)) )* 0.5
                 faceZ(HY_YMOM_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,i,j,kL)*Uin(VELY_VAR,i,j,kL) &
                      + vR * Uin(DENS_VAR,i,j,kR)*Uin(VELY_VAR,i,j,kR) &
                      - sMax*(Uin(DENS_VAR,i,j,kR)*Uin(VELY_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(VELY_VAR,i,j,kL)) )* 0.5
                 faceZ(HY_ZMOM_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,i,j,kL)*Uin(VELZ_VAR,i,j,kL) &
                      + vR * Uin(DENS_VAR,i,j,kR)*Uin(VELZ_VAR,i,j,kR) &
                      - sMax*(Uin(DENS_VAR,i,j,kR)*Uin(VELZ_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(VELZ_VAR,i,j,kL)) )* 0.5
                 faceZ(HY_ZMOM_FLUX,i,j,k) = faceZ(HY_ZMOM_FLUX,i,j,k) &
                      + (Uin(PRES_VAR,i,j,kL) + Uin(PRES_VAR,i,j,kR)) * 0.5
                 faceZ(HY_ENER_FLUX,i,j,k) = (  vL * Uin(DENS_VAR,i,j,kL)*Uin(ENER_VAR,i,j,kL) &
                      + vR * Uin(DENS_VAR,i,j,kR)*Uin(ENER_VAR,i,j,kR) &
                      - sMax*(Uin(DENS_VAR,i,j,kR)*Uin(ENER_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(ENER_VAR,i,j,kL)))* 0.5
                 faceZ(HY_ENER_FLUX,i,j,k) = faceZ(HY_ENER_FLUX,i,j,k) &
                      + (vL * Uin(PRES_VAR,i,j,kL) + vR * Uin(PRES_VAR,i,j,kR)) * 0.5

                 faceZ(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) = faceZ(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) * dtdz
              end do
           end do
        end do
#endif

        deallocate(auxC)


     !! ************************************************************************
     !! Unsplit update for conservative variables from n to n+1 time step

     !  instead of  call hy_hllUnsplitUpdate(blockID,dt,dtOld,del,datasize,blkLimits, ...)

        do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
           do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
              do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
                 Uout(VELX_VAR,i,j,k) = Uin(VELX_VAR,i,j,k) * Uin(DENS_VAR,i,j,k)
                 Uout(VELY_VAR,i,j,k) = Uin(VELY_VAR,i,j,k) * Uin(DENS_VAR,i,j,k)
                 Uout(VELZ_VAR,i,j,k) = Uin(VELX_VAR,i,j,k) * Uin(DENS_VAR,i,j,k)
                 Uout(ENER_VAR,i,j,k) = Uin(ENER_VAR,i,j,k) * Uin(DENS_VAR,i,j,k)
                 Uout(DENS_VAR,i,j,k) = Uin(DENS_VAR,i,j,k) + faceX(HY_DENS_FLUX,i,j,k) - faceX(HY_DENS_FLUX,i+1,j,k)
                 if (NDIM > 1) Uout(DENS_VAR,i,j,k) = Uout(DENS_VAR,i,j,k) + faceY(HY_DENS_FLUX,i,j,k) - faceY(HY_DENS_FLUX,i,j+1,k)
                 if (NDIM > 2) Uout(DENS_VAR,i,j,k) = Uout(DENS_VAR,i,j,k) + faceZ(HY_DENS_FLUX,i,j,k) - faceZ(HY_DENS_FLUX,i,j,k+1)

              end do
           end do
        end do
        do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
           do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
              do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
                 Uout(VELX_VAR,i,j,k) = Uout(VELX_VAR,i,j,k) + faceX(HY_XMOM_FLUX,i,j,k) - faceX(HY_XMOM_FLUX,i+1,j,k)
                 Uout(VELY_VAR,i,j,k) = Uout(VELY_VAR,i,j,k) + faceX(HY_YMOM_FLUX,i,j,k) - faceX(HY_YMOM_FLUX,i+1,j,k)
                 Uout(VELZ_VAR,i,j,k) = Uout(VELZ_VAR,i,j,k) + faceX(HY_ZMOM_FLUX,i,j,k) - faceX(HY_ZMOM_FLUX,i+1,j,k)
                 Uout(ENER_VAR,i,j,k) = Uout(ENER_VAR,i,j,k) + faceX(HY_ENER_FLUX,i,j,k) - faceX(HY_ENER_FLUX,i+1,j,k)
              end do
           end do
        end do
#if NDIM > 1
        do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
           do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
              do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
                 Uout(VELX_VAR,i,j,k) = Uout(VELX_VAR,i,j,k) + faceY(HY_XMOM_FLUX,i,j,k) - faceY(HY_XMOM_FLUX,i,j+1,k)
                 Uout(VELY_VAR,i,j,k) = Uout(VELY_VAR,i,j,k) + faceY(HY_YMOM_FLUX,i,j,k) - faceY(HY_YMOM_FLUX,i,j+1,k)
                 Uout(VELZ_VAR,i,j,k) = Uout(VELZ_VAR,i,j,k) + faceY(HY_ZMOM_FLUX,i,j,k) - faceY(HY_ZMOM_FLUX,i,j+1,k)
                 Uout(ENER_VAR,i,j,k) = Uout(ENER_VAR,i,j,k) + faceY(HY_ENER_FLUX,i,j,k) - faceY(HY_ENER_FLUX,i,j+1,k)
              end do
           end do
        end do
#endif
#if NDIM > 2
        do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
           do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
              do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
                 Uout(VELX_VAR,i,j,k) = Uout(VELX_VAR,i,j,k) + faceZ(HY_XMOM_FLUX,i,j,k) - faceZ(HY_XMOM_FLUX,i,j,k+1)
                 Uout(VELY_VAR,i,j,k) = Uout(VELY_VAR,i,j,k) + faceZ(HY_YMOM_FLUX,i,j,k) - faceZ(HY_YMOM_FLUX,i,j,k+1)
                 Uout(VELZ_VAR,i,j,k) = Uout(VELZ_VAR,i,j,k) + faceZ(HY_ZMOM_FLUX,i,j,k) - faceZ(HY_ZMOM_FLUX,i,j,k+1)
                 Uout(ENER_VAR,i,j,k) = Uout(ENER_VAR,i,j,k) + faceZ(HY_ENER_FLUX,i,j,k) - faceZ(HY_ENER_FLUX,i,j,k+1)
              end do
           end do
        end do
#endif
        do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
           do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
              do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
                 invNewDens = 1.0 / Uout(DENS_VAR,i,j,k)
                 Uout(VELX_VAR,i,j,k) = Uout(VELX_VAR,i,j,k) * invNewDens
                 Uout(VELY_VAR,i,j,k) = Uout(VELY_VAR,i,j,k) * invNewDens
                 Uout(VELZ_VAR,i,j,k) = Uout(VELZ_VAR,i,j,k) * invNewDens
                 Uout(ENER_VAR,i,j,k) = Uout(ENER_VAR,i,j,k) * invNewDens
              end do
           end do
        end do


     !! Correct energy if necessary
     !  instead of  call hy_uhd_energyFix(blockID,blkLimits,dt,del,hy_unsplitEosMode)

#ifdef EINT_VAR
        do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
           do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
              do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
                 invNewDens = 1.0 / Uout(DENS_VAR,i,j,k)
                 Uout(EINT_VAR,i,j,k) = Uout(ENER_VAR,i,j,k) &
                      - 0.5 * dot_product(Uout(VELX_VAR:VELZ_VAR,i,j,k),Uout(VELX_VAR:VELZ_VAR,i,j,k))
              end do
           end do
        end do
#endif

     
        deallocate(faceX)
        if (NDIM > 1) then 
           deallocate(faceY)
        end if
        if (NDIM > 2) then 
           deallocate(faceZ)
        end if

        !! Call to Eos - note this is a variant where we pass a buffer not a blockID.
        call Eos_wrapped(hy_eosModeAfter, tileLimits, Uout)

        call Grid_releaseTileVarPtrs(tileID,gridDataStruct=CENTER, &
             inPtr=Uin, &
             outPtr=Uout,&
             tilingContext=tilingCtx)
     end do
     !$omp end do

     !$omp end parallel

  end do
  !! End of leaf block do-loop - no flux conserve call

  call Grid_endLoopTiling(tilingCtx)

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
#endif

End Subroutine hy_llfUnsplit

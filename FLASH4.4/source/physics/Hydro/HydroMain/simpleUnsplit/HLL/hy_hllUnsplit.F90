!!****if* source/physics/Hydro/HydroMain/simpleUnsplit/HLL/hy_hllUnsplit
!!
!! NAME
!!
!!  hy_hllUnsplit
!!
!! SYNOPSIS
!!
!!  call hy_hllUnsplit( integer (IN) :: blockCount,
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

! Note: the following arrays need to be spelled exactly like this in the code below,
!       preserving case.
!!REORDER(4): Uin, Uout, face[XYZ], auxC

#ifdef DEBUG_ALL
#define DEBUG_UHD
#endif

#include "Flash.h"
#include "constants.h"
#include "Eos.h"
#include "UHD.h"

Subroutine hy_hllUnsplit ( tileLimits, Uin, plo, Uout, del, dt )
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
                         hy_threadBlockList
#ifdef FLASH_USM_MHD
  use Hydro_data, ONLY : hy_E_upwind
#endif

  use Driver_interface,  ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampVarMask

  implicit none

  !! ---- Argument List ----------------------------------
  integer, intent(IN)  :: tileLimits(LOW:HIGH, 1:MDIM)
  integer, intent(IN)  :: plo(*)
  real,    intent(IN)  :: UIN(plo(1):,plo(2):,plo(3):,plo(4):)  !CAPITALIZATION INTENTIONAL!
  real,    intent(OUT) :: UOUT(plo(1):,plo(2):,plo(3):,plo(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(IN)  :: del(1:MDIM)
  real,    intent(IN)  :: dt
  !! -----------------------------------------------------

  integer :: ib, i,j,k,blockID
  integer :: is,js,ks
  integer :: ix,iy,iz
  integer :: iL,iR, jL, jR, kL, kR
  real :: dtdx, dtdy, dtdz, vn, invNewDens
  real :: c, sL, sR
  real :: sRsL, vL, vR
  integer, dimension(2,MDIM) :: eosRange
  integer :: t, tileID

  real, pointer, dimension(:,:,:,:) :: faceX, faceY, faceZ, auxC

  integer :: tileCount
  integer :: tileList(1024)

  integer,parameter,dimension(HY_VARINUM4) :: &
       outVarList=(/DENS_VAR,&
       VELX_VAR,VELY_VAR,VELZ_VAR,&
       PRES_VAR,&
       GAMC_VAR,GAME_VAR,EINT_VAR,&
       ENER_VAR/)

  nullify(faceX)
  nullify(faceY)
  nullify(faceZ)
  nullify(auxC)

  !! End of data declaration ***********************************************
#ifdef DEBUG_UHD
98 format(A4,'(',I3,':   ,',   I3,':   ,',   I3,':   ,',   I3,':   )')
99 format(A4,'(',I3,':',I3,',',I3,':',I3,',',I3,':',I3,',',I3,':',I3,')')
  print *, "plo" ,plo(1:MDIM+1)
  print 98,"Uin" ,(plo(i),i=1,4)
  print 99,"Uin" ,(lbound(Uin ,i),ubound(Uin ,i),i=1,4)
  print 99,"Uout",(lbound(Uout,i),ubound(Uout,i),i=1,4)
  print*,'tileLim:',tileLimits
#endif

#ifdef FLASH_GRID_PARAMESH2
  call Driver_abortFlash("The unsplit Hydro solver only works with PARAMESH 3 or 4!")
#endif


#ifdef FLASH_GRID_PARAMESH3OR4
  if (hy_fluxCorrect) then
     call Driver_abortFlash("hy_hllUnsplit: flux correction is not implemented!")
  end if
#endif

  if (hy_useGravity) then
     call Driver_abortFlash("hy_hllUnsplit: support for gravity not implemented!")
  end if

  if (.NOT.hy_updateHydroFluxes) then
     return
  end if

  !! ***************************************************************************
  !! There is only one overall loop in this simplified advancement             *
  !! ***************************************************************************
  !! Loop over the blocks
     
  dtdx = dt / del(IAXIS)
  if (NDIM > 1) dtdy = dt / del(JAXIS)
  if (NDIM > 2) dtdz = dt / del(KAXIS)

  ! Note: Not handling gravity.

  allocate(auxC(1,tileLimits(LOW,IAXIS)-1  :tileLimits(HIGH,IAXIS)+1  , &
                  tileLimits(LOW,JAXIS)-K2D:tileLimits(HIGH,JAXIS)+K2D, &
                  tileLimits(LOW,KAXIS)-K3D:tileLimits(HIGH,KAXIS)+K3D) )

  allocate(faceX(5,tileLimits(LOW,IAXIS):tileLimits(HIGH,IAXIS)+1  , &
                  tileLimits(LOW,JAXIS):tileLimits(HIGH,JAXIS), &
                  tileLimits(LOW,KAXIS):tileLimits(HIGH,KAXIS)) )

  if (NDIM > 1) then 
     allocate(faceY(5,tileLimits(LOW,IAXIS):tileLimits(HIGH,IAXIS)  , &
                      tileLimits(LOW,JAXIS):tileLimits(HIGH,JAXIS)+1, &
                      tileLimits(LOW,KAXIS):tileLimits(HIGH,KAXIS)) )
  end if
  if (NDIM > 2) then 
     allocate(faceZ(5,tileLimits(LOW,IAXIS):tileLimits(HIGH,IAXIS)  , &
                      tileLimits(LOW,JAXIS):tileLimits(HIGH,JAXIS), &
                      tileLimits(LOW,KAXIS):tileLimits(HIGH,KAXIS)+1) )
  end if

  !! ************************************************************************
  !! Calculate Riemann (interface) states
  !  No equivalent really to  call hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del)

  !! calculate sound speed
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
           sL = min(Uin(VELX_VAR,i-1,j,k)-auxC(1, i-1,j,k), Uin(VELX_VAR,i,j,k)-auxC(1, i,j,k))
           sR = max(Uin(VELX_VAR,i-1,j,k)+auxC(1, i-1,j,k), Uin(VELX_VAR,i,j,k)+auxC(1, i,j,k))
           sRsL = sR - sL
           if (sL > 0.0) then
              vn = Uin(VELX_VAR,i-1,j,k)
              is = i-1
              iL = i-1; iR=i-1
           else if (sR < 0.0) then
              vn = Uin(VELX_VAR,i,j,k)
              is = i
              iL = i; iR=i
           else
              vn = (Uin(VELX_VAR,i-1,j,k)+Uin(VELX_VAR,i,j,k)) * 0.5
              is = i
              iL = i-1; iR=i
              if (vn>0.0) is = is-1
           end if
           vL = Uin(VELX_VAR,iL,j,k);  vR = Uin(VELX_VAR,iR,j,k)
           if (iL==iR) then
              faceX(HY_DENS_FLUX,i,j,k) = vn * Uin(DENS_VAR,is,j,k)
              faceX(HY_XMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,is,j,k) * Uin(VELX_VAR,is,j,k)
              faceX(HY_XMOM_FLUX,i,j,k) = faceX(HY_XMOM_FLUX,i,j,k) &
                   + Uin(PRES_VAR,is,j,k)
              faceX(HY_YMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,is,j,k) * Uin(VELY_VAR,is,j,k)
              faceX(HY_ZMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,is,j,k) * Uin(VELZ_VAR,is,j,k)
              faceX(HY_ENER_FLUX,i,j,k) = vn * Uin(DENS_VAR,is,j,k) * Uin(ENER_VAR,is,j,k)
              faceX(HY_ENER_FLUX,i,j,k) = faceX(HY_ENER_FLUX,i,j,k) &
                   + vn * Uin(PRES_VAR,is,j,k)
           else
              faceX(HY_DENS_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,iL,j,k) - sL * vR * Uin(DENS_VAR,iR,j,k) &
                   + sR*sL*(Uin(DENS_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k))) / sRsL
              faceX(HY_XMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,iL,j,k)*Uin(VELX_VAR,iL,j,k) &
                   - sL * vR * Uin(DENS_VAR,iR,j,k)*Uin(VELX_VAR,iR,j,k) &
                   + sR*sL*(Uin(DENS_VAR,iR,j,k)*Uin(VELX_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(VELX_VAR,iL,j,k)) )/sRsL
              faceX(HY_XMOM_FLUX,i,j,k) = faceX(HY_XMOM_FLUX,i,j,k) &
                   + (sR * Uin(PRES_VAR,iL,j,k) - sL * Uin(PRES_VAR,iR,j,k)) / sRsL
              faceX(HY_YMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,iL,j,k)*Uin(VELY_VAR,iL,j,k) &
                   - sL * vR * Uin(DENS_VAR,iR,j,k)*Uin(VELY_VAR,iR,j,k) &
                   + sR*sL*(Uin(DENS_VAR,iR,j,k)*Uin(VELY_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(VELY_VAR,iL,j,k)) )/sRsL
              faceX(HY_ZMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,iL,j,k)*Uin(VELZ_VAR,iL,j,k) &
                   - sL * vR * Uin(DENS_VAR,iR,j,k)*Uin(VELZ_VAR,iR,j,k) &
                   + sR*sL*(Uin(DENS_VAR,iR,j,k)*Uin(VELZ_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(VELZ_VAR,iL,j,k)) )/sRsL
              faceX(HY_ENER_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,iL,j,k)*Uin(ENER_VAR,iL,j,k) &
                   - sL * vR * Uin(DENS_VAR,iR,j,k)*Uin(ENER_VAR,iR,j,k) &
                   + sR*sL*(Uin(DENS_VAR,iR,j,k)*Uin(ENER_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(ENER_VAR,iL,j,k)))/sRsL
              faceX(HY_ENER_FLUX,i,j,k) = faceX(HY_ENER_FLUX,i,j,k) &
                   + (sR * vL * Uin(PRES_VAR,iL,j,k) - sL * vR * Uin(PRES_VAR,iR,j,k)) / sRsL
           end if
           faceX(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) = faceX(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) * dtdx
        end do
     end do
  end do
#if NDIM > 1
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)+1
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           sL = min(Uin(VELY_VAR,i,j-1,k)-auxC(1, i,j-1,k), Uin(VELY_VAR,i,j,k)-auxC(1, i,j,k))
           sR = max(Uin(VELY_VAR,i,j-1,k)+auxC(1, i,j-1,k), Uin(VELY_VAR,i,j,k)+auxC(1, i,j,k))
           sRsL = sR - sL
           if (sL > 0.0) then
              vn = Uin(VELY_VAR,i,j-1,k)
              js = j-1
              jL = j-1; jR=j-1
           else if (sR < 0.0) then
              vn = Uin(VELY_VAR,i,j,k)
              js = j
              jL = j; jR=j
           else
              vn = (Uin(VELY_VAR,i,j-1,k)+Uin(VELY_VAR,i,j,k)) * 0.5
              js = j
              jL = j-1; jR=j
              if (vn>0.0) js = js-1
           end if
           vL = Uin(VELY_VAR,i,jL,k);  vR = Uin(VELY_VAR,i,jR,k)
           if (jL==jR) then
              faceY(HY_DENS_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,js,k)
              faceY(HY_XMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,js,k) * Uin(VELX_VAR,i,js,k)
              faceY(HY_YMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,js,k) * Uin(VELY_VAR,i,js,k)
              faceY(HY_YMOM_FLUX,i,j,k) = faceY(HY_YMOM_FLUX,i,j,k) &
                   + Uin(PRES_VAR,i,js,k)
              faceY(HY_ZMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,js,k) * Uin(VELZ_VAR,i,js,k)
              faceY(HY_ENER_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,js,k) * Uin(ENER_VAR,i,js,k)
              faceY(HY_ENER_FLUX,i,j,k) = faceY(HY_ENER_FLUX,i,j,k) &
                   + vn * Uin(PRES_VAR,i,js,k)
           else
              faceY(HY_DENS_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,jL,k) - sL * vR * Uin(DENS_VAR,i,jR,k) &
                   + sR*sL*(Uin(DENS_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k))) / sRsL
              faceY(HY_XMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,jL,k)*Uin(VELX_VAR,i,jL,k) &
                   - sL * vR * Uin(DENS_VAR,i,jR,k)*Uin(VELX_VAR,i,jR,k) &
                   + sR*sL*(Uin(DENS_VAR,i,jR,k)*Uin(VELX_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(VELX_VAR,i,jL,k)) )/sRsL
              faceY(HY_YMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,jL,k)*Uin(VELY_VAR,i,jL,k) &
                   - sL * vR * Uin(DENS_VAR,i,jR,k)*Uin(VELY_VAR,i,jR,k) &
                   + sR*sL*(Uin(DENS_VAR,i,jR,k)*Uin(VELY_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(VELY_VAR,i,jL,k)) )/sRsL
              faceY(HY_YMOM_FLUX,i,j,k) = faceY(HY_YMOM_FLUX,i,j,k) &
                   + (sR * Uin(PRES_VAR,i,jL,k) - sL * Uin(PRES_VAR,i,jR,k)) / sRsL
              faceY(HY_ZMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,jL,k)*Uin(VELZ_VAR,i,jL,k) &
                   - sL * vR * Uin(DENS_VAR,i,jR,k)*Uin(VELZ_VAR,i,jR,k) &
                   + sR*sL*(Uin(DENS_VAR,i,jR,k)*Uin(VELZ_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(VELZ_VAR,i,jL,k)) )/sRsL
              faceY(HY_ENER_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,jL,k)*Uin(ENER_VAR,i,jL,k) &
                   - sL * vR * Uin(DENS_VAR,i,jR,k)*Uin(ENER_VAR,i,jR,k) &
                   + sR*sL*(Uin(DENS_VAR,i,jR,k)*Uin(ENER_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(ENER_VAR,i,jL,k)))/sRsL
              faceY(HY_ENER_FLUX,i,j,k) = faceY(HY_ENER_FLUX,i,j,k) &
                   + (sR * vL * Uin(PRES_VAR,i,jL,k) - sL * vR * Uin(PRES_VAR,i,jR,k)) / sRsL
           end if
           faceY(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) = faceY(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) * dtdy
        end do
     end do
  end do
#endif
#if NDIM > 2
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)+1
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           sL = min(Uin(VELZ_VAR,i,j,k-1)-auxC(1, i,j,k-1), Uin(VELZ_VAR,i,j,k)-auxC(1, i,j,k))
           sR = max(Uin(VELZ_VAR,i,j,k-1)+auxC(1, i,j,k-1), Uin(VELZ_VAR,i,j,k)+auxC(1, i,j,k))
           sRsL = sR - sL
           if (sL > 0.0) then
              vn = Uin(VELZ_VAR,i,j,k-1)
              ks = k-1
              kL = k-1; kR = k-1
           else if (sR < 0.0) then
              vn = Uin(VELZ_VAR,i,j,k)
              ks = k
              kL = k; kR = k
           else
              vn = (Uin(VELZ_VAR,i,j,k-1)+Uin(VELZ_VAR,i,j,k)) * 0.5
              ks = k
              kL = k-1; kR = k
              if (vn>0.0) ks = ks-1
           end if
           vL = Uin(VELZ_VAR,i,j,kL);  vR = Uin(VELZ_VAR,i,j,kR)
           if (kL==kR) then
              faceZ(HY_DENS_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,j,ks)
              faceZ(HY_XMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,j,ks) * Uin(VELX_VAR,i,j,ks)
              faceZ(HY_YMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,j,ks) * Uin(VELY_VAR,i,j,ks)
              faceZ(HY_ZMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,j,ks) * Uin(VELZ_VAR,i,j,ks)
              faceZ(HY_ZMOM_FLUX,i,j,k) = faceZ(HY_ZMOM_FLUX,i,j,k) &
                   + Uin(PRES_VAR,i,j,ks)
              faceZ(HY_ENER_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,j,ks) * Uin(ENER_VAR,i,j,ks)
              faceZ(HY_ENER_FLUX,i,j,k) = faceZ(HY_ENER_FLUX,i,j,k) &
                   + vn * Uin(PRES_VAR,i,j,ks)
           else

              faceZ(HY_DENS_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,j,kL) - sL * vR * Uin(DENS_VAR,i,j,kR) &
                   + sR*sL*(Uin(DENS_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL))) / sRsL
              faceZ(HY_XMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,j,kL)*Uin(VELX_VAR,i,j,kL) &
                   - sL * vR * Uin(DENS_VAR,i,j,kR)*Uin(VELX_VAR,i,j,kR) &
                   + sR*sL*(Uin(DENS_VAR,i,j,kR)*Uin(VELX_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(VELX_VAR,i,j,kL)) )/sRsL
              faceZ(HY_YMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,j,kL)*Uin(VELY_VAR,i,j,kL) &
                   - sL * vR * Uin(DENS_VAR,i,j,kR)*Uin(VELY_VAR,i,j,kR) &
                   + sR*sL*(Uin(DENS_VAR,i,j,kR)*Uin(VELY_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(VELY_VAR,i,j,kL)) )/sRsL
              faceZ(HY_ZMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,j,kL)*Uin(VELZ_VAR,i,j,kL) &
                   - sL * vR * Uin(DENS_VAR,i,j,kR)*Uin(VELZ_VAR,i,j,kR) &
                   + sR*sL*(Uin(DENS_VAR,i,j,kR)*Uin(VELZ_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(VELZ_VAR,i,j,kL)) )/sRsL
              faceZ(HY_ZMOM_FLUX,i,j,k) = faceZ(HY_ZMOM_FLUX,i,j,k) &
                   + (sR * Uin(PRES_VAR,i,j,kL) - sL * Uin(PRES_VAR,i,j,kR)) / sRsL
              faceZ(HY_ENER_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,j,kL)*Uin(ENER_VAR,i,j,kL) &
                   - sL * vR * Uin(DENS_VAR,i,j,kR)*Uin(ENER_VAR,i,j,kR) &
                   + sR*sL*(Uin(DENS_VAR,i,j,kR)*Uin(ENER_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(ENER_VAR,i,j,kL)))/sRsL
              faceZ(HY_ENER_FLUX,i,j,k) = faceZ(HY_ENER_FLUX,i,j,k) &
                   + (sR * vL * Uin(PRES_VAR,i,j,kL) - sL * vR * Uin(PRES_VAR,i,j,kR)) / sRsL
           end if
           faceZ(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) = faceZ(HY_DENS_FLUX:HY_ENER_FLUX,i,j,k) * dtdz
        end do
     end do
  end do
#endif

  deallocate(auxC)

  !! ************************************************************************
  !! Unsplit update for conservative variables from n to n+1 time step
  !  instead of  call hy_hllUnsplitUpdate(blockID,dt,dtOld,del,datasize,blkLimits, ...)
  !! this section starts the update
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           Uout(VELX_VAR,i,j,k) = Uin(VELX_VAR,i,j,k) * Uin(DENS_VAR,i,j,k)
           Uout(VELY_VAR,i,j,k) = Uin(VELY_VAR,i,j,k) * Uin(DENS_VAR,i,j,k)
           Uout(VELZ_VAR,i,j,k) = Uin(VELZ_VAR,i,j,k) * Uin(DENS_VAR,i,j,k)
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

End Subroutine hy_hllUnsplit

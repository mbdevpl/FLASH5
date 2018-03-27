!!****if* source/physics/Hydro/HydroMain/unsplit/hy_unitConvert
!!
!! NAME
!!
!!  hy_unitConvert
!!
!! SYNOPSIS
!!
!!  hy_unitConvert( integer (IN) :: blockID, 
!!                      integer (IN) :: convertDir)                    
!!
!! DESCRIPTION
!!
!!  This routine converts physical unit of measure (CGS, SI, or NONE) 
!!  for the unsplit MHD/Hydro units.
!!
!! ARGUMENTS
!!
!!  blockID    - local block ID
!!  convertDir - a direction of conversion
!!
!!***

!!REORDER(4): U,B[xyz]

Subroutine hy_unitConvert(U,blkLimitsGC,convertDir)

  use Hydro_data,     ONLY : hy_dref, hy_eref, hy_pref, &
                             hy_vref, hy_bref
  use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Argument list -------------------------
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
  integer, intent(IN) :: convertDir
  real, dimension(:,:,:,:) :: U
  !! ---------------------------------------

  integer :: ib,ie,jb,je,kb,ke

  ib=blkLimitsGC(LOW,IAXIS)
  ie=blkLimitsGC(HIGH,IAXIS)  
  jb=blkLimitsGC(LOW,JAXIS)
  je=blkLimitsGC(HIGH,JAXIS)  
  kb=blkLimitsGC(LOW,KAXIS)
  ke=blkLimitsGC(HIGH,KAXIS)  
#ifdef FLASH_USM_MHD
#if NFACE_VARS > 0
#if NDIM > 1
  real, pointer, dimension(:,:,:,:) :: Bx,By,Bz
#endif
#endif
#endif


!!$#ifdef FLASH_USM_MHD
!!$#if NFACE_VARS > 0
!!$#if NDIM > 1
!!$  call Grid_getBlkPtr(blockID,Bx,FACEX)
!!$  call Grid_getBlkPtr(blockID,By,FACEY)
!!$  if (NDIM == 3) call Grid_getBlkPtr(blockID,Bz,FACEZ)
!!$#endif
!!$#endif
!!$#endif

  select case(convertDir)
  case(FWDCONVERT)
     U(DENS_VAR,ib:ie,jb:je,kb:ke) = U(DENS_VAR,ib:ie,jb:je,kb:ke)/hy_dref
     U(ENER_VAR,ib:ie,jb:je,kb:ke) = U(ENER_VAR,ib:ie,jb:je,kb:ke)/hy_eref
     U(PRES_VAR,ib:ie,jb:je,kb:ke) = U(PRES_VAR,ib:ie,jb:je,kb:ke)/hy_pref
     U(VELX_VAR:VELZ_VAR,ib:ie,jb:je,kb:ke) = U(VELX_VAR:VELZ_VAR,ib:ie,jb:je,kb:ke)/hy_vref

!!$#ifdef FLASH_USM_MHD
!!$     U(MAGX_VAR:MAGZ_VAR,ib:ie,jb:je,kb:ke) = U(MAGX_VAR:MAGZ_VAR,ib:ie,jb:je,kb:ke)/hy_bref
!!$#if NFACE_VARS > 0
!!$#if NDIM > 1
!!$     if (NDIM >= 2) then
!!$        Bx(MAG_FACE_VAR,ib:ie,jb:je,kb:ke) = Bx(MAG_FACE_VAR,ib:ie,jb:je,kb:ke)/hy_bref
!!$        By(MAG_FACE_VAR,ib:ie,jb:je,kb:ke) = By(MAG_FACE_VAR,ib:ie,jb:je,kb:ke)/hy_bref
!!$        if (NDIM==3) Bz(MAG_FACE_VAR,ib:ie,jb:je,kb:ke) = Bz(MAG_FACE_VAR,ib:ie,jb:je,kb:ke)/hy_bref
!!$     endif
!!$#endif
!!$#endif
!!$#endif

  case(BWDCONVERT)
     U(DENS_VAR,ib:ie,jb:je,kb:ke) = U(DENS_VAR,ib:ie,jb:je,kb:ke)*hy_dref
     U(ENER_VAR,ib:ie,jb:je,kb:ke) = U(ENER_VAR,ib:ie,jb:je,kb:ke)*hy_eref
     U(VELX_VAR:VELZ_VAR,ib:ie,jb:je,kb:ke) = U(VELX_VAR:VELZ_VAR,ib:ie,jb:je,kb:ke)*hy_vref

!!$#ifdef FLASH_USM_MHD
!!$     U(MAGX_VAR:MAGZ_VAR,ib:ie,jb:je,kb:ke) = U(MAGX_VAR:MAGZ_VAR,ib:ie,jb:je,kb:ke)*hy_bref
!!$#if NFACE_VARS > 0
!!$#if NDIM > 1
!!$     if (NDIM >= 2) then
!!$        Bx(MAG_FACE_VAR,ib:ie,jb:je,kb:ke) = Bx(MAG_FACE_VAR,ib:ie,jb:je,kb:ke)*hy_bref
!!$        By(MAG_FACE_VAR,ib:ie,jb:je,kb:ke) = By(MAG_FACE_VAR,ib:ie,jb:je,kb:ke)*hy_bref
!!$        if (NDIM==3) Bz(MAG_FACE_VAR,ib:ie,jb:je,kb:ke) = Bz(MAG_FACE_VAR,ib:ie,jb:je,kb:ke)*hy_bref
!!$     endif
!!$#endif
!!$#endif
!!$#endif

  end select

!!$#ifdef FLASH_USM_MHD
!!$#if NFACE_VARS > 0
!!$#if NDIM > 1
!!$  call Grid_releaseBlkPtr(blockID,Bx,FACEX)
!!$  call Grid_releaseBlkPtr(blockID,By,FACEY)
!!$  if (NDIM==3) call Grid_releaseBlkPtr(blockID,Bz,FACEZ)
!!$#endif
!!$#endif
!!$#endif

End Subroutine hy_unitConvert

  

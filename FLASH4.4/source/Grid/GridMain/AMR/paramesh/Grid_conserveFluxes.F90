!!****if* source/Grid/GridMain/paramesh/Grid_conserveFluxes
!!
!! NAME
!!  Grid_conserveFluxes
!!
!! SYNOPSIS
!!
!!  call Grid_conserveFluxes(integer(IN) :: axis,
!!                           integer(IN) :: coarse_level)
!!  
!! DESCRIPTION 
!!  
!!  Flux conservation is necessary when 2 blocks of differing
!!  levels (meaning having different grid spacings) border 
!!  one another. 
!!  
!!  This routine can perform flux conservation on the finest
!!  blocks, the most typical usage for the Paramesh Grid or on
!!  blocks of a certain level.
!!  
!!  The routine overwrites the flux arrays maintained by the Grid
!!  
!! ARGUMENTS 
!!
!!
!!  axis - conserve fluxes in just one direction if 
!!         IAXIS, JAXIS, KAXIS, or in all directions if ALLDIR.
!!         These constants are defined in constants.h.
!!
!!  coarse_level - refinement level. Selects the level (coarse level) for
!!          which fluxes are updated.
!!          Can be UNSPEC_LEVEL for all levels (except, as an
!!          optimizing shortcut, the highest possible one).
!!
!!***
!!REORDER(5): flux_[xyz], gr_[xyz]flx
!!REORDER(4): fluxx,fluxy,fluxz
!!REORDER(5): gr_xflx_[yz]face, gr_yflx_[xz]face, gr_zflx_[xy]face
#include "Flash.h"
#include "constants.h"
subroutine Grid_conserveFluxes( axis, coarse_level)
  use paramesh_interfaces, ONLY : amr_flux_conserve
  use physicaldata, ONLY : flux_x, flux_y, flux_z, nfluxes
  use tree, ONLY : surr_blks, nodetype, lrefine_max
  use gr_specificData, ONLY : gr_xflx, gr_yflx, gr_zflx, gr_flxx, gr_flxy, gr_flxz
  use gr_specificData, ONLY : gr_iloFl, gr_jloFl, gr_kloFl
  use leaf_iterator, ONLY : leaf_iterator_t
  use block_metadata, ONLY : block_metadata_t
  use Grid_interface, ONLY : Grid_getLeafIterator, Grid_releaseLeafIterator

#ifndef FLASH_GRID_PARAMESH2
  use physicaldata, ONLY: no_permanent_guardcells
#endif
  use Grid_data, ONLY : gr_meshMe,gr_maxRefine
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >=2
  use gr_specificData, ONLY : gr_xflx_yface, gr_yflx_xface
#if NDIM == 3
  use gr_specificData, ONLY : gr_xflx_zface, gr_yflx_zface, gr_zflx_xface, gr_zflx_yface
#endif
#endif
#endif

  implicit none
  integer, intent(in) ::  axis, coarse_level
  integer :: gridDataStruct
  integer :: presVar, np
  integer,save,dimension(1),target :: presDefault = (/-1/)
  integer,pointer,dimension(:) :: presP
  integer :: sx,ex,sy,ey,sz,ez
  real,pointer, dimension(:,:,:,:) :: fluxx,fluxy,fluxz

  type(leaf_iterator_t)  :: itor
  type(block_metadata_t) :: blockDesc
  integer :: blockID
  logical :: xtrue,ytrue,ztrue
  integer, dimension(MDIM) :: datasize


  !! Dev:  Disabled immediate return for now, let the caller control this. - KW
  !if(coarse_level > 1) return

  if (axis == ALLDIR) then
     call amr_flux_conserve(gr_meshMe, 0, 0)
  else
     call amr_flux_conserve(gr_meshMe, 0, axis)     
  endif


#ifndef FLASH_GRID_PARAMESH2
  gridDataStruct = CENTER
#if NFACE_VARS > 0
  gridDataStruct = CENTER_FACES
#endif
  !! Dev -- AD I have no idea why the following code is there at all
  !! but keeping it around caused crash in parallel mode

  if (no_permanent_guardcells) then
     call gr_commSetup(gridDataStruct)
  else
     call gr_freeCommRecvBuffer
  end if
#endif
  dataSize(IAXIS)=NXB+2*NGUARD
  dataSize(JAXIS)=NYB+2*NGUARD*K2D
  dataSize(KAXIS)=NZB+2*NGUARD*K3D
  
  sx = NGUARD+1
  sy = NGUARD*K2D+1
  sz = NGUARD*K3D+1
  ex = dataSize(1)-NGUARD
  ey = dataSize(2)-NGUARD*K2D
  ez = dataSize(3)-NGUARD*K3D
  
  xtrue=.true.
  ytrue= (NDIM>1)
  ztrue= (NDIM>2)

  if(axis /= ALLDIR)then
     xtrue = (axis==IAXIS)
     ytrue = (axis==JAXIS)     
     ztrue = (axis==KAXIS)
  end if

  call Grid_getLeafIterator(itor, level=coarse_level)
  do while(itor%is_valid())
     call itor%blkMetaData(blockDesc)
     if ((coarse_level == UNSPEC_LEVEL) .AND. (blockDesc%level == lrefine_max)) then
        call itor%next()
        CYCLE !Skip blocks at highest level.
     end if
     blockID=blockDesc%id
     
     fluxx(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => gr_flxx(:,:,:,:,blockID)
     fluxy(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => gr_flxy(:,:,:,:,blockID)
     fluxz(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => gr_flxz(:,:,:,:,blockID)
!!$     if (present(pressureSlots)) then
!!$        presP => pressureSlots
!!$     else
!!$        presP => presDefault
!!$     end if
     
     if(xtrue) then
!!$#ifdef FLASH_HYDRO_UNSPLIT
!!$#if NDIM >= 2
!!$        fluxx(:,sx+1:ex,sy,sz:ez) = gr_xflx_yface(:,1,:,:,blockID)
!!$        fluxx(:,sx+1:ex,ey,sz:ez) = gr_xflx_yface(:,2,:,:,blockID)
!!$#if NDIM == 3
!!$        fluxx(:,sx+1:ex,sy:ey,sz) = gr_xflx_zface(:,:,1,:,blockID)
!!$        fluxx(:,sx+1:ex,sy:ey,ez) = gr_xflx_zface(:,:,2,:,blockID)
!!$#endif
!!$#endif
!!$#endif
        
        if (surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
            surr_blks(3,1,1+K2D,1+K3D,blockID) == PARENT_BLK) &
                   fluxx(:,sx,  sy:ey,sz:ez) = flux_x(:nfluxes,1,:,:,blockID)
        if (surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
            surr_blks(3,3,1+K2D,1+K3D,blockID) == PARENT_BLK) &
                   fluxx(:,ex+1,sy:ey,sz:ez) = flux_x(:nfluxes,2,:,:,blockID)
!!$        fluxx(:,sx+1,sy:ey,sz:ez) = gr_xflx(1,:,:,:,blockID)
!!$        fluxx(:,ex,  sy:ey,sz:ez) = gr_xflx(2,:,:,:,blockID)
        
!!$        do np = 1,size(presP,1)
!!$           presVar = presP(np)
!!$           if (presVar > 0) then
!!$              if (.NOT.(surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
!!$                   surr_blks(3,1,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
!!$                 where (areaLeft(sx,sy:ey,sz:ez).NE.0.0) &
!!$                      fluxx(:,sx,sy:ey,sz:ez,presVar) = fluxx(:,sx,sy:ey,sz:ez,presVar) / areaLeft(sx,sy:ey,sz:ez)
!!$              end if
!!$              if (.NOT.(surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
!!$                   surr_blks(3,3,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
!!$                 fluxx(:,ex+1,sy:ey,sz:ez,presVar) = fluxx(:,ex+1,sy:ey,sz:ez,presVar) / areaLeft(ex+1,sy:ey,sz:ez)
!!$              end if
!!$           end if
!!$        end do
     end if
     
     if(ytrue) then
!!$#ifdef FLASH_HYDRO_UNSPLIT
!!$#if NDIM >= 2
!!$        fluxy(:,sx,sy+1:ey,sz:ez) = gr_yflx_xface(1,:,:,:,blockID)
!!$        fluxy(:,ex,sy+1:ey,sz:ez) = gr_yflx_xface(2,:,:,:,blockID)
!!$#if NDIM == 3
!!$        fluxy(:,sx:ex,sy+1:ey,sz) = gr_yflx_zface(:,:,1,:,blockID)
!!$        fluxy(:,sx:ex,sy+1:ey,ez) = gr_yflx_zface(:,:,2,:,blockID)
!!$#endif
!!$#endif
!!$#endif
        if (surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
            surr_blks(3,2,1,1+K3D,blockID) == PARENT_BLK) &
                   fluxy(:,sx:ex,sy,  sz:ez) = flux_y(:nfluxes,:,1,:,blockID)
        if (surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
            surr_blks(3,2,3,1+K3D,blockID) == PARENT_BLK) &
                   fluxy(:,sx:ex,ey+1,sz:ez) = flux_y(:nfluxes,:,2,:,blockID)
!!$        fluxy(:,sx:ex,sy+1,sz:ez) = gr_yflx(:,1,:,:,blockID)
!!$        fluxy(:,sx:ex,ey,  sz:ez) = gr_yflx(:,2,:,:,blockID)
        
        
#if NDIM > 1
!!$        do np = 1,size(presP,1)
!!$           presVar = presP(np)
!!$           if (presVar > 0) then
!!$              if (.NOT.(surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
!!$                   surr_blks(3,2,1,1+K3D,blockID) == nodetype(blockID))) then
!!$                 where (areaLeft(sx:ex,sy,sz:ez).NE.0.0) &
!!$                      fluxy(sx:ex,sy,sz:ez,presVar) = fluxy(sx:ex,sy,sz:ez,presVar) / areaLeft(sx:ex,sy,sz:ez)
!!$              end if
!!$              if (.NOT.(surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
!!$                   surr_blks(3,2,3,1+K3D,blockID) == nodetype(blockID))) then
!!$                 where (areaLeft(sx:ex,ey+1,sz:ez).NE.0.0) &
!!$                      fluxy(sx:ex,ey+1,sz:ez,presVar) = fluxy(sx:ex,ey+1,sz:ez,presVar) / areaLeft(sx:ex,ey+1,sz:ez)
!!$              end if
!!$           end if
!!$        end do
#endif
     end if
     
     if(ztrue) then
!!$#ifdef FLASH_HYDRO_UNSPLIT
!!$#if NDIM == 3
!!$        fluxz(:,sx,sy:ey,sz+1:ez) = gr_zflx_xface(1,:,:,:,blockID)
!!$        fluxz(:,ex,sy:ey,sz+1:ez) = gr_zflx_xface(2,:,:,:,blockID)
!!$        fluxz(:,sx:ex,sy,sz+1:ez) = gr_zflx_yface(:,1,:,:,blockID)
!!$        fluxz(:,sx:ex,ey,sz+1:ez) = gr_zflx_yface(:,2,:,:,blockID)
!!$#endif
!!$#endif
        if (surr_blks(1,2,2,1,blockID) > 0 .AND. &
            surr_blks(3,2,2,1,blockID) == PARENT_BLK) &
                   fluxz(:,sx:ex,sy:ey,sz  ) = flux_z(:nfluxes,:,:,1,blockID)
        if (surr_blks(1,2,2,3,blockID) > 0 .AND. &
            surr_blks(3,2,2,3,blockID) == PARENT_BLK) &
                   fluxz(:,sx:ex,sy:ey,ez+1) = flux_z(:nfluxes,:,:,2,blockID)
!!$        fluxz(:,sx:ex,sy:ey,sz+1) = gr_zflx(:,:,1,:,blockID)
!!$        fluxz(:,sx:ex,sy:ey,ez  ) = gr_zflx(:,:,2,:,blockID)
#if NDIM > 2
!!$        do np = 1,size(presP,1)
!!$           presVar = presP(np)
!!$           if (presVar > 0) then
!!$        where (areaLeft(sx:ex,sy:ey,sz).NE.0.0) &  ! should not happen in any supported geometry
!!$              if (.NOT.(surr_blks(1,2,2,1,blockID) > 0 .AND. &
!!$                   surr_blks(3,2,2,1,blockID) == nodetype(blockID))) then
!!$                 fluxz(:,sx:ex,sy:ey,sz,presVar) = fluxz(:,sx:ex,sy:ey,sz,presVar) / areaLeft(sx:ex,sy:ey,sz)
!!$              end if
!!$              if (.NOT.(surr_blks(1,2,2,3,blockID) > 0 .AND. &
!!$                   surr_blks(3,2,2,3,blockID) == nodetype(blockID))) then
!!$                 fluxz(:,sx:ex,sy:ey,ez+1,presVar) = fluxz(:,sx:ex,sy:ey,ez+1,presVar) / areaLeft(sx:ex,sy:ey,ez+1)
!!$              end if
!!$           end if
!!$     end do
#endif
     end if
     call itor%next()
  end do
  call Grid_releaseLeafIterator(itor)
  return
end subroutine Grid_conserveFluxes
   



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
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator, &
                             Grid_getCellFaceAreas

#ifndef FLASH_GRID_PARAMESH2
  use physicaldata, ONLY: no_permanent_guardcells
#endif
  use Grid_data, ONLY : gr_meshMe,gr_maxRefine, &
                        gr_geometry
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

  type(Grid_iterator_t)  :: itor
  type(Grid_tile_t) :: tileDesc
  integer :: blockID
  logical :: xtrue,ytrue,ztrue
  logical :: divideFluxx,divideFluxy,divideFluxz !whether to divide by area
  integer, dimension(MDIM) :: datasize
  integer, dimension(MDIM) :: offs
  integer                  :: offx, offy, offz
  integer                  :: blkLev
  logical                  :: loCase4, hiCase4
  real,allocatable,target :: faceAreas(:,:,:)
  real,pointer            :: areaLeft(:,:,:)


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

  select case (gr_geometry)
  case(CARTESIAN)
     divideFluxx = .false.
     divideFluxy = .false.
     divideFluxz = .false.
  case default
     divideFluxx = .TRUE.
     divideFluxy = .TRUE.
     divideFluxz = .TRUE.
  end select


  call Grid_getTileIterator(itor, LEAF, level=coarse_level, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     if ((coarse_level == UNSPEC_LEVEL) .AND. (tileDesc%level == lrefine_max)) then
        call itor%next()
        CYCLE !Skip blocks at highest level.
     end if
     blockID=tileDesc%id
     blkLev =tileDesc%level
     
     fluxx(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => gr_flxx(:,:,:,:,blockID)
     fluxy(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => gr_flxy(:,:,:,:,blockID)
     fluxz(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => gr_flxz(:,:,:,:,blockID)
!!$     if (present(pressureSlots)) then
!!$        presP => pressureSlots
!!$     else
!!$        presP => presDefault
!!$     end if
     
     offs(:) = tileDesc%blkLimitsGC(LOW,1:MDIM) - 1
     offx = offs(IAXIS); offy = offs(JAXIS); offz = offs(KAXIS)

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
        
        ! With PARAMESH, there are four cases for the (same-level) neighbor in
        ! a certain direction of a *LEAF* block at refinement level blkLev:
        !
        ! Case | surr_blks(1,...) | surr_blks(3,...)  || Description: what's there?
        ! ======================= | ================= || ============================
        ! i.   !      <= -20      |  [ LEAF ? ]       || domain boundary (non-PERIODIC)
        ! ii.  |        -1        |  [  ignored  ]    || coarser block
        ! iii. | neighBlkID > 0   |  1  (LEAF)        || same refinement leaf block
        ! iv.  | neighBlkID > 0   |  2  (PARENT)      || finer blocks

        ! * We copy from flux_{x,y,z} in case iv.
        ! * We divide by face areas if
        !   - other conditions are satisfied (geometry and direction), and
        !   - face area is nonzero, and
        !   - case is i, ii, or iv.
        !     (Note - This could probably be reduced to case iv only!)

        loCase4 = (surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
            surr_blks(3,1,1+K2D,1+K3D,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
            surr_blks(3,3,1+K2D,1+K3D,blockID) == PARENT_BLK)

        if (loCase4) &
                   fluxx(:,sx,  sy:ey,sz:ez) = flux_x(:nfluxes,1,:,:,blockID)
        if (hiCase4) &
                   fluxx(:,ex+1,sy:ey,sz:ez) = flux_x(:nfluxes,2,:,:,blockID)
!!$        fluxx(:,sx+1,sy:ey,sz:ez) = gr_xflx(1,:,:,:,blockID)
!!$        fluxx(:,ex,  sy:ey,sz:ez) = gr_xflx(2,:,:,:,blockID)

        if (divideFluxx) then
           allocate(faceAreas(offx+sx:offx+ex+1,offy+sy:offy+ey,offz+sz:offz+ez))
           areaLeft(sx:,sy:,sz:)  => faceAreas
           
           call Grid_getCellFaceAreas(IAXIS, blkLev, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)
           do presVar = 1,nfluxes
              if (loCase4) then
                 where (areaLeft(sx,sy:ey,sz:ez).NE.0.0) &
                      fluxx(presVar,sx,sy:ey,sz:ez) = fluxx(presVar,sx,sy:ey,sz:ez) / areaLeft(sx,sy:ey,sz:ez)
              end if
              if (hiCase4) then
                 fluxx(presVar,ex+1,sy:ey,sz:ez) = fluxx(presVar,ex+1,sy:ey,sz:ez) / areaLeft(ex+1,sy:ey,sz:ez)
              end if
           end do
           deallocate(faceAreas)
        end if
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
        loCase4 = (surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
            surr_blks(3,2,1,1+K3D,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
            surr_blks(3,2,3,1+K3D,blockID) == PARENT_BLK)

        if (loCase4) &
                   fluxy(:,sx:ex,sy,  sz:ez) = flux_y(:nfluxes,:,1,:,blockID)
        if (hiCase4) &
                   fluxy(:,sx:ex,ey+1,sz:ez) = flux_y(:nfluxes,:,2,:,blockID)
!!$        fluxy(:,sx:ex,sy+1,sz:ez) = gr_yflx(:,1,:,:,blockID)
!!$        fluxy(:,sx:ex,ey,  sz:ez) = gr_yflx(:,2,:,:,blockID)
        
        
#if NDIM > 1
        if (divideFluxy) then
           allocate(faceAreas(offx+sx:offx+ex,offy+sy:offy+ey+1,offz+sz:offz+ez))
           areaLeft(sx:,sy:,sz:)  => faceAreas
           
           call Grid_getCellFaceAreas(JAXIS, blkLev, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)
           do presVar = 1,nfluxes
              if (loCase4) then
                 where (areaLeft(sx:ex,sy,sz:ez).NE.0.0) &
                      fluxy(presVar,sx:ex,sy,sz:ez) = fluxy(presVar,sx:ex,sy,sz:ez) / areaLeft(sx:ex,sy,sz:ez)
              end if
              if (hiCase4) then
                 where (areaLeft(sx:ex,ey+1,sz:ez).NE.0.0) &
                      fluxy(presVar,sx:ex,ey+1,sz:ez) = fluxy(presVar,sx:ex,ey+1,sz:ez) / areaLeft(sx:ex,ey+1,sz:ez)
              end if
           end do
           deallocate(faceAreas)
        end if
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
        loCase4 = (surr_blks(1,2,2,1,blockID) > 0 .AND. &
            surr_blks(3,2,2,1,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,2,2,3,blockID) > 0 .AND. &
            surr_blks(3,2,2,3,blockID) == PARENT_BLK)

        if (loCase4) &
                   fluxz(:,sx:ex,sy:ey,sz  ) = flux_z(:nfluxes,:,:,1,blockID)
        if (hiCase4) &
                   fluxz(:,sx:ex,sy:ey,ez+1) = flux_z(:nfluxes,:,:,2,blockID)
!!$        fluxz(:,sx:ex,sy:ey,sz+1) = gr_zflx(:,:,1,:,blockID)
!!$        fluxz(:,sx:ex,sy:ey,ez  ) = gr_zflx(:,:,2,:,blockID)
#if NDIM > 2
        if (divideFluxz) then
           allocate(faceAreas(offx+sx:offx+ex,offy+sy:offy+ey,offz+sz:offz+ez+1))
           areaLeft(sx:,sy:,sz:)  => faceAreas
           
           call Grid_getCellFaceAreas(KAXIS, blkLev, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)
           do presVar = 1,nfluxes
!        where (areaLeft(sx:ex,sy:ey,sz).NE.0.0) &  ! should not happen in any supported geometry
              if (loCase4) then
                 fluxz(presVar,sx:ex,sy:ey,sz) = fluxz(presVar,sx:ex,sy:ey,sz) / areaLeft(sx:ex,sy:ey,sz)
              end if
              if (hiCase4) then
                 fluxz(presVar,sx:ex,sy:ey,ez+1) = fluxz(presVar,sx:ex,sy:ey,ez+1) / areaLeft(sx:ex,sy:ey,ez+1)
              end if
           end do
           deallocate(faceAreas)
        end if
#endif
     end if
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
  return
end subroutine Grid_conserveFluxes
   



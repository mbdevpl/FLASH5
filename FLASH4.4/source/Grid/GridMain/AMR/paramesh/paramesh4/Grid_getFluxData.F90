!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_getFluxData
!!
!! NAME
!!  Grid_getFluxData
!!
!! SYNOPSIS
!!
!!
!!  call Grid_getFluxData(integer(IN) :: blockID,
!!                   integer(IN) :: axis,
!!                   real(INOUT) :: fluxes(NFLUXES,dataSize(1),dataSize(2),dataSize(3)),
!!                   integer(IN) :: dataSize(3),
!!          OPTIONAL,integer(IN) :: pressureSlots,
!!          OPTIONAL,real(IN)    :: areaLeft(:,:,:))
!!
!! DESCRIPTION 
!!
!!  Get the fluxes in a direction specified by axis for boundary cells
!!  for block blockID. This routine needs to be used when using adaptive mesh
!!  since fluxes calculated by the two blocks that are at a fine/coarse
!!  boundary have different accuracy.
!!  
!!  This should be called after Grid_conserveFluxes, which gets the 
!!  fluxes updated for cells at fine/coarse boundaries and makes them consistent.
!!
!! ARGUMENTS
!!
!!  blockID : The local blockid
!!
!!
!!  axis : integer value specifying on which cell faces to get fluxes. 
!!         The options are IAXIS, JAXIS, or KAXIS defined in constants.h
!!
!!
!!  fluxes :  real array with space for fluxes, through one axis, 
!!            for all cells of a block and for all flux variables.
!!            fluxes(VAR, i, j, k) is VAR's flux through 
!!            the left cell face for cell i, j, k.  The 
!!            fluxes of the boundary cells of coarse blocks
!!            will have been appropriately changed, if the flux matching
!!            between fine and coarse boundaries (normally by a call to
!!            Grid_conserveFluxes), has occured since the last call
!!            of Grid_putFluxData.
!!
!!
!!  dataSize : integer array specifying the dimensions for fluxes
!!
!!             dataSize (1) holds the number of cells returned in the i direction
!!
!!             dataSize (2) holds the number of cells returned in the j direction
!!                          if 1 d problem, set datasize(2) = 1
!!
!!             dataSize (3) holds the number of cells returned in the k direction
!!                          if 1 or 2 d problem, set datasize(3) = 1
!!
!!             fluxes should contain space for fluxes of all cells in the block, 
!!             including guardcells.
!!
!!  pressureSlots: an optional array of integer values indicating fluxes that are to be
!!                 handled like pressures, rather than genuine fluxes.  If present, each
!!                 element greater than zero indicates one flux variable in the fluxes
!!                 array that may need special handling because it really does not scale
!!                 like a flux. Normally this would be pressure, or (e.g., with multiTemp
!!                 physics) several partial pressures given at cell interfaces. The
!!                 pressureSlots array could also list other flux variables that the
!!                 caller keeps in flux density form.
!!
!!  areaLeft :     areas of left and right faces, only used if special scaling is
!!                 requested with the pressureSlots argument.
!!
!! NOTES 
!!
!!   Any code calling this subroutine needs to know the explicit interface,
!!   since this interface contains optional dummy arguments and assumed-shape
!!   dummy arrays. Calling FORTRAN units should therefore contain a line like
!!       use Grid_interface, ONLY: Grid_getFluxData
!!
!!   This implementation is specific to Paramesh 4.
!!
!! SEE ALSO
!!
!!   Grid_putFluxData
!!   Grid_conserveFluxes
!!   hy_sweep
!!***

!!REORDER(5): flux_[xyz], gr_[xyz]flx
!!REORDER(4): fluxx,fluxy,fluxz
!!REORDER(5): gr_xflx_[yz]face, gr_yflx_[xz]face, gr_zflx_[xy]face
#include "Flash.h"
subroutine Grid_getFluxData(block, fluxx, fluxy, fluxz, dataSize,axis,  pressureSlots, areaLeft)

  use physicaldata, ONLY : flux_x, flux_y, flux_z, nfluxes
  use tree, ONLY : surr_blks, nodetype
  use gr_specificData, ONLY : gr_xflx, gr_yflx, gr_zflx
  use block_metadata, ONLY : block_metadata_t

#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >=2
  use gr_specificData, ONLY : gr_xflx_yface, gr_yflx_xface
#if NDIM == 3
  use gr_specificData, ONLY : gr_xflx_zface, gr_yflx_zface, gr_zflx_xface, gr_zflx_yface
#endif
#endif
#endif
  implicit none

!#include "Flash.h"
#include "constants.h"

  type(block_metadata_t), intent(IN) :: block
  integer, intent(IN),optional :: axis
  integer, intent(IN), dimension(3) :: dataSize
  real, intent(INOUT), dimension(NFLUXES,dataSize(1),dataSize(2),dataSize(3)) :: fluxx,fluxy,fluxz
  integer, intent(IN), OPTIONAL,target :: pressureSlots(:)
  real, intent(IN), OPTIONAL :: areaLeft(:,:,:)
  integer :: blockID
  logical :: xtrue, ytrue, ztrue

  integer :: presVar, np
  integer,save,dimension(1),target :: presDefault = (/-1/)
  integer,pointer,dimension(:) :: presP
  integer :: sx,ex,sy,ey,sz,ez

  if (NFLUXES > 0) return
  if (present(pressureSlots)) then
     presP => pressureSlots
  else
     presP => presDefault
  end if

  sx = NGUARD+1
  sy = NGUARD*K2D+1
  sz = NGUARD*K3D+1
  ex = dataSize(1)-NGUARD
  ey = dataSize(2)-NGUARD*K2D
  ez = dataSize(3)-NGUARD*K3D

  blockID=block%id
  xtrue=.true.
  ytrue= (NDIM>1)
  ztrue= (NDIM>2)

  if(present(axis))then
     xtrue = (axis==IAXIS)
     ytrue = (axis==JAXIS)     
     ztrue = (axis==KAXIS)
  end if

  if(xtrue) then
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >= 2
     fluxx(:,sx+1:ex,sy,sz:ez) = gr_xflx_yface(:,:,1,:,blockID)
     fluxx(:,sx+1:ex,ey,sz:ez) = gr_xflx_yface(:,:,2,:,blockID)
#if NDIM == 3
     fluxx(:,sx+1:ex,sy:ey,sz) = gr_xflx_zface(:,:,:,1,blockID)
     fluxx(:,sx+1:ex,sy:ey,ez) = gr_xflx_zface(:,:,:,2,blockID)
#endif
#endif
#endif
     
     fluxx(:,sx,  sy:ey,sz:ez) = flux_x(:nfluxes,1,:,:,blockID) 
     fluxx(:,ex+1,sy:ey,sz:ez) = flux_x(:nfluxes,2,:,:,blockID)
     fluxx(:,sx+1,sy:ey,sz:ez) = gr_xflx(:,1,:,:,blockID) 
     fluxx(:,ex,  sy:ey,sz:ez) = gr_xflx(:,2,:,:,blockID)
     
     do np = 1,size(presP,1)
        presVar = presP(np)
        if (presVar > 0) then
           if (.NOT.(surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
                surr_blks(3,1,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
              where (areaLeft(sx,sy:ey,sz:ez).NE.0.0) &
                   fluxx(presVar,sx,sy:ey,sz:ez) = fluxx(presVar,sx,sy:ey,sz:ez) / areaLeft(sx,sy:ey,sz:ez)
           end if
           if (.NOT.(surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
                surr_blks(3,3,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
              fluxx(presVar,ex+1,sy:ey,sz:ez) = fluxx(presVar,ex+1,sy:ey,sz:ez) / areaLeft(ex+1,sy:ey,sz:ez)
           end if
        end if
     end do
  end if

  if(ytrue) then
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >= 2
     fluxy(:,sx,sy+1:ey,sz:ez) = gr_yflx_xface(:,1,:,:,blockID)
     fluxy(:,ex,sy+1:ey,sz:ez) = gr_yflx_xface(:,2,:,:,blockID)
#if NDIM == 3
     fluxy(:,sx:ex,sy+1:ey,sz) = gr_yflx_zface(:,:,:,1,blockID)
     fluxy(:,sx:ex,sy+1:ey,ez) = gr_yflx_zface(:,:,:,2,blockID)
#endif
#endif
#endif
     fluxy(:,sx:ex,sy,  sz:ez) = flux_y(:nfluxes,:,1,:,blockID) 
     fluxy(:,sx:ex,ey+1,sz:ez) = flux_y(:nfluxes,:,2,:,blockID) 
     fluxy(:,sx:ex,sy+1,sz:ez) = gr_yflx(:,:,1,:,blockID) 
     fluxy(:,sx:ex,ey,  sz:ez) = gr_yflx(:,:,2,:,blockID)


#if NDIM > 1
     do np = 1,size(presP,1)
        presVar = presP(np)
        if (presVar > 0) then
           if (.NOT.(surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,2,1,1+K3D,blockID) == nodetype(blockID))) then
              where (areaLeft(sx:ex,sy,sz:ez).NE.0.0) &
                fluxy(presVar,sx:ex,sy,sz:ez) = fluxy(presVar,sx:ex,sy,sz:ez) / areaLeft(sx:ex,sy,sz:ez)
           end if
           if (.NOT.(surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,2,3,1+K3D,blockID) == nodetype(blockID))) then
              where (areaLeft(sx:ex,ey+1,sz:ez).NE.0.0) &
                   fluxy(presVar,sx:ex,ey+1,sz:ez) = fluxy(presVar,sx:ex,ey+1,sz:ez) / areaLeft(sx:ex,ey+1,sz:ez)
           end if
        end if
     end do
#endif
  end if

  if(ztrue) then
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM == 3
     fluxz(:,sx,sy:ey,sz+1:ez) = gr_zflx_xface(:,1,:,:,blockID)
     fluxz(:,ex,sy:ey,sz+1:ez) = gr_zflx_xface(:,2,:,:,blockID)
     fluxz(:,sx:ex,sy,sz+1:ez) = gr_zflx_yface(:,:,1,:,blockID)
     fluxz(:,sx:ex,ey,sz+1:ez) = gr_zflx_yface(:,:,2,:,blockID)
#endif
#endif
     fluxz(:,sx:ex,sy:ey,sz  ) = flux_z(:nfluxes,:,:,1,blockID) 
     fluxz(:,sx:ex,sy:ey,ez+1) = flux_z(:nfluxes,:,:,2,blockID) 
     fluxz(:,sx:ex,sy:ey,sz+1) = gr_zflx(:,:,:,1,blockID) 
     fluxz(:,sx:ex,sy:ey,ez  ) = gr_zflx(:,:,:,2,blockID)
#if NDIM > 2
     do np = 1,size(presP,1)
        presVar = presP(np)
        if (presVar > 0) then
!!$        where (areaLeft(sx:ex,sy:ey,sz).NE.0.0) &  ! should not happen in any supported geometry
           if (.NOT.(surr_blks(1,2,2,1,blockID) > 0 .AND. &
                surr_blks(3,2,2,1,blockID) == nodetype(blockID))) then
              fluxz(presVar,sx:ex,sy:ey,sz) = fluxz(presVar,sx:ex,sy:ey,sz) / areaLeft(sx:ex,sy:ey,sz)
           end if
           if (.NOT.(surr_blks(1,2,2,3,blockID) > 0 .AND. &
                surr_blks(3,2,2,3,blockID) == nodetype(blockID))) then
              fluxz(presVar,sx:ex,sy:ey,ez+1) = fluxz(presVar,sx:ex,sy:ey,ez+1) / areaLeft(sx:ex,sy:ey,ez+1)
           end if
        end if
     end do
#endif
  end if

  return
end subroutine Grid_getFluxData



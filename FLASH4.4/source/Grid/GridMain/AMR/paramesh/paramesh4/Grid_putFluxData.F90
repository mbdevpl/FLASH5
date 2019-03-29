!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_putFluxData
!!
!! NAME
!!  Grid_putFluxData
!!
!! SYNOPSIS
!!
!!
!!  call Grid_putFluxData(integer(IN) :: blockID,
!!                   integer(IN) :: axis,
!!                   real(IN)    :: fluxes(NFLUXES,dataSize(1),dataSize(2),dataSize(3)),
!!                   integer(IN) :: dataSize(3),
!!          OPTIONAL,integer(IN) :: pressureSlots(:),
!!          OPTIONAL,real(IN)    :: areaLeft(:,:,:))
!!
!! DESCRIPTION 
!!
!!
!!  Put the fluxes in a direction specified by axis for boundary cells
!!  for block blockID. This routine needs to be used with adaptive mesh
!!  since fluxes calculated by the two blocks that are at fine/coarse boundary have 
!!  different accuracy. The fluxes calculated by individual blocks are reported to 
!!  the Grid through this call. Once that is done, a call to Grid_conserveFluxes 
!!  applies the flux conservation algorithm to make it consistent across the fine/coarse 
!!  boundaries.
!!
!! ARGUMENTS
!!
!!  blockID : The local blockid
!!
!!
!!  axis : integer value specifying on which cell faces to put fluxes. 
!!         The options are IAXIS, JAXIS, or KAXIS defined in constants.h
!!
!!
!!  fluxes :  real array with space for fluxes, through one axis, 
!!            for all cells of a block and for all flux variables.
!!            fluxes(VAR, i, j, k) is VAR's flux through 
!!            the left cell face for cell i, j, k.
!!
!!
!!  dataSize : integer array specifying the dimensions for the array, fluxes
!!
!!             dataSize (1) holds the number of cells provided in the i direction
!!
!!             dataSize (2) holds the number of cells provided in the j direction
!!                          if 1 d problem, set datasize(2) = 1
!!
!!             dataSize (3) holds the number of cells provided in the k direction
!!                          if 1 or 2 d problem, set datasize(3) = 1
!!
!!             fluxes should contain space for fluxes of all cells in the block, 
!!             including guardcells, and the  fluxes must be correct for 
!!             the interior cells of the block, as this interface does not know which 
!!             cell fluxes the Grid will need to store.
!!
!!  pressureSlots: If present and greater than zero, this indicates one or more flux
!!                 variables in the fluxes array that may need special handling because
!!                 they really scale like flux densities when other flux variables scale
!!                 like fluxes.  For the split PPM Hydro implementation, for example,
!!                 this is normally used for the pressure "flux" variable, but
!!                 it can be applied to other flux variables that the caller keeps in
!!                 flux density form.
!!                 The special handling consists in multiplying the flux variables with
!!                 the appropriate face areas, which are taken from the areaLeft array
!!                 argument, before storing them in the arrays on which the Grid unit
!!                 acts.
!!                 Special handling should only be requested if the Grid units handles
!!                 flux variable "as fluxes" (or else if it does nto matter anyway,
!!                 as is the case in 1D).
!!                 The pressureSlots argument in the corresponding Grid_getFluxData
!!                 call should generally match the one in the Grid_putFluxData call.
!!
!!  areaLeft :     areas of left and right faces, only used if special scaling is
!!                 requested with the pressureSlot argument.
!!                 The areaLeft argument in the corresponding Grid_getFluxData
!!                 call should generally match the one in the Grid_putFluxData call.
!!
!! NOTES 
!!
!!   Any code calling this subroutine needs to know the explicit interface,
!!   since this interface contains optional dummy arguments and assumed-shape
!!   dummy arrays. Calling FORTRAN units should therefore contain a line like
!!       use Grid_interface, ONLY: Grid_putFluxData
!!
!!   This implementation is specific to Paramesh 4.
!!
!! SEE ALSO
!!
!!   Grid_getFluxData
!!   Grid_conserveFluxes
!!   hy_sweep
!!***

!!REORDER(5): flux_[xyz],gr_[xyz]flx
!!REORDER(4): fluxx,fluxy,fluxz
!!REORDER(5): gr_xflx_[yz]face, gr_yflx_[xz]face, gr_zflx_[xy]face
#include "Flash.h"
subroutine Grid_putFluxData(level, axis, pressureSlots, areaLeftIGNORE)

  use physicaldata, ONLY : flux_x, flux_y, flux_z, nfluxes
  use tree, ONLY : surr_blks, nodetype
  use gr_specificData, ONLY : gr_xflx, gr_yflx, gr_zflx, gr_flxx, gr_flxy, gr_flxz
  use gr_specificData, ONLY : gr_iloFl, gr_jloFl, gr_kloFl
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator, &
                             Grid_getCellFaceAreas
  use Grid_data,      ONLY : gr_geometry
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >=2
  use gr_specificData, ONLY : gr_xflx_yface, gr_yflx_xface
#if NDIM == 3
  use gr_specificData, ONLY : gr_xflx_zface, gr_yflx_zface, gr_zflx_xface, gr_zflx_yface
#endif
#endif
#endif

  implicit none

#include "constants.h"

  integer, intent(IN) :: level
  integer, intent(IN),optional :: axis
  integer, intent(IN), OPTIONAL,target :: pressureSlots(:)
  real, intent(IN), OPTIONAL :: areaLeftIGNORE(:,:,:)
  real,pointer, dimension(:,:,:,:) :: fluxx,fluxy,fluxz

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc

  integer :: blockID

  integer :: presVar, np
  integer,save,dimension(1),target :: presDefault = (/-1/)
  integer,pointer,dimension(:) :: presP
  integer :: sx,ex,sy,ey,sz,ez
  logical :: xtrue,ytrue,ztrue
  logical :: multFluxx,multFluxy,multFluxz !whether to multiply by area
  integer, dimension(MDIM) :: datasize
  integer, dimension(MDIM) :: offs
  integer                  :: offx, offy, offz
  integer                  :: blkLev
  real,allocatable,target :: faceAreas(:,:,:)
  real,pointer            :: areaLeft(:,:,:)

  !! Dev - AD for AMReX, it should be if((level==1).or.(NFLUXES<1))return
  if(NFLUXES<1)return
  
  dataSize(IAXIS)=NXB+2*NGUARD
  dataSize(JAXIS)=NYB+2*NGUARD*K2D
  dataSize(KAXIS)=NZB+2*NGUARD*K3D

!! This will need to be changed so that they can be computed here  
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

  xtrue=.true.
  ytrue= (NDIM>1)
  ztrue= (NDIM>2)
  
  if(present(axis))then
     xtrue = (axis==IAXIS)
     ytrue = (axis==JAXIS)     
     ztrue = (axis==KAXIS)
  end if
  
  select case (gr_geometry)
  case(CARTESIAN)
     multFluxx = .false.
     multFluxy = .false.
     multFluxz = .false.
  case default
     multFluxx = .TRUE.
     multFluxy = .TRUE.
     multFluxz = .TRUE.
  end select

  call Grid_getTileIterator(itor, LEAF, level=level, tiling=.FALSE.)

  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     if ((level == UNSPEC_LEVEL) .AND. (tileDesc%level == 1)) then
        call itor%next()
        CYCLE !Skip blocks at level 1.
     end if
     blockID=tileDesc%id
     blkLev =tileDesc%level
     fluxx(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => gr_flxx(:,:,:,:,blockID)
     fluxy(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => gr_flxy(:,:,:,:,blockID)
     fluxz(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => gr_flxz(:,:,:,:,blockID)

     offs(:) = tileDesc%blkLimitsGC(LOW,1:MDIM) - 1
     offx = offs(IAXIS); offy = offs(JAXIS); offz = offs(KAXIS)

     if(xtrue) then
        flux_x(:nfluxes,1,:,:,blockID) = fluxx(:,sx,sy:ey,sz:ez) 
        flux_x(:nfluxes,2,:,:,blockID) = fluxx(:,ex+1,sy:ey,sz:ez) 
!!$        gr_xflx(:,1,:,:,blockID) = fluxx(:,sx+1,sy:ey,sz:ez)
!!$        gr_xflx(:,2,:,:,blockID) = fluxx(:,ex,sy:ey,sz:ez)
!!$#ifdef FLASH_HYDRO_UNSPLIT
!!$        !! Store transverse components for the faces in global scratch arrays.
!!$#if NDIM >=2
!!$        gr_xflx_yface(:,:,1,:,blockID) = fluxx(:,sx+1:ex,sy,sz:ez)
!!$        gr_xflx_yface(:,:,2,:,blockID) = fluxx(:,sx+1:ex,ey,sz:ez)
!!$#if NDIM>2
!!$        gr_xflx_zface(:,:,:,1,blockID) = fluxx(:,sx+1:ex,sy:ey,sz)
!!$        gr_xflx_zface(:,:,:,2,blockID) = fluxx(:,sx+1:ex,sy:ey,ez)
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

        ! * We copy to flux_{x,y,z} in cases i,ii,iii,iv (e.g., above)
        ! * We multiply with face areas if
        !   - other conditions are satisfied (geometry and direction), and
        !   - face area is nonzero, and
        !   - case is i, ii, or iv.

        if (multFluxx) then

           allocate(faceAreas(offx+sx:offx+ex+1,offy+sy:offy+ey,offz+sz:offz+ez))
           areaLeft(sx:,sy:,sz:)  => faceAreas
           
           call Grid_getCellFaceAreas(IAXIS, blkLev, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)
           if (blkLev == 4 .AND. blockID == 28 .OR. &
               blkLev == 3 .AND. blockID == 70) then
           print*,'SHAPE (fluxx) is',SHAPE (fluxx),' for blockID,blkLev=',blockID, blkLev
           print*,'LBOUND(fluxx) is',LBOUND(fluxx)
           print*,'UBOUND(fluxx) is',UBOUND(fluxx)
!!$           print*,'SHAPE (fluxy) is',SHAPE (fluxy),' for blockID,blkLev=',blockID, blkLev
!!$           print*,'LBOUND(fluxy) is',LBOUND(fluxy)
!!$           print*,'UBOUND(fluxy) is',UBOUND(fluxy)
           print*,'SHAPE (faceA) is',SHAPE (faceAreas),' for blockID,blkLev=',blockID, blkLev
           print*,'LBOUND(faceA) is',LBOUND(faceAreas)
           print*,'UBOUND(faceA) is',UBOUND(faceAreas)
           print*,'SHAPE (areaL) is',SHAPE (areaLeft),' for blockID,blkLev=',blockID, blkLev
           print*,'LBOUND(areaL) is',LBOUND(areaLeft)
           print*,'UBOUND(areaL) is',UBOUND(areaLeft)
           print*,'areaL sx  :',areaLeft(sx,sy,sz)  ,fluxx(1,sx,sy,sz),fluxx(2,sx,sy,sz)
           print*,'areaL ex+1:',areaLeft(ex+1,sy,sz),fluxx(1,ex+1,sy,sz),fluxx(2,ex+1,sy,sz)
           end if
           
           do presVar = 1,nfluxes
              if (.NOT.(surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,1,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
!!$                 print*,'flux_x L bef:',blockID,blkLev,presVar,flux_x(presVar,1,:,:,blockID)
                 where (areaLeft(sx,sy:ey,sz:ez).NE.0.0)
                    flux_x(presVar,1,:,:,blockID) = flux_x(presVar,1,:,:,blockID) * areaLeft(sx,sy:ey,sz:ez)
                 end where
!!$                 print*,'flux_x L aft:',blockID,blkLev,presVar,flux_x(presVar,1,:,:,blockID)
              end if
              if (.NOT.(surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,3,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
!!$                 print*,'flux_x R bef:',blockID,blkLev,presVar,flux_x(presVar,2,:,:,blockID)
                 flux_x(presVar,2,:,:,blockID) = flux_x(presVar,2,:,:,blockID) * areaLeft(ex+1,sy:ey,sz:ez)
!!$                 print*,'flux_x R aft:',blockID,blkLev,presVar,flux_x(presVar,2,:,:,blockID)
              end if
           end do
           if (blkLev == 4 .AND. blockID == 28 .OR. &
               blkLev == 3 .AND. blockID == 70) then
!!$              print*,'AREAL SX  :',areaLeft(sx,sy,sz)  ,fluxx(1,sx,sy,sz),fluxx(2,sx,sy,sz)
!!$              print*,'AREAL EX+1:',areaLeft(ex+1,sy,sz),fluxx(1,ex+1,sy,sz),fluxx(2,ex+1,sy,sz)
              print*,'AREAL SX  :',areaLeft(sx,sy,sz)  ,flux_x(1,1,sy,sz,blockID),flux_x(2,1,sy,sz,blockID)
              print*,'AREAL EX+1:',areaLeft(ex+1,sy,sz),flux_x(1,2,sy,sz,blockID),flux_x(2,2,sy,sz,blockID)
           end if
           deallocate(faceAreas)
        end if
        
     end if
     
#if NDIM>1
     if(ytrue) then
        flux_y(:nfluxes,:,1,:,blockID)  = fluxy(:,sx:ex,sy,sz:ez)
        flux_y(:nfluxes,:,2,:,blockID)  = fluxy(:,sx:ex,ey+1,sz:ez) 
!!$        gr_yflx(:,:,1,:,blockID) =  fluxy(:,sx:ex,sy+1,sz:ez)
!!$        gr_yflx(:,:,2,:,blockID) =  fluxy(:,sx:ex,ey,sz:ez)
!!$#ifdef FLASH_HYDRO_UNSPLIT
!!$        !! Store transverse components for the faces in global scratch arrays.
!!$        gr_yflx_xface(:,1,:,:,blockID) = fluxy(:,sx,sy+1:ey,sz:ez)
!!$        gr_yflx_xface(:,2,:,:,blockID) = fluxy(:,ex,sy+1:ey,sz:ez)
!!$#if NDIM>2
!!$        gr_yflx_zface(:,:,:,1,blockID) = fluxy(:,sx:ex,sy+1:ey,sz)
!!$        gr_yflx_zface(:,:,:,2,blockID) = fluxy(:,sx:ex,sy+1:ey,ez)
!!$#endif
!!$#endif
        if (multFluxy) then
           allocate(faceAreas(offx+sx:offx+ex,offy+sy:offy+ey+1,offz+sz:offz+ez))
           areaLeft(sx:,sy:,sz:)  => faceAreas
           
           call Grid_getCellFaceAreas(JAXIS, blkLev, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)

           do presVar = 1,nfluxes
              if (.NOT.(surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,2,1,1+K3D,blockID) == nodetype(blockID))) then
                 where (areaLeft(sx:ex,sy,sz:ez).NE.0.0)
                    flux_y(presVar,:,1,:,blockID) = flux_y(presVar,:,1,:,blockID) * areaLeft(sx:ex,sy,sz:ez)
                 end where
              end if
              if (.NOT.(surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,2,3,1+K3D,blockID) == nodetype(blockID))) then
                 where (areaLeft(sx:ex,ey+1,sz:ez).NE.0.0)
                    flux_y(presVar,:,2,:,blockID) = flux_y(presVar,:,2,:,blockID) * areaLeft(sx:ex,ey+1,sz:ez)
                    
                 end where
              end if
           end do
           deallocate(faceAreas)
        end if
     end if
#endif
     
#if NDIM>2
     if(ztrue) then
        flux_z(:nfluxes,:,:,1,blockID) = fluxz(:,sx:ex,sy:ey,sz) 
        flux_z(:nfluxes,:,:,2,blockID) = fluxz(:,sx:ex,sy:ey,ez+1) 
!!$        gr_zflx(:,:,:,1,blockID) = fluxz(:,sx:ex,sy:ey,sz+1)
!!$        gr_zflx(:,:,:,2,blockID) = fluxz(:,sx:ex,sy:ey,ez)
!!$#ifdef FLASH_HYDRO_UNSPLIT
!!$        !! Store transverse components for the faces in global scratch arrays.
!!$        gr_zflx_xface(:,1,:,:,blockID) = fluxz(:,sx,sy:ey,sz+1:ez)
!!$        gr_zflx_xface(:,2,:,:,blockID) = fluxz(:,ex,sy:ey,sz+1:ez)
!!$        
!!$        gr_zflx_yface(:,:,1,:,blockID) = fluxz(:,sx:ex,sy,sz+1:ez)
!!$        gr_zflx_yface(:,:,2,:,blockID) = fluxz(:,sx:ex,ey,sz+1:ez)
!!$#endif
        if (multFluxz) then
           allocate(faceAreas(offx+sx:offx+ex,offy+sy:offy+ey,offz+sz:offz+ez+1))
           areaLeft(sx:,sy:,sz:)  => faceAreas
           
           call Grid_getCellFaceAreas(KAXIS, blkLev, &
                                  lbound(faceAreas), ubound(faceAreas), &
                                  faceAreas)
           do presVar = 1,nfluxes
              if (.NOT.(surr_blks(1,2,2,1,blockID) > 0 .AND. &
                   surr_blks(3,2,2,1,blockID) == nodetype(blockID))) then
                 flux_z(presVar,:,:,1,blockID) = flux_z(presVar,:,:,1,blockID) * areaLeft(sx:ex,sy:ey,sz)
              end if
              if (.NOT.(surr_blks(1,2,2,3,blockID) > 0 .AND. &
                   surr_blks(3,2,2,3,blockID) == nodetype(blockID))) then
                 flux_z(presVar,:,:,2,blockID) = flux_z(presVar,:,:,2,blockID) * areaLeft(sx:ex,sy:ey,ez+1)
              end if
           end do
           deallocate(faceAreas)
        end if
     end if
#endif
     nullify(fluxx)
     nullify(fluxy)
     nullify(fluxz)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
  return
end subroutine Grid_putFluxData


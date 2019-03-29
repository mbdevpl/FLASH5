!!****if* source/Grid/GridMain/Chombo/gr_getCellFaceArea
!!
!! NAME
!!  gr_getCellFaceArea
!!
!! SYNOPSIS
!!
!!  gr_getCellFaceArea(integer(IN) :: xb,
!!                     integer(IN) :: xe,
!!                     integer(IN) :: yb,
!!                     integer(IN) :: ye,
!!                     integer(IN) :: zb,
!!                     integer(IN) :: ze,
!!                     integer(IN) :: face,
!!                     integer(IN) :: blockID,
!!                     real(xb:xe,yb:ye,zb:ze),(OUT) :: dataBlock)
!!  
!! DESCRIPTION 
!!  
!!  This routine calculate the face area of the specified face on
!!  the specified cells (bounded by xb:xe,yb:ye,zb:ze) of a given block.
!!
!! ARGUMENTS 
!!
!!  xb        : starting cell along IAXIS
!!  xe        : last cell along IAXIS
!!  yb        : starting cell along JAXIS
!!  ye        : last cell along JAXIS
!!  zb        : starting cell along KAXIS
!!  ze        : last cell along KAXIS
!!  face      : specification of the face on which the area must be computed
!!               valid  values are ILO_FACE,IHI_FACE, JLO_FACE etc.
!!  blockID   : my block number
!!  dataBlock : storage for returning calculated values
!!  beginCount : How index argument values are intended.
!!               EXTERIOR or INTERIOR  relative to the block
!!               GLOBALIDX1            global 1-based per-level indexing
!!
!! NOTES
!!
!!   Currently this routine handles only cartesian geometry
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_getCellFaceArea(xb,xe,yb,ye,zb,ze,face,block,dataBlock,beginCount)
  use Driver_interface, ONLY : Driver_abortFlash
!  use Grid_interface,   ONLY : Grid_getDeltas, &
!                               Grid_getCellCoords, &
!                               Grid_getGeometry
  use block_metadata,   ONLY : block_metadata_t

  implicit none

  type(block_metadata_t), intent(IN) :: block
  integer,                intent(IN) :: xb,xe,yb,ye,zb,ze,face
  real,dimension(xb:xe,yb:ye,zb:ze),intent(OUT)::dataBlock
  real,dimension(MDIM) :: del
  integer,intent(IN) :: beginCount

  dataBlock(:,:,:) = 0.0
  call Driver_abortFlash("[gr_getCellFaceArea] DEPRECATED.  Use Grid_getCellFaceAreas")

!  integer :: geometry
!  real :: areaFactor_1, areaFactor_2, areaFactor_3
!  integer :: i,j,k
!  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
!  real,dimension(:),allocatable :: xCoordLeft, xCoordRight, yCoordLeft, &
!       yCoordRight, zCoordLeft, zCoordRight
!  logical, parameter :: gCell = .true.
!
!  if (beginCount .NE. GLOBALIDX1 .AND. beginCount .NE. DEFAULTIDX) then
!     call Driver_abortFlash("gr_getCellFaceArea: invalid beginCount: only global indexing is supported in Amrex mode")
!  end if
!
!  call Grid_getGeometry(geometry)
!  if (geometry==CARTESIAN) then
!     if(NDIM==1) then
!        dataBlock = 1.00
!     else
!        call Grid_getDeltas(block%level, del)
!        select case(face)
!        case(ILO_FACE,IHI_FACE)
!           dataBlock=del(JAXIS)
!           if(NDIM==3)dataBlock=dataBlock*del(KAXIS)
!        case(JLO_FACE,JHI_FACE)
!           dataBlock=del(IAXIS)
!           if(NDIM==3)dataBlock=dataBlock*del(KAXIS)
!        case(KLO_FACE,KHI_FACE)
!           if(NDIM==3)then
!              dataBlock=del(JAXIS)*del(IAXIS)
!           else
!              call Driver_abortFlash("getCellFaceArea : this face doesn't exit")
!           end if
!        end select
!     end if
!
!  else
!
!     blkLimitsGC = block%limitsGC
!
!     allocate(xCoordLeft (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)),&
!              xCoordRight(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
!     allocate(yCoordLeft (blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)),&
!              yCoordRight(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)))
!     allocate(zCoordLeft (blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
!              zCoordRight(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
!
!     call Grid_getCellCoords(IAXIS, LEFT_EDGE,  block%level, blkLimitsGC(LOW, :), blkLimitsGC(HIGH, :), xCoordLeft)
!     call Grid_getCellCoords(IAXIS, RIGHT_EDGE, block%level, blkLimitsGC(LOW, :), blkLimitsGC(HIGH, :), xCoordRight)
!     call Grid_getCellCoords(JAXIS, LEFT_EDGE,  block%level, blkLimitsGC(LOW, :), blkLimitsGC(HIGH, :), yCoordLeft)
!     call Grid_getCellCoords(JAXIS, RIGHT_EDGE, block%level, blkLimitsGC(LOW, :), blkLimitsGC(HIGH, :), yCoordRight)
!     call Grid_getCellCoords(KAXIS, LEFT_EDGE,  block%level, blkLimitsGC(LOW, :), blkLimitsGC(HIGH, :), zCoordLeft)
!     call Grid_getCellCoords(KAXIS, RIGHT_EDGE, block%level, blkLimitsGC(LOW, :), blkLimitsGC(HIGH, :), zCoordRight)
!
!#define ICOORD_LEFT(I)  xCoordLeft(I)
!#define ICOORD_RIGHT(I) xCoordRight(I)
!#define JCOORD_LEFT(I)  yCoordLeft(I)
!#define JCOORD_RIGHT(I) yCoordRight(I)
!#define KCOORD_LEFT(I)  zCoordLeft(I)
!#define KCOORD_RIGHT(I) zCoordRight(I)
!
!     ! The following code is based on Paramesh3's amr_block_geometry.
!     if (face .eq. ILO_FACE .or. face .eq. IHI_FACE) then
!
!        ! compute cell area of faces perpendicular to first coord axis
!        do k = zb,ze
!           do j = yb,ye
!              do i = xb,xe
!                 if (face .eq. ILO_FACE) then
!                    areaFactor_1 = ICOORD_LEFT(i)
!                 else
!                    areaFactor_1 = ICOORD_RIGHT(i)
!                 end if
!                 areaFactor_2 = 1.
!                 areaFactor_3 = 1.
!
!                 if (geometry .eq. SPHERICAL) then
!                    areaFactor_1 = areaFactor_1**2
!                    if (NDIM.ge.2) then
!                       areaFactor_2 = cos( JCOORD_LEFT(j)   ) - & 
!                            &         cos( JCOORD_RIGHT(j) )
!                    else
!                       areaFactor_2 = 2
!                    end if
!                    if(NDIM.eq.3) then
!                       areaFactor_3 = KCOORD_RIGHT(k) - & 
!                            &         KCOORD_LEFT(k)
!                    else
!                       areaFactor_3 = 2*PI
!                    end if
!
!                 else if (geometry .eq. POLAR) then
!                    if(NDIM.ge.2) then
!                       areaFactor_2 =  JCOORD_RIGHT(j) - & 
!                            &          JCOORD_LEFT(j)
!                    else
!                       areaFactor_2 = 2*PI
!                    end if
!                    if(NDIM.eq.3) then
!                       areaFactor_3 = KCOORD_RIGHT(k) - & 
!                            &         KCOORD_LEFT(k)
!                    end if
!                 else if (geometry .eq. CYLINDRICAL) then ! perp to r
!                    ! INT(dz) * r*INT(d theta)
!                    if(NDIM.ge.2) & 
!                         &       areaFactor_2 =  JCOORD_RIGHT(j) - & 
!                         &                       JCOORD_LEFT(j)
!                    if (NDIM.eq.3) then
!                       areaFactor_3 =  KCOORD_RIGHT(k) - & 
!                            &                          KCOORD_LEFT(k)
!                    else
!                       areaFactor_3 = 2*PI
!                    end if
!
!                 end if
!
!                 dataBlock(i,j,k) = abs(areaFactor_1 * areaFactor_2  & 
!                      &                                              * areaFactor_3)
!
!              enddo
!           enddo
!        enddo
!
!
!     else if (face .eq. JLO_FACE .or. face .eq. JHI_FACE) then
!
!        ! compute cell area of faces perpendicular to second coord axis
!        do k = zb,ze
!           do j = yb,ye
!              do i = xb,xe
!                 areaFactor_1 = 1.
!                 areaFactor_2 = 1.
!                 areaFactor_3 = 1.
!
!                 if (geometry .eq. SPHERICAL) then
!                    areaFactor_1 = (ICOORD_RIGHT(i)-ICOORD_LEFT(i)) & 
!                         &            *(ICOORD_LEFT(i)+ICOORD_RIGHT(i))*.5
!                    if (face .eq. JLO_FACE) then
!                       areaFactor_2 = sin( JCOORD_LEFT(j) )
!                    else
!                       areaFactor_2 = sin( JCOORD_RIGHT(j) )
!                    end if
!                    if(NDIM.eq.3) then
!                       areaFactor_3 = KCOORD_RIGHT(k) - & 
!                            KCOORD_LEFT(k)
!                    else
!                       areaFactor_3 = 2*PI
!                    end if
!
!                 else if (geometry .eq. POLAR) then
!                    areaFactor_1 =  ICOORD_RIGHT(i) - & 
!                         &                     ICOORD_LEFT(i)
!
!                 else if (geometry .eq. CYLINDRICAL) then ! perp to z
!                    ! INT(rdr) * INT(d theta)
!                    areaFactor_1 = ( ICOORD_RIGHT(i)**2 -  & 
!                         &                      ICOORD_LEFT(i)**2 )*.5
!                    if(NDIM.eq.3) then
!                       areaFactor_3 =  KCOORD_RIGHT(k) - & 
!                            KCOORD_LEFT(k)
!                    else
!                       areaFactor_3 = 2*PI
!                    end if
!                 endif
!
!                 dataBlock(i,j,k) = abs(areaFactor_1 * areaFactor_2  & 
!                      &                                              * areaFactor_3)
!
!              enddo
!           enddo
!        enddo
!
!
!     else if (face .eq. KLO_FACE .or. face .eq. KHI_FACE) then
!        ! compute cell area of faces perpendicular to third coord axis
!        do k = zb,ze
!           do j = yb,ye
!              do i = xb,xe
!                 areaFactor_1 = 1.
!                 areaFactor_2 = 1.
!                 areaFactor_3 = 1.
!
!                 if (geometry .eq. SPHERICAL) then
!                    areaFactor_1 = (ICOORD_RIGHT(i)-ICOORD_LEFT(i)) & 
!                         &            *(ICOORD_LEFT(i)+ICOORD_RIGHT(i))*.5
!                    if(NDIM.ge.2) then
!                       areaFactor_2 = JCOORD_RIGHT(j) - & 
!                            JCOORD_LEFT(j)
!                    else
!                       areaFactor_2 = 2*PI
!                    end if
!
!                 else if (geometry .eq. POLAR) then
!                    areaFactor_1 = ( ICOORD_RIGHT(i)**2 -  & 
!                         &                      ICOORD_LEFT(i)**2 )*.5
!                    if(NDIM.ge.2) then
!                       areaFactor_2 =  JCOORD_RIGHT(j) - & 
!                            JCOORD_LEFT(j)
!                    else
!                       areaFactor_2 = 2*PI
!                    end if
!
!                 else if (geometry .eq. CYLINDRICAL) then  ! perp to theta
!                    ! INT(dr) * INT(dz)
!                    areaFactor_1 =  ICOORD_RIGHT(i) - & 
!                         &                     ICOORD_LEFT(i)
!                    if(NDIM.ge.2) & 
!                         &       areaFactor_2 =  JCOORD_RIGHT(j) - & 
!                         &                       JCOORD_LEFT(j)
!                 end if
!
!                 dataBlock(i,j,k) = abs(areaFactor_1 * areaFactor_2  & 
!                      &                                              * areaFactor_3)
!
!              enddo
!           enddo
!        enddo
!     end if
!     !-----------------------
!
!     deallocate(xCoordLeft,xCoordRight)
!     deallocate(yCoordLeft,yCoordRight)
!     deallocate(zCoordLeft,zCoordRight)
!
!  end if
end subroutine gr_getCellFaceArea

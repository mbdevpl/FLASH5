!!****if* source/Grid/GridMain/gr_getCellFaceArea
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
!!
!!***
subroutine gr_getCellFaceArea(xb,xe,yb,ye,zb,ze,face,block,dataBlock,beginCount)
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getDeltas
  use block_metadata, ONLY : block_metadata_t
#include "Flash.h"
  use Grid_data, ONLY : gr_geometry
  use Grid_data, ONLY : gr_lrefineMax, gr_delta, gr_globalDomain
#ifdef FLASH_GRID_UG
  use Grid_data, ONLY : gr_iCoords,gr_jCoords,gr_kCoords,gr_iloGc,gr_ihiGc,&
                         gr_jloGc,gr_jhiGc,gr_kloGc,gr_khiGc
#else
  use gr_specificData, ONLY : gr_oneBlock,gr_iloGc,gr_ihiGc,&
                         gr_jloGc,gr_jhiGc,gr_kloGc,gr_khiGc
#endif
  implicit none
#include "constants.h"

  type(block_metadata_t),intent(IN) :: block
  integer,intent(IN) :: xb,xe,yb,ye,zb,ze,face
  real,dimension(xb:xe,yb:ye,zb:ze),intent(OUT)::dataBlock
  real,dimension(MDIM) :: del
  integer,intent(IN) :: beginCount

  real :: areaFactor_1, areaFactor_2, areaFactor_3
  integer :: blockID
  integer :: level
  integer :: i,j,k

  blockID = block%id
  level   = block%level

  if (gr_geometry==CARTESIAN) then
     if(NDIM==1) then
        dataBlock = 1.00
     else
        call Grid_getDeltas(level,del)
        select case(face)
        case(ILO_FACE,IHI_FACE)
           dataBlock=del(JAXIS)
           if(NDIM==3)dataBlock=dataBlock*del(KAXIS)
        case(JLO_FACE,JHI_FACE)
           dataBlock=del(IAXIS)
           if(NDIM==3)dataBlock=dataBlock*del(KAXIS)
        case(KLO_FACE,KHI_FACE)
           if(NDIM==3)then
              dataBlock=del(JAXIS)*del(IAXIS)
           else
              call Driver_abortFlash("getCellFaceArea : this face doesn't exit")
           end if
        end select
     end if
  else

#ifdef FLASH_GRID_UG

#define ICOORD_LEFT(I)  gr_iCoords(LEFT_EDGE ,I,1)
#define ICOORD_RIGHT(I) gr_iCoords(RIGHT_EDGE,I,1)
#define JCOORD_LEFT(I)  gr_jCoords(LEFT_EDGE ,I,1)
#define JCOORD_RIGHT(I) gr_jCoords(RIGHT_EDGE,I,1)
#define KCOORD_LEFT(I)  gr_kCoords(LEFT_EDGE ,I,1)
#define KCOORD_RIGHT(I) gr_kCoords(RIGHT_EDGE,I,1)

#else

#define ICOORD_OB_LEFT(I) (gr_oneBlock(blockID)%firstAxisCoords(LEFT_EDGE,I))
#define ICOORD_OB_RIGHT(I) (gr_oneBlock(blockID)%firstAxisCoords(RIGHT_EDGE,I))
#define JCOORD_OB_LEFT(I) (gr_oneBlock(blockID)%secondAxisCoords(LEFT_EDGE,I))
#define JCOORD_OB_RIGHT(I) (gr_oneBlock(blockID)%secondAxisCoords(RIGHT_EDGE,I))
#define KCOORD_OB_LEFT(I) (gr_oneBlock(blockID)%thirdAxisCoords(LEFT_EDGE,I))
#define KCOORD_OB_RIGHT(I) (gr_oneBlock(blockID)%thirdAxisCoords(RIGHT_EDGE,I))

#endif

! The following code is based on Paramesh3's amr_block_geometry.
     if (face .eq. ILO_FACE .or. face .eq. IHI_FACE) then
        
        ! compute cell area of faces perpendicular to first coord axis
        do k = zb,ze
           do j = yb,ye
              do i = xb,xe
                 if (face .eq. ILO_FACE) then
                    areaFactor_1 = ICOORD_LEFT(i)
                 else
                    areaFactor_1 = ICOORD_RIGHT(i)
                 end if
                 areaFactor_2 = 1.
                 areaFactor_3 = 1.

                 if (gr_geometry .eq. SPHERICAL) then
                    areaFactor_1 = areaFactor_1**2
                    if (NDIM.ge.2) then
                       areaFactor_2 = cos( JCOORD_LEFT(j)   ) - & 
                            &         cos( JCOORD_RIGHT(j) )
                    else
                       areaFactor_2 = 2.
                    end if
                    if(NDIM.eq.3) then
                       areaFactor_3 = KCOORD_RIGHT(k) - & 
                            &         KCOORD_LEFT(k)
                    else
                       areaFactor_3 = 2.*PI
                    end if

                 else if (gr_geometry .eq. POLAR) then
                    if(NDIM.ge.2) then
                       areaFactor_2 =  JCOORD_RIGHT(j) - & 
                            &          JCOORD_LEFT(j)
                    else
                       areaFactor_2 = 2.*PI
                    end if
                    if(NDIM.eq.3) then
                       areaFactor_3 = KCOORD_RIGHT(k) - & 
                            &         KCOORD_LEFT(k)
                    end if
                 else if (gr_geometry .eq. CYLINDRICAL) then ! perp to r
                    ! INT(dz) * r*INT(d theta)
                    if(NDIM.ge.2) & 
                         &       areaFactor_2 =  JCOORD_RIGHT(j) - & 
                         &                       JCOORD_LEFT(j)
                    if (NDIM.eq.3) then
                       areaFactor_3 =  KCOORD_RIGHT(k) - & 
                            &                          KCOORD_LEFT(k)
                    else
                       areaFactor_3 = 2.*PI
                    end if

                 end if

                 dataBlock(i,j,k) = abs(areaFactor_1 * areaFactor_2  & 
                      &                                              * areaFactor_3)

              enddo
           enddo
        enddo


     else if (face .eq. JLO_FACE .or. face .eq. JHI_FACE) then

        ! compute cell area of faces perpendicular to second coord axis
        do k = zb,ze
           do j = yb,ye
              do i = xb,xe
                 areaFactor_1 = 1.
                 areaFactor_2 = 1.
                 areaFactor_3 = 1.

                 if (gr_geometry .eq. SPHERICAL) then
                    areaFactor_1 = (ICOORD_RIGHT(i)-ICOORD_LEFT(i)) & 
                         &            *(ICOORD_LEFT(i)+ICOORD_RIGHT(i))*.5
                    if (face .eq. JLO_FACE) then
                       areaFactor_2 = sin( JCOORD_LEFT(j) )
                    else
                       areaFactor_2 = sin( JCOORD_RIGHT(j) )
                    end if
                    if(NDIM.eq.3) then
                       areaFactor_3 = KCOORD_RIGHT(k) - & 
                                      KCOORD_LEFT(k)
                    else
                       areaFactor_3 = 2.*PI
                    end if

                 else if (gr_geometry .eq. POLAR) then
                    areaFactor_1 =  ICOORD_RIGHT(i) - & 
                         &                     ICOORD_LEFT(i)

                 else if (gr_geometry .eq. CYLINDRICAL) then ! perp to z
                    ! INT(rdr) * INT(d theta)
                    areaFactor_1 = ( ICOORD_RIGHT(i)**2 -  & 
                         &                      ICOORD_LEFT(i)**2 )*.5
                    if(NDIM.eq.3) then
                       areaFactor_3 =  KCOORD_RIGHT(k) - & 
                                       KCOORD_LEFT(k)
                    else
                       areaFactor_3 = 2.*PI
                    end if
                 endif

                 dataBlock(i,j,k) = abs(areaFactor_1 * areaFactor_2  & 
                      &                                              * areaFactor_3)

              enddo
           enddo
        enddo


     else if (face .eq. KLO_FACE .or. face .eq. KHI_FACE) then
        ! compute cell area of faces perpendicular to third coord axis
        do k = zb,ze
           do j = yb,ye
              do i = xb,xe
                 areaFactor_1 = 1.
                 areaFactor_2 = 1.
                 areaFactor_3 = 1.

                 if (gr_geometry .eq. SPHERICAL) then
                    areaFactor_1 = (ICOORD_RIGHT(i)-ICOORD_LEFT(i)) & 
                         &            *(ICOORD_LEFT(i)+ICOORD_RIGHT(i))*.5
                    if(NDIM.ge.2) then
                       areaFactor_2 = JCOORD_RIGHT(j) - & 
                                      JCOORD_LEFT(j)
                    else
                       areaFactor_2 = 2.*PI
                    end if

                 else if (gr_geometry .eq. POLAR) then
                    areaFactor_1 = ( ICOORD_RIGHT(i)**2 -  & 
                         &                      ICOORD_LEFT(i)**2 )*.5
                    if(NDIM.ge.2) then
                       areaFactor_2 =  JCOORD_RIGHT(j) - & 
                                       JCOORD_LEFT(j)
                    else
                       areaFactor_2 = 2.*PI
                    end if

                 else if (gr_geometry .eq. CYLINDRICAL) then  ! perp to theta
                    ! INT(dr) * INT(dz)
                    areaFactor_1 =  ICOORD_RIGHT(i) - & 
                         &                     ICOORD_LEFT(i)
                    if(NDIM.ge.2) & 
                         &       areaFactor_2 =  JCOORD_RIGHT(j) - & 
                         &                       JCOORD_LEFT(j)
                 end if

                 dataBlock(i,j,k) = abs(areaFactor_1 * areaFactor_2  & 
                      &                                              * areaFactor_3)

              enddo
           enddo
        enddo
     end if
     !-----------------------


  end if
contains
  pure real function ICOORD_left(i)
    integer,VALUE,INTENT(IN) :: i
    integer       :: i0
    if ((beginCount==EXTERIOR .OR. beginCount==INTERIOR) .AND. blockID > 0) then
       ICOORD_left = ICOORD_OB_LEFT(i)
    else
       i0 = i-1
       if (beginCount==EXTERIOR .OR. beginCount==INTERIOR) i0 = i + block%limitsGC(LOW,IAXIS) - 2
       ICOORD_left = gr_globalDomain(LOW,IAXIS) + i0 * gr_delta(IAXIS,level)
    end if
  end function ICOORD_left
  pure real function ICOORD_right(i)
    integer,VALUE,INTENT(IN) :: i
    integer       :: i1
    if ((beginCount==EXTERIOR .OR. beginCount==INTERIOR) .AND. blockID > 0) then
       ICOORD_right = ICOORD_OB_RIGHT(i)
    else
       i1 = i
       if (beginCount==EXTERIOR .OR. beginCount==INTERIOR) i1 = i + block%limitsGC(LOW,IAXIS) - 1
       ICOORD_right = gr_globalDomain(LOW,IAXIS) + i1 * gr_delta(IAXIS,level)
    end if
  end function ICOORD_right
  pure real function JCOORD_left(i)
    integer,VALUE,INTENT(IN) :: i
    integer       :: i0
    if ((beginCount==EXTERIOR .OR. beginCount==INTERIOR) .AND. blockID > 0) then
       JCOORD_left = JCOORD_OB_LEFT(i)
    else
       i0 = i-1
       if (beginCount==EXTERIOR .OR. beginCount==INTERIOR) i0 = i + block%limitsGC(LOW,JAXIS) - 2
       JCOORD_left = gr_globalDomain(LOW,JAXIS) + i0 * gr_delta(JAXIS,level)
    end if
  end function JCOORD_left
  pure real function JCOORD_right(i)
    integer,VALUE,INTENT(IN) :: i
    integer       :: i1
    if ((beginCount==EXTERIOR .OR. beginCount==INTERIOR) .AND. blockID > 0) then
       JCOORD_right = JCOORD_OB_RIGHT(i)
    else
       i1 = i
       if (beginCount==EXTERIOR .OR. beginCount==INTERIOR) i1 = i + block%limitsGC(LOW,JAXIS) - 1
       JCOORD_right = gr_globalDomain(LOW,JAXIS) + i1 * gr_delta(JAXIS,level)
    end if
  end function JCOORD_right
  pure real function KCOORD_left(i)
    integer,VALUE,INTENT(IN) :: i
    integer       :: i0
    if ((beginCount==EXTERIOR .OR. beginCount==INTERIOR) .AND. blockID > 0) then
       KCOORD_left = KCOORD_OB_LEFT(i)
    else
       i0 = i-1
       if (beginCount==EXTERIOR .OR. beginCount==INTERIOR) i0 = i + block%limitsGC(LOW,KAXIS) - 2
       KCOORD_left = gr_globalDomain(LOW,KAXIS) + i0 * gr_delta(KAXIS,level)
    end if
  end function KCOORD_left
  pure real function KCOORD_right(i)
    integer,VALUE,INTENT(IN) :: i
    integer       :: i1
    if ((beginCount==EXTERIOR .OR. beginCount==INTERIOR) .AND. blockID > 0) then
       KCOORD_right = KCOORD_OB_RIGHT(i)
    else
       i1 = i
       if (beginCount==EXTERIOR .OR. beginCount==INTERIOR) i1 = i + block%limitsGC(LOW,KAXIS) - 1
       KCOORD_right = gr_globalDomain(LOW,KAXIS) + i1 * gr_delta(KAXIS,level)
    end if
  end function KCOORD_right
end subroutine gr_getCellFaceArea

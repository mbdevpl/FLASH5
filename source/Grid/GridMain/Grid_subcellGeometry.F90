!!****if* source/Grid/GridMain/Grid_subcellGeometry
!!
!! NAME
!!
!!  Grid_subcellGeometry
!!
!! SYNOPSIS
!!
!!  call Grid_subcellGeometry(integer,VALUE(in) :: nsubi,
!!                            integer,VALUE(in) :: nsubj,
!!                            integer,VALUE(in) :: nsubk,
!!                            real(in) :: dvcell,
!!                            real, intent(OUT)  :: dvsub(nsubI,nsubJ),
!!                            real,OPTIONAL(in) :: xl,
!!                            real,OPTIONAL(in) :: xr,
!!                            real,OPTIONAL(in) :: yl,
!!                            real,OPTIONAL(in) :: yr,
!!                            integer,OPTIONAL(in) :: pos(*),
!!                            integer,OPTIONAL(in) :: blockID)
!!
!! DESCRIPTION
!!
!!  Geometrically correct computation of the volumes of subcells.
!!
!!
!! ARGUMENTS
!!
!!   nsubi :     Number of subcell lengths per cell length in X direction.
!!
!!   nsubj :     Number of subcell lengths per cell length in Y direction.
!!
!!   nsubk :     Number of subcell lengths per cell length in Z direction.
!!
!!   dvcell :    Volume of the whole cell.
!!
!!   dvsub :     Volumes of subcells. Note that this is a 2-dimensional array.
!!               For the geometries that FLASH supports, subcells volumes
!!               never depend on the third coordinate.
!!
!!   xl :        X-coordinate of left cell face.
!!
!!   xr :        X-coordinate of right cell face.
!!
!!   yl :        Y-coordinate of lower cell face.
!!
!!   yr :        Y-coordinate of upper cell face.
!!
!!   pos :       currently unused.
!!
!!   blockID : ID of block in current processor, currently unused
!!
!!
!!***

#include "FortranLangFeatures.fh"
#include "constants.h"
#include "Flash.h"

subroutine Grid_subcellGeometry(nsubI,nsubJ,nsubK, &
     dvCell, dvSub, xL,xR, yL,yR, pos, blockID)
  use Grid_interface, ONLY : Grid_getGeometry
  implicit none
  integer, VALUE_INTENT(IN) :: nsubI, nsubJ, nsubK
  real,    intent(in)  :: dvCell
  real,    intent(OUT) :: dvSub(nsubI, nsubJ)
  real,OPTIONAL,intent(in) :: xL, xR
  real,OPTIONAL,intent(in) :: yL, yR
  integer,OPTIONAL, intent(in) :: blockID
  integer,OPTIONAL, intent(in) :: pos(*)

  integer              :: geometry
  real                 :: dvs
  real                 :: xccN2inv, xcsN
  real                 :: x2ccNinv3, x32csN
  real                 :: xLs, xRs
  real                 :: yLs, yRs
  integer              :: i, j, effGeometry

  call Grid_getGeometry(geometry)
  effGeometry = geometry
  if (nsubI == 1) then
     if (geometry .NE. SPHERICAL) effGeometry = CARTESIAN
     if (nsubJ == 1) effGeometry = CARTESIAN
  end if
  select case (effGeometry)
     case(CARTESIAN)
        dvs          = dvCell
        dvs          = dvs   / real(nsubI)
        dvs          = dvs   / real(nsubJ)
        dvSub(:,:)   = dvs   / real(nsubK)
     case(CYLINDRICAL,POLAR)
!        xcc          =           (xL+xR) * 0.5
!        xccN         =           (xL+xR) * 0.5 * real(nsubI)
!        xccNinv      = 2.0 / (real(nsubI) * (xL+xR))
        xccN2inv     = 2.0 / (real(nsubI)*real(nsubI) * (xL+xR))
!        dxiInv = real(nsubI) / (xR - xL)
!        dxsInv = real(nsubI) / (xR - xL)
        do i=1,nsubI
!           xLs =    ( xL * (nsubI+1-i)  +  xR * (i-1) ) * dxsInv
!           xRs =    ( xL * (nsubI  -i)  +  xR * (i  ) ) * dxsInv
!!$           xcs =    ( xL * (real(nsubI-i)+0.5)  +  xR * (real(i)-0.5) ) / real(nsubI)
!!$           dvSub(i,1) =dvCell *  xcs / xcc / real(nsubI)
!           xcs =    ( xL * (real(nsubI-i)+0.5)  +  xR * (real(i)-0.5) ) / real(nsubI)
           xcsN =    ( xL * (real(nsubI-i)+0.5)  +  xR * (real(i)-0.5) )
           dvSub(i,1) =dvCell *  xcsN * xccN2inv
           dvSub(i,1) = dvSub(i,1) / real(nsubJ)
           dvSub(i,:) = dvSub(i,1) / real(nsubK)
        end do
     case(SPHERICAL)
!        dxiInv = real(nsubI) / (xR - xL)
!        dxsInv = real(nsubI) / (xR - xL)
!!$           x2cc          =           ((xL+xR)**2-xL*xR) * 0.5
!!$           x2ccN         =           ((xL+xR)**2-xL*xR) * 0.5 * real(nsubI)
        x2ccNinv3      = 1.0 / (real(nsubI) * ((xL+xR)**2-xL*xR))
!!$           x2ccN2inv3     = 1.0 / (real(nsubI)*real(nsubI) * ((xL+xR)**2-xL*xR))
        do i=1,nsubI
           xLs =    ( xL * (nsubI+1-i)  +  xR * (i-1) ) / real(nsubI)
           xRs =    ( xL * (nsubI  -i)  +  xR * (i  ) ) / real(nsubI)
!!$              xLsn =    ( xL * (nsubI+1-i)  +  xR * (i-1) ) * dxsInv
!!$              xRsn =    ( xL * (nsubI  -i)  +  xR * (i  ) ) * dxsInv
!           xcs =    ( xL * (real(nsubI-i)+0.5)  +  xR * (real(i)-0.5) ) / real(nsubI)

           x32csN =    ((xLs+xRs)**2-xLs*xRs)
           dvSub(i,1) =dvCell *  x32csN * x2ccNinv3
        end do
        if(NDIM > 1 .OR. nsubJ > 1) then
           do j=1,nsubJ
              yLs =    ( yL * (nsubJ+1-j)  +  yR * (j-1) ) / real(nsubJ)
              yRs =    ( yL * (nsubJ  -j)  +  yR * (j  ) ) / real(nsubJ)

              dvSub(:,j) = dvSub(:,1) * (cos(yLs)-cos(yRs)) / ((cos(yL)-cos(yR)) * real(nsubK))
           end do
        else
           dvSub(:,1) = dvSub(:,1) / real(nsubK)
        end if


     case default
        call Driver_abortFlash('[Grid_subcellGeometry] Not implemented for current geometry')
  end select


end subroutine Grid_subcellGeometry

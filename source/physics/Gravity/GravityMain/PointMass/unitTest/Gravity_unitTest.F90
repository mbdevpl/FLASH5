!!****if* source/physics/Gravity/GravityMain/PointMass/unitTest/Gravity_unitTest
!! NAME
!!
!!  Gravity_unitTest
!!
!! SYNOPSIS
!!
!!  call Gravity_unitTest(integer(IN) :: fileUnit,
!!                       logical(OUT) :: perfect)
!!
!! DESCRIPTION
!!
!! This function is a generic unit test for the Gravity unit that is
!! suitable for constant PointMass gravity. It is invoked in
!! the setup unitTest/Gravity/PointMass.
!!
!!  ARGUMENTS
!!
!!
!!   fileUnit : unit number for file opened by the unitTest/Gravity setup
!!              in which to write results of the test
!!
!!   perfect : indicates test ran without error is true.
!!
!!  PARAMETERS
!!
!!   ptxpos  -  X-coordinate of the point mass (IAXIS);
!!              geometrically, this is R or r in non-Cartesian geometries.
!!   ptypos  -  Y-coordinate of the point mass (JAXIS);
!!              geometrically, this is z or theta or phi in
!!              cylindrical or spherical or polar geometries, respectively.
!!   ptzpos  -  Z-coordinate of the point mass (KAXIS)
!!
!!  NOTES
!!
!!   This test (and the PointMass implementation of Gravity) is
!!   probably only meaningful under the following circumstances:
!!   * 3D Cartesian;
!!   * 2D Cartesian, with ptzpos ==(forced) 0.0;
!!   * 1D Cartesian, with ptzpos as above and ptypos ==(forced) 0.0;
!!   * 2D cylindrical, with ptxpos == 0.0 and ptzpos ignored;
!!   * ptxpos == ptypos == ptzpos in all other cases.
!!***

!!REORDER(4): Uin

#include "constants.h"
#include "Flash.h"

subroutine Gravity_unitTest( fileUnit, perfect)

  use Logfile_interface, ONLY:  Logfile_stampMessage
  use Grid_iterator,       ONLY : Grid_iterator_t
  use Grid_tile,         ONLY:  Grid_tile_t
  use Grid_interface,    ONLY:  Grid_getCellCoords,  &
                                  Grid_getTileIterator, &
                                  Grid_releaseTileIterator, &
                                  Grid_getGeometry
  use Gravity_interface, ONLY:  Gravity_accelOneRow

  use Gravity_data,      ONLY:  grv_factor,          &
                                grv_ptxpos,          &
                                grv_ptypos,          &
                                grv_ptzpos

  implicit none

  integer, intent(in) :: fileUnit
  logical, intent(out) :: perfect

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc
  real,pointer, dimension(:,:,:,:) :: Uin
  real,allocatable, dimension(:)   :: xCell, yCell, zCell, xAccel
  real,             dimension(1)   ::                              yAccel, zAccel
  real                             :: ovec(MDIM), dist
  real                             :: accelVec(MDIM), accelAbs
  real                             :: expectedAccelVec(MDIM) !, impliedFactor

  integer                  :: geometry
  integer                  :: i, j, k
  integer,dimension(MDIM)  :: lo, hi
  integer                  :: maxAccelDir

  perfect = .TRUE.

  call Grid_getGeometry(geometry)

  maxAccelDir = NDIM
  if (geometry == CYLINDRICAL) then
     maxAccelDir = min(2,maxAccelDir)
     perfect = (grv_ptxpos == 0.0)
     if (maxAccelDir == 1) &
          perfect = (perfect .AND. (grv_ptypos == 0.0))
  else if (geometry .NE. CARTESIAN) then
     maxAccelDir = 1
     perfect = (grv_ptxpos == 0.0)
  else                            !CARTESIAN
     if (maxAccelDir == 1) &
          perfect = (perfect .AND. (grv_ptypos == 0.0))
     if (maxAccelDir < 3) &
          perfect = (perfect .AND. (grv_ptzpos == 0.0))
  end if



  call Grid_getTileIterator(itor, LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     call tileDesc%getDataPtr(Uin, CENTER)

     ! Get the coordinate information for the current tile

     lo=tileDesc%limits(LOW,:)
     hi=tileDesc%limits(HIGH,:)

     allocate(xCell (lo(IAXIS):hi(IAXIS)))
     allocate(xAccel(lo(IAXIS):hi(IAXIS)))
     allocate(yCell (lo(JAXIS):hi(JAXIS)))
     allocate(zCell (lo(KAXIS):hi(KAXIS)))
     if (NDIM == 3) then
        call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, lo, hi, zCell)
     endif
     if (NDIM >= 2) then
        call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, lo, hi, yCell)
     endif
     call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xCell)


     do       k = lo(KAXIS), hi(KAXIS)
        do    j = lo(JAXIS), hi(JAXIS)
           call Gravity_accelOneRow((/j,k/), SWEEP_X, tileDesc, &
                lo(IAXIS), hi(IAXIS), xAccel, &
                Uin, -1, (/-1,-1,-1/))

           do i = lo(IAXIS), hi(IAXIS)

              call Gravity_accelOneRow((/i,k/), SWEEP_Y, tileDesc, &
                   j, j, yAccel, &
                   Uin)
              call Gravity_accelOneRow((/i,j/), SWEEP_Z, tileDesc, &
                   k, k, zAccel)

              accelVec = (/xAccel(i), yAccel(1), zAccel(1)/)
              accelAbs = NORM2(accelVec)
              ovec     = (/xCell(i) - grv_ptxpos, &
                           yCell(j) - grv_ptypos, &
                           zCell(k) - grv_ptzpos/)
              dist     = NORM2(ovec(:maxAccelDir))
              expectedAccelVec = grv_factor * ovec / dist**3

              if (ANY(abs(accelVec(:maxAccelDir) - expectedAccelVec(:maxAccelDir)) &
                      > 1.e-15 * abs(grv_factor /(max(dist**2,1.e-200))))) then
                 print*,'Large error:'
                 perfect = .FALSE.
!                 impliedFactor = accelAbs*dist**2
999              format(1x,3(i3),2x,4(0Pf11.8,1x,3(1PG15.6,:1x)))
                 print 999,i,j,k,dist,ovec,accelAbs,accelVec,grv_factor/dist**2, expectedAccelVec, &
                                 -accelAbs - grv_factor / dist**2, accelVec - expectedAccelVec
              end if

           end do
        end do
     end do

     deallocate(xCell)
     deallocate(xAccel)
     deallocate(yCell)
     deallocate(zCell)
     call tileDesc%releaseDataPtr(Uin, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)




  if (.NOT. perfect) then
     write(*,*)"[Gravity_unitTest] Failure.  Gravity_accelOneRow was NOT tested successfully."
     write(fileUnit,*)"[Gravity_unitTest] Failure.  Gravity_accelOneRow was NOT tested successfully."
     call Logfile_stampMessage("[Gravity_unitTest] Failure.  Gravity_accelOneRow was NOT tested successfully.")
  end if

  return
end subroutine Gravity_unitTest





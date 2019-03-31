!!****f* source/Grid/GridMain/Grid_getCellFaceAreas
!!
!! NAME
!!  Grid_getCellFaceAreas
!!
!! NOTES
!!  Please refer to the stub file for complete documentation of this routine.
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Grid_getCellFaceAreas(axis, level, lo, hi, areas)
   use Driver_interface, ONLY : Driver_abortFlash
   use Grid_interface,   ONLY : Grid_getDeltas, &
                                Grid_getCellCoords
   use Grid_data,        ONLY : gr_geometry

   integer, intent(IN)  :: axis
   integer, intent(IN)  :: level
   integer, intent(IN)  :: lo(1:MDIM)
   integer, intent(IN)  :: hi(1:MDIM)
   real,    intent(OUT) :: areas(lo(IAXIS):hi(IAXIS), &
                                 lo(JAXIS):hi(JAXIS), &
                                 lo(KAXIS):hi(KAXIS))

   real    :: deltas(1:MDIM)
   integer :: loCell(1:MDIM)
   integer :: hiCell(1:MDIM)

   real, allocatable :: faceCoords(:), rf(:)
#if NDIM >= 2
   real, allocatable :: thf(:)
#endif

   real    :: area, facebase
   integer :: i, j, k

   if (      (gr_geometry /= CARTESIAN) &
       .AND. (gr_geometry /= SPHERICAL   .OR. NDIM > 2) &
       .AND. (gr_geometry /= CYLINDRICAL .OR. NDIM /= 2)) then
     areas(:, :, :) = 0.0
     call Driver_abortFlash("[Grid_getCellFaceAreas] Not tested yet")
   else if ((axis /= IAXIS) .AND. (axis /= JAXIS) .AND. (axis /= KAXIS)) then
     call Driver_abortFlash("[Grid_getCellFaceAreas] Invalid axis")
   end if

   call Grid_getDeltas(level, deltas)

   select case (gr_geometry)
   case (CARTESIAN)
      associate(dx => deltas(IAXIS), &
                dy => deltas(JAXIS), &
                dz => deltas(KAXIS))
#if   NDIM == 1
         area = 1.0
#elif NDIM == 2
         if      (axis == IAXIS) then
            area = dy
         else if (axis == JAXIS) then
            area = dx
         else
            ! DEV: TODO Should this set area to 1?
            call Driver_abortFlash("[Grid_getCellFaceAreas] Invalid axis for 2D")
         end if
#elif NDIM == 3
         if      (axis == IAXIS) then
            area = dy * dz
         else if (axis == JAXIS) then
            area = dx * dz
         else if (axis == KAXIS) then
            area = dx * dy
         end if
#endif

         do       k = lo(KAXIS), hi(KAXIS)
            do    j = lo(JAXIS), hi(JAXIS)
               do i = lo(IAXIS), hi(IAXIS)
                  areas(i, j, k) = area
               end do
            end do
         end do
      end associate
   case (CYLINDRICAL)
      if      (axis == IAXIS) then
         ! Get radial coordinate of face centers
         allocate(faceCoords(lo(IAXIS):hi(IAXIS)))

         ! Convert face indices to indices of cells associated with faces 
         loCell(:) = lo(:)
         hiCell(:) = hi(:)
         hiCell(axis) = hiCell(axis) - 1
         call Grid_getCellCoords(IAXIS, FACES, level, loCell, hiCell, faceCoords)

         associate(dz   => deltas(JAXIS), &
                   dPhi => deltas(KAXIS), &
                   r    => faceCoords)
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)
#if   NDIM == 1
                     areas(i, j, k) = 2.0 * PI * ABS(r(i))
#elif NDIM == 2
                     areas(i, j, k) = 2.0 * PI * ABS(r(i)) * dz
#elif NDIM == 3
                     areas(i, j, k) = ABS(r(i)) * dz * dPhi
#endif
                  end do
               end do
            end do
         end associate
         deallocate(faceCoords)
      else if (axis == JAXIS) then
         ! Get radial coordinate of face centers
         allocate(faceCoords(lo(IAXIS):hi(IAXIS)+1))

         ! Convert face indices to indices of cells associated with faces 
         loCell(:) = lo(:)
         hiCell(:) = hi(:)
         hiCell(axis) = hiCell(axis) - 1
         call Grid_getCellCoords(IAXIS, FACES, level, loCell, hiCell, faceCoords)

         associate(dPhi => deltas(KAXIS), &
                   r    => faceCoords)
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)
                     ! DEV: TODO These can be done more simply using
                     ! the radii of the cell centers (see cell volumes)
#if   NDIM == 1
                     areas(i, j, k) = PI * (r(i+1)**2 - r(i)**2)
#elif NDIM == 2
                     areas(i, j, k) = PI * (r(i+1)**2 - r(i)**2)
#elif NDIM == 3
                     areas(i, j, k) = 0.5 * (r(i+1)**2 - r(i)**2) * dPhi
#endif
                  end do
               end do
            end do
         end associate
         deallocate(faceCoords)
      else if (axis == KAXIS) then
         associate(dr   => deltas(IAXIS), &
                   dz   => deltas(JAXIS))
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)
#if   NDIM == 1
                     areas(i, j, k) = dr
#elif NDIM == 2
                     areas(i, j, k) = dr * dz
#elif NDIM == 3
                     areas(i, j, k) = dr * dz
#endif
                  end do
               end do
            end do
         end associate
      end if
   case (SPHERICAL)
      if      (axis == IAXIS) then
         ! Get coordinates of faces
         allocate(faceCoords(lo(axis):hi(axis)))

         ! Convert face indices to indices of cells associated with faces
         loCell(:) = lo(:)
         hiCell(:) = hi(:)
         hiCell(axis) = hiCell(axis) - 1
         call Grid_getCellCoords(axis, FACES, level, loCell, hiCell, faceCoords)
#if NDIM >= 2
         allocate(thf         (lo(JAXIS):hi(JAXIS)+1))
         call Grid_getCellCoords(JAXIS, FACES, level, lo, hi, thf )
#endif

         associate(dPhi => deltas(KAXIS), &
                   r    => faceCoords)
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)
                     facebase = r(i) * r(i)
#if   NDIM == 1
                     areas(i, j, k) = facebase * 4.0 * PI
#elif NDIM == 2
                     areas(i, j, k) = facebase * ( cos(thf(i)) - cos(thf(i+1)) ) * 2.0 * PI
#elif NDIM == 3
                     areas(i, j, k) = facebase * ( cos(thf(i)) - cos(thf(i+1)) ) * dPhi
#endif
                  end do
               end do
            end do
         end associate
         deallocate(faceCoords)
#if NDIM >= 2
         deallocate(thf)
#endif
      else if (axis == JAXIS) then
         ! Get r-coordinates of r-faces
         allocate(rf          (lo(IAXIS):hi(IAXIS)+1))
         call Grid_getCellCoords(IAXIS, FACES, level, lo, hi, rf  )
         ! Get coordinates of faces
         allocate(faceCoords(lo(axis):hi(axis)))

         ! Convert face indices to indices of cells associated with faces 
         loCell(:) = lo(:)
         hiCell(:) = hi(:)
         hiCell(axis) = hiCell(axis) - 1
         call Grid_getCellCoords(axis, FACES, level, loCell, hiCell, faceCoords)

         associate(dPhi => deltas(KAXIS), &
                   thf  => faceCoords)
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)
                     facebase = (rf(i)+rf(i+1))*(rf(i+1)-rf(i))*0.5
#if   NDIM == 1
                     areas(i, j, k) = facebase * 2.0 * PI
#elif NDIM == 2
                     areas(i, j, k) = facebase * sin(thf(j)) * 2.0 * PI
#elif NDIM == 3
                     areas(i, j, k) = facebase * sin(thf(j)) * dPhi
#endif
                  end do
               end do
            end do
         end associate
         deallocate(rf)
         deallocate(faceCoords)
      else if (axis == KAXIS) then
         ! Get r-coordinates of faces
         allocate(faceCoords(lo(IAXIS):hi(IAXIS)))
         call Grid_getCellCoords(axis, FACES, level, lo, hi, faceCoords)
#if NDIM >= 2
         allocate(thf         (lo(JAXIS):hi(JAXIS)+1))
         call Grid_getCellCoords(JAXIS, FACES, level, lo, hi, thf )
#endif
         associate(dr   => deltas(IAXIS), &
                   rf   => faceCoords)
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)
                     facebase = dr *  &
                        ( rf(i  ) *  rf(i  )  +  &
                          rf(i  ) *  rf(i+1)  +  &
                          rf(i+1) *  rf(i+1) )
#if   NDIM == 1
                     areas(i, j, k) = facebase * 2./3.
#elif NDIM >= 2
                     areas(i, j, k) = facebase * ( cos(thf(i)) - cos(thf(i+1)) ) / 3.0
#endif
                  end do
               end do
            end do
         end associate
         deallocate(facecoords)
#if NDIM >= 2
         deallocate(thf)
#endif
      end if
   end select
end subroutine Grid_getCellFaceAreas


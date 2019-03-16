!!****f* source/Grid/GridMain/Grid_getCellVolumes
!!
!! NAME
!!  Grid_getCellVolumes
!!
!! NOTES
!!  Please refer to the stub implementation of this routine for the full
!!  documentation of this routine.
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Grid_getCellVolumes(level, lo, hi, volumes)
   use Driver_interface, ONLY : Driver_abortFlash
   use Grid_interface,   ONLY : Grid_getDeltas, &
                                Grid_getCellCoords
   use Grid_data,        ONLY : gr_geometry
 
   integer, intent(IN)  :: level
   integer, intent(IN)  :: lo(1:MDIM)
   integer, intent(IN)  :: hi(1:MDIM)
   real,    intent(OUT) :: volumes(lo(IAXIS):hi(IAXIS), &
                                   lo(JAXIS):hi(JAXIS), &
                                   lo(KAXIS):hi(KAXIS))

   real, allocatable :: centerCoords(:)

   real    :: deltas(1:MDIM)
   integer :: i, j, k

   if (      (gr_geometry /= CARTESIAN) &
       .AND. (gr_geometry /= CYLINDRICAL .OR. NDIM /= 2)) then
     volumes(:, :, :) = 0.0
     call Driver_abortFlash("[Grid_getCellVolumes] Not tested yet")
   end if

   call Grid_getDeltas(level, deltas)

   select case (gr_geometry)
   case (CARTESIAN)
      associate(dx => deltas(IAXIS), &
                dy => deltas(JAXIS), &
                dz => deltas(KAXIS))
         do       k = lo(KAXIS), hi(KAXIS)
            do    j = lo(JAXIS), hi(JAXIS)
               do i = lo(IAXIS), hi(IAXIS)
#if   NDIM == 1
                  volumes(i, j, k) = dx
#elif NDIM == 2
                  volumes(i, j, k) = dx * dy
#elif NDIM == 3
                  volumes(i, j, k) = dx * dy * dz
#endif
               end do
            end do
         end do
      end associate
   case (CYLINDRICAL)
      allocate(centerCoords(lo(IAXIS):hi(IAXIS)))
      call Grid_getCellCoords(IAXIS, CENTER, level, lo, hi, centerCoords)

      associate(dr   => deltas(IAXIS), &
                dz   => deltas(JAXIS), &
                dPhi => deltas(KAXIS), &
                r    => centerCoords)
         do       k = lo(KAXIS), hi(KAXIS)
            do    j = lo(JAXIS), hi(JAXIS)
               do i = lo(IAXIS), hi(IAXIS)
#if   NDIM == 1
                  volumes(i, j, k) = 2.0 * PI * ABS(r(i)) * dr
#elif NDIM == 2
                  volumes(i, j, k) = 2.0 * PI * ABS(r(i)) * dr * dz
#elif NDIM == 3
                  volumes(i, j, k) = ABS(r(i)) * dr * dz * dPhi
#endif
               end do
            end do
         end do
      end associate

      deallocate(centerCoords)
   end select
end subroutine Grid_getCellVolumes


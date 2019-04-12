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

   real, allocatable :: centerCoords(:), rf(:)
#if NDIM >= 2
   real, allocatable :: thf(:)
#endif

   real    :: deltas(1:MDIM)
   real    :: cellvolume
   integer :: i, j, k

   if (.NOT.(     (gr_geometry == CARTESIAN)                   &
             .OR. (gr_geometry == SPHERICAL   .AND. NDIM < 3)  &
             .OR. (gr_geometry == CYLINDRICAL .AND. NDIM == 2) &
             )     ) then
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

   case (SPHERICAL)
      allocate(rf          (lo(IAXIS):hi(IAXIS)+1))
      call Grid_getCellCoords(IAXIS, FACES, level, lo, hi, rf  )
#if NDIM >= 2
      allocate(thf         (lo(JAXIS):hi(JAXIS)+1))
      call Grid_getCellCoords(JAXIS, FACES, level, lo, hi, thf )
#endif

      associate(dr   => deltas(IAXIS), &
                dPhi => deltas(KAXIS))
         do       k = lo(KAXIS), hi(KAXIS)
            do    j = lo(JAXIS), hi(JAXIS)
               do i = lo(IAXIS), hi(IAXIS)
                  cellvolume = dr *  &
                        ( rf(i  ) *  rf(i  )  +  &
                          rf(i  ) *  rf(i+1)  +  &
                          rf(i+1) *  rf(i+1) )
#if   NDIM == 1
                  volumes(i, j, k) = cellvolume * 4.*PI/3.
#elif NDIM == 2
                  volumes(i, j, k) = cellvolume * ( cos(thf(j)) - cos(thf(j+1)) ) * 2.*PI/3.
#elif NDIM == 3
                  volumes(i, j, k) = cellvolume * ( cos(thf(j)) - cos(thf(j+1)) ) *   &
                                     dPhi / 3.0
#endif
               end do
            end do
         end do
      end associate

      deallocate(rf)
#if NDIM >= 2
      deallocate(thf)
#endif
   end select
end subroutine Grid_getCellVolumes


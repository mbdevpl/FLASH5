!!****f* source/Grid/Grid_getCellFaceAreas
!!
!! NAME
!!  Grid_getCellFaceAreas
!!
!! SYNOPSIS
!!  call Grid_getCellFaceAreas(integer(IN) : axis,
!!                             integer(IN) : edge,
!!                             integer(IN) : level,
!!                             integer(IN) : lo(1:MDIM),
!!                             integer(IN) : hi(1:MDIM),
!!                             real(OUT)   : areas) 
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!   axis - specifies the face-centered index space 
!!   level - 
!!   lo/hi - the lower-left and upper-right points in the global face-centered
!!           index space that specify the region of faces whose areas are 
!!           requested
!!   areas - the array in which the requested area values will be stored
!!
!!***

#include "constants.h"

subroutine Grid_getCellFaceAreas(axis, level, lo, hi, areas)
   integer, intent(IN)  :: axis
   integer, intent(IN)  :: level
   integer, intent(IN)  :: lo(1:MDIM)
   integer, intent(IN)  :: hi(1:MDIM)
   real,    intent(OUT) :: areas(:, :, :)

   areas(:,:,:) = 0.0
end subroutine Grid_getCellFaceAreas


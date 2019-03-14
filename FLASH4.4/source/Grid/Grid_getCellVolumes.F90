!!****f* source/Grid/Grid_getCellVolumes
!!
!! NAME
!!  Grid_getCellVolumes
!!
!! SYNOPSIS
!!  call Grid_getCellVolumes(integer(IN) : level,
!!                           integer(IN) : lo(1:MDIM),
!!                           integer(IN) : hi(1:MDIM),
!!                           real(OUT)   : volumes) 
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!   level - 
!!   lo/hi - the lower-left and upper-right points in the global cell-centered
!!           index space that specify the region of cells whose volumes are
!!           requested
!!   volumes - the array in which the requested volume values will be stored
!!
!!***

#include "constants.h"

subroutine Grid_getCellVolumes(level, lo, hi, volumes)
   integer, intent(IN)  :: level
   integer, intent(IN)  :: lo(1:MDIM)
   integer, intent(IN)  :: hi(1:MDIM)
   real,    intent(OUT) :: volumes(:, :, :)

   volumes(:,:,:) = 0.0
end subroutine Grid_getCellVolumes


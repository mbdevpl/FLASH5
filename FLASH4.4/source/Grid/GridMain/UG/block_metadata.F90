!!****ih* source/Grid/GridMain/AMR/paramesh/block_metadata
!!
!!
!!
!!****

#include "constants.h"
#include "Flash.h"

module block_metadata

    implicit none

    private

    type, public :: block_metadata_t
        integer :: id
        integer :: cid(MDIM)
        integer :: stride(MDIM)
        integer :: level
        integer :: limits(LOW:HIGH, MDIM)
        integer :: limitsGC(LOW:HIGH, MDIM)
        integer :: localLimits(LOW:HIGH, MDIM)
        integer :: localLimitsGC(LOW:HIGH, MDIM)
    end type block_metadata_t

end module block_metadata


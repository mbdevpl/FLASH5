!!****ih* source/Grid/GridMain/AMR/Amrex/block_metadata
!!
!!
!!
!!****

module block_metadata

    implicit none

#include "constants.h"

    private

    type, public :: block_metadata_t
        integer :: level
        integer :: grid_index
        integer :: limits(LOW:HIGH, MDIM)
        integer :: limitsGC(LOW:HIGH, MDIM)
    end type block_metadata_t

end module block_metadata


!!****ih* source/Grid/GridMain/AMR/Amrex/block_metadata
!!
!!
!!
!!****

module block_metadata

    use amrex_box_module, ONLY : amrex_box

    implicit none

#include "constants.h"

    private

    type, public :: block_metadata_t
        integer         :: level
        integer         :: grid_index
        type(amrex_box) :: box
    end type block_metadata_t

end module block_metadata


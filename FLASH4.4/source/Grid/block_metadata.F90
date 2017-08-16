!!****ih* source/Grid/GridMain/AMR/paramesh/block_metadata
!!
!!
!!
!!****

module block_metadata

    implicit none

#include "constants.h"

    private

    type, public :: block_metadata_t
        integer :: id
        integer :: cid(MDIM)
        integer :: stride(MDIM)
        integer :: level
        integer :: limits(LOW:HIGH, MDIM)
        integer :: limitsGC(LOW:HIGH, MDIM)
        real    :: del(MDIM)
    end type block_metadata_t

end module block_metadata


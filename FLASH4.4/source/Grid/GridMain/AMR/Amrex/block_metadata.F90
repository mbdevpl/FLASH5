!!****ih* source/Grid/GridMain/AMR/Amrex/block_metadata
!!
!!
!!
!!****

module block_metadata

    implicit none

#include "constants.h"

    private

    public :: bmd_print

    type, public :: block_metadata_t
        integer :: level
        integer :: grid_index
        integer :: limits(LOW:HIGH, MDIM)
        integer :: limitsGC(LOW:HIGH, MDIM)
        integer :: localLimits(LOW:HIGH, MDIM)
        integer :: localLimitsGC(LOW:HIGH, MDIM)
    end type block_metadata_t

contains

    subroutine bmd_print(block)
        type(block_metadata_t), intent(IN) :: block

        write(*,*) "Level      = ", block%level
        write(*,*) "Grid Index = ", block%grid_index
        write(*,*) "Limits     = ", block%limits(LOW, :)
        write(*,*) "             ", block%limits(HIGH, :)
        write(*,*) "LimitsGC   = ", block%limitsGC(LOW, :)
        write(*,*) "             ", block%limitsGC(HIGH, :)
    end subroutine bmd_print

end module block_metadata


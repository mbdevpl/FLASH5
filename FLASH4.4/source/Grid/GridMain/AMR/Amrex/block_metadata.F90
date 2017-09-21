!!****ih* source/Grid/GridMain/AMR/Amrex/block_metadata
!!
!!  NOTES
!!    The data contained in the fields of this derived type are
!!
!!    level - the refinement level at which the block resides.  Levels are
!!            referenced with the index set {1, 2, ..., gr_lRefineMax}, where
!!            1 refers to the coarsest level and mesh refinement increases with
!!            increasing level index.  The value of gr_lRefineMax is fixed to
!!            the value of the runtime parameter lrefine_max.
!!
!!    grid_index - used internally by AMReX.  Client code should never make
!!                 direct use of the value stored in this field.
!!
!!    limits - a two-dimensional array where limits(LOW, :) is the array
!!             containing the (i, j, k) coordinate in index space of the 
!!             lower-leftmost corner cell of the block.  Similarly, 
!!             limits(HIGH, :) is the array containing the (i, j, k) 
!!             coordinate of the upper-rightmost corner cell of the block.
!!             The two coordinates are specified for the cell-center index space
!!             defined so that 
!!                - i is an integer in {1, 2, ..., n_x_cells},
!!                - j is an integer in {1, 2, ..., n_y_cells, and
!!                - k is an integer in {1, 2, ..., n_z_cells},
!!             where for FLASH, 
!!                  n_[xyz]_cells = 2^(level - 1) * N[XYZ]B * nblock[xyz].
!!             This index space is the global index space so that (1, 1, 1) is
!!             the lower-leftmost cell in the full domain of the problem.
!!
!!    limitsGC - a two-dimensional array identical to limits except that the
!!               two cell coordinates in index space refer to the lower-leftmost
!!               and the upper-rightmost guardcells for the block.  The
!!               coordinates are specified w.r.t. an expanded version of the
!!               global cell-center index space used for limits so that
!!                - i is an integer in {1-gr_iguard, ..., n_x_cells+gr_iguard},
!!                - j is an integer in {1-gr_jguard, ..., n_y_cells+gr_jguard},
!!                - k is an integer in {1-gr_kguard, ..., n_z_cells+gr_kguard}.
!!               Hence, (1, 1, 1) still refers to the lower-leftmost guardcell for
!!               the full domain of the problem.  If gr_[ijk]guard=4, then 
!!               (-3, -3, -3) indexes the lower-leftmost guardcell in the
!!               full domain.
!!
!!    localLimits/localLimitsGC - TODO: Fill these is with gory detail.
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


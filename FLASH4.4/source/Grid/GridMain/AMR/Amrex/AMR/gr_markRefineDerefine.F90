!!****if* source/Grid/GridMain/Amrex/AMR/gr_markRefineDerefine
!!
!! NAME
!!  gr_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  gr_markRefineDerefine(integer(IN) :: iref,
!!                        real(IN) :: refine_cutoff,
!!                        real(IN) :: derefine_cutoff,
!!                        real(IN) :: refine_filter)
!!  
!!  DESCRIPTION
!!  
!!    Blocks are marked for refining or derefining.
!!    This version uses the second derivative calculations on the specified variable to 
!!    determine if the block needs more resoultion (refine) or less resolution (derefine)
!!    de/refine_cutoff are the thresholds for triggering the corresponding action.
!!    Once the blocks have been marked, the control is passed to Paramesh to update refinement.
!!
!!  ARGUMENTS 
!!
!!    iref - index of the refinement variable in data structure "unk"
!!
!!    refine_cutoff - the threshold value for triggering refinement 
!!
!!    derefine_cutoff - the threshold for triggereing derefinement
!!
!!    refine_filter - makes sure that error calculations to determine refinement
!!                    don't diverge numerically 
!! 
!!  NOTES
!!  
!!    See Grid_markRefineDerefine
!!
!!  SEE ALSO
!!  
!!    Grid_markRefineDerefine
!!
!!***

#include "constants.h"
#include "flash_bool.h"
#include "Flash.h"

!Just in case.  There should be no reference to these defines!
#undef NXB
#undef NYB
#undef NZB

subroutine gr_markRefineDerefine(iref, refine_cutoff, &
                                 derefine_cutoff, refine_filter)
  
  implicit none
  include "Flash_mpi.h"

  integer, intent(IN) :: iref
  real,    intent(IN) :: refine_cutoff, derefine_cutoff, refine_filter

  ! DEVNOTE: Implement for AMReX

end subroutine gr_markRefineDerefine


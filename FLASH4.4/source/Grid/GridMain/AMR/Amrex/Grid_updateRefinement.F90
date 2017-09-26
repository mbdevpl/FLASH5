!!****if* source/Grid/GridMain/paramesh/Grid_updateRefinement
!!
!! NAME
!!
!!  Grid_updateRefinement
!!
!!
!! SYNOPSIS
!!  
!!  call Grid_updateRefinement(integer(IN) :: nstep,
!!                             real(IN)    :: time,
!!                    OPTIONAL,logical(OUT):: gridChanged)
!!
!!
!! DESCRIPTION
!!
!!  Apply the user-defined refinment critera to determine which blocks need 
!!  to be refined and derefined.  Once the blocks are marked, call 
!!  gr_updateGridRefinement to carry out the rest of the routine.
!!  The internal routine does the refinements (amr_refine_derefine)
!!  During this
!!  stage, the blocks are redistributed across processors (if needed).  
!!  After the refinement, the newly created child blocks are filled via
!!  prolongation from the coarse parents.  This prolongation step can use
!!  prolongation routines from paramesh, or defined by the user
!!  Once the prolongation is done, the guardcells are filled.  Finally, the
!!  EOS is called on the block interiors to make them thermodynamically
!!  consistent. The internal routine also calls Particles_updateRefinement to
!!  move the particles to the correct block after the grid refines.
!!
!!
!!
!! ARGUMENTS
!!
!!  nstep : current step number
!!  time  : current evolution time
!!  gridChanged : returns TRUE if grid may actually have changed.
!!
!!***

subroutine Grid_updateRefinement(nstep, time, gridChanged)
  use amrex_amrcore_module, ONLY : amrex_regrid

  use Grid_interface,       ONLY : Grid_fillGuardCells
  use Grid_data,            ONLY : gr_nrefs, gr_maxRefine
  use Timers_interface,     ONLY : Timers_start, Timers_stop
 
  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in)            :: nstep
  real,    intent(in)            :: time
  logical, intent(out), OPTIONAL :: gridChanged

  ! We only consider refinements every nrefs timesteps.
  if (mod(nstep, gr_nrefs) == 0) then
     call Timers_start("Grid_updateRefinement")

     ! AMReX uses 0-based level index set
     call Grid_fillGuardCells(CENTER, ALLDIR)
     write(*,'(A,I4,A,E9.3)') "[Grid_updateRefinement] AMReX Regridding @ step=", & 
                              nstep, " / time = ", time
     call amrex_regrid(0, time)
 
     call Timers_stop("Grid_updateRefinement")

     if (present(gridChanged)) then
        ! DEV: FIXME: Shouldn't this actually check if AMReX
        ! decided to regrid based on the contents of the physical quantities?
        gridChanged = .TRUE.
     end if
  else
     if (present(gridChanged)) then
        gridChanged = .FALSE.
     end if
  end if

end subroutine Grid_updateRefinement


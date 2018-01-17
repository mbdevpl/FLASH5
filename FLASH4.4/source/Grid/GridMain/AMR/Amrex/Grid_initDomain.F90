!!****if* source/Grid/GridMain/AMR/Amrex/Grid_initDomain
!!
!! NAME
!!
!!  Grid_initDomain
!!
!! SYNOPSIS
!!
!!  Grid_initDomain(logical(IN)  :: restart,
!!                  logical(INOUT) :: particlesInitialized)
!!
!! DESCRIPTION
!!
!!  Create the mesh, initialize all the mesh data structures
!!  and apply initial conditions
!!
!!  Initially very few blocks are created (number supplied at runtime).
!!  then user-defined refinment critera is applied to determine the 
!!  blocks that need to be refined and derefined.  
!!
!!  After the refinement, the newly created child blocks are filled via
!!  prolongation from the coarse parents.  This prolongation step can use
!!  prolongation routine supplied with paramesh or defined by the user.
!!
!!  Once the prolongation is done, the guardcells are filled.  Finally, the
!!  EOS is called on the block interiors to make them thermodynamically
!!  consistent.
!!
!!  In simulations with particles, under certain conditions particle
!!  positions will also be initialized.  Currently this is the case
!!  if and only if the runtime parameter refine_on_particle_count is
!!  true.
!!
!! ARGUMENTS
!!
!!  restart : is true if the execution is starting from a checkpoint
!!            file, otherwise false.
!!  particlesInitialized : is true if particle positions were initialized before returning
!!                         from this routine
!!
!! NOTES
!!  When restarting from a checkpoint file, block interiors are assumed to
!!  have been filled when this interface is called. The EOS is not called on
!!  the block interiors in this implementation for use with Paramesh. It is
!!  assumed that data is already thermodynamically consistent, because
!!  that is how data are written to checkpoint files.
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Grid_initDomain(restart,particlesInitialized)
  use amrex_fort_module,    ONLY : wp => amrex_real
  use amrex_amr_module,     ONLY : amrex_init_from_scratch, &
                                   amrex_max_level

  use gr_physicalMultifabs, ONLY : unk, &
                                   facevarx, facevary, facevarz

  use Grid_data, ONLY : gr_eosMode, gr_eosModeNow

  implicit none

  logical, intent(IN)    :: restart
  logical, intent(INOUT) :: particlesInitialized

  real(wp), parameter :: T_INIT = 0.0_wp

  ! DEV: TODO Implement parameters

  !!!!!----- ALLOCATE DATA STRUCTURES
  ! multifabs 
  !
  ! NOTE: We implement these with the 0-based level indexing scheme native to
  ! AMReX instead of the 1-based level indexing scheme of FLASH.
  !   => all code dealing with multifabs arrays must consider the need for 
  !      index translation
  allocate(unk     (0:amrex_max_level))
  allocate(facevarx(0:amrex_max_level))
  allocate(facevary(0:amrex_max_level))
  allocate(facevarz(0:amrex_max_level))

  ! Setup grids and initialize the data
  call amrex_init_from_scratch(T_INIT)

  gr_eosModeNow = gr_eosMode !may be different from gr_eosModeInit
end subroutine Grid_initDomain


!!****if* source/Grid/GridMain/AMR/Amrex/Grid_initDomain
!!
!! NAME
!!
!!  Grid_initDomain
!!
!! SYNOPSIS
!!
!!  Grid_initDomain(logical(IN)    :: restart,
!!                  logical(INOUT) :: particlesInitialized)
!!
!! DESCRIPTION
!!
!!  Create the coarsest mesh, initialize all the mesh data structures,
!!  apply initial conditions, and run EoS on interior and guardcells to make
!!  them thermodynamically consistent.
!!
!!  User-defined refinement critera is applied to determine the 
!!  blocks that require refinement.  All new child blocks are filled
!!  with the initial conditions and EoS is run on interior and guardcells.
!!
!!  Please see the documentation for gr_initNewLevelCallback for more
!!  information regarding how the EoS runs are done.
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
  use Driver_interface,     ONLY : Driver_abortFlash

  implicit none

  logical, intent(IN)    :: restart
  logical, intent(INOUT) :: particlesInitialized

  real(wp), parameter :: T_INIT = 0.0_wp

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

  ! DEV: TODO Implement parameters
  if (.NOT. restart) then
    !  This creates all refinement levels needed based on ICs,
    !  runs EoS on interiors, fills GCs, and runs EoS on GCs.
    !  All this is done through the callback functions.
    call amrex_init_from_scratch(T_INIT)
  else 
    call Driver_abortFlash("[Grid_initDomain] restarts not yet implemented")
  end if

  particlesInitialized = .FALSE.
end subroutine Grid_initDomain


!!****if* source/Particles/ParticlesMain/Particles_finalize
!!
!! NAME
!!    Particles_finalize
!!
!! SYNOPSIS
!!    Particles_finalize( )
!!
!! DESCRIPTION
!!
!!    Finalize routine for the particle unit.  Removes memory usage
!!      set in Particles_init.
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!
!!***

#include "Flash.h"

subroutine Particles_finalize()
  use pt_interface, ONLY : pt_initFinalize
  use Particles_data, ONLY : particles,useParticles

#ifdef FLASH_GRID_AMREX
  use Particles_data, ONLY : pt_containers
  use amrex_particlecontainer_module, ONLY: amrex_particlecontainer_destroy
#endif
  implicit none


  integer :: i
  if(useParticles) deallocate(particles)
#ifdef FLASH_GRID_AMREX
  if(useParticles) then
    do i=1,NPART_TYPES
       call amrex_particlecontainer_destroy(pt_containers(i))
    end do
  end if
#endif
  call pt_initFinalize()
end subroutine Particles_finalize

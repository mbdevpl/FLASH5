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
  
subroutine Particles_finalize()
  use Particles_data, ONLY : particles,useParticles, pt_containers
  use pt_interface, ONLY : pt_initFinalize
use amrex_particlecontainer_module, ONLY: amrex_particlecontainer_destroy
  implicit none

#include "Flash.h"

  integer :: i
  if(useParticles) deallocate(particles)
  if(useParticles) then
    do i=1,NPART_TYPES
       call amrex_particlecontainer_destroy(pt_containers(i))
    end do
  end if
  call pt_initFinalize()
end subroutine Particles_finalize

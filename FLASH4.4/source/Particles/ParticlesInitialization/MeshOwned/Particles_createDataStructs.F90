#include "Flash.h"

subroutine Particles_createDataStructs
  use Particles_data, ONLY :  pt_containers, useParticles
  use amrex_amr_module, ONLY : amrex_get_amrcore
  use amrex_particlecontainer_module, ONLY : amrex_particlecontainer_build

  implicit none
    integer :: i
    if(useParticles) then
        do i=1, NPART_TYPES
        call amrex_particlecontainer_build(pt_containers(i), amrex_get_amrcore())
        end do
    endif

end subroutine Particles_createDataStructs


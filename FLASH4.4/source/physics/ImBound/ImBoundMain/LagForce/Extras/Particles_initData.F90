!!****if* source/physics/ImBound/ImBoundMain/LagForce/Extras/Particles_initData
!!
!! NAME
!!
!!  Particles_initData
!!
!! SYNOPSIS
!!
!!  call Particles_initData(logical(IN) :: restart)
!!
!! DESCRIPTION
!!
!!  Initialization of particle attributes after they have been
!!  dropped into their positions.
!!
!!  This interface is called during Flash initialization after
!!   o  particles positions,
!!   o  mass properties for active particles,  (NOTE: it is true in
!!                                              most cases, but not a 
!!                                              requirement)
!!   o  and particle tags
!!  have been initialized (or possibly read from a checkpoint, if
!!  restart is true).
!!
!!  This is where velocity components are initialized.
!!  For passive particles, this is done by mapping from
!!  the fluid velocity field. The sequence of initialization
!!  is such that fluid fields are guaranteed to have been
!!  initialized to their final starting values, and the
!!  configuration of the grid has been finalized including
!!  creation of blocks up to the maximum refinement level,
!!  when the Particles_initData interface is called.
!!
!!  The default implementation in the Particles main subunit
!!  calls Simulation_initParticleAttrib before returning, so
!!  if custom initialization of additional particle properties
!!  is needed, that's the place to implement it.
!!
!! ARGUMENTS
!!
!!  restart - true if restarting from a checkpoint, false otherwise.
!!
!!
!! SEE ALSO
!!
!!  Driver_intFlash
!!  Simulation_initParticleAttrib
!!
!!***

subroutine Particles_initData(restart, partPosInitialized)
  use Simulation_interface,ONLY : Simulation_initParticleAttrib
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,ONLY : Grid_mapMeshToParticles, Grid_sortParticles
  use Particles_data,ONLY : particles, pt_numLocal, pt_maxPerProc, &
       pt_posAttrib, pt_velNumAttrib, pt_velAttrib, pt_velInitialized,&
       useParticles, pt_meshMe, pt_meshNumProcs,pt_typeInfo,pt_numLocal
  use Particles_interface, ONLY : Particles_initPositions
  use pt_interface, ONLY : pt_updateTypeDS
  implicit none
#include "Flash.h"
#include "constants.h"
#include "Particles.h"
  logical, intent(IN) :: restart
  logical, intent(INOUT) :: partPosInitialized
  integer :: part_props=NPART_PROPS
  logical :: updateRefine, needVel
  integer :: i,p_begin,p_count,p_end
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk

  ! Return immediately if useParticles is false.
  if (.NOT. useParticles) then
     return
  end if

  updateRefine = .false.

  pt_velInitialized = .TRUE.

  ! Call Simulation_initParticleAttrib, normally just a stub, to allow setups some customization
  call Simulation_initParticleAttrib(restart)


end subroutine Particles_initData

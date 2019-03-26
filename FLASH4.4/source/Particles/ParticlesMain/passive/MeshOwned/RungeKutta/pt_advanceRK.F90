!!****if* source/Particles/ParticlesMain/Amrex/passive/RungeKutta/pt_advanceRK_pc
!!
!! NAME
!!
!!  pt_advanceRK_pc
!!
!! SYNOPSIS
!!
!!  call pt_advanceRK(real(in)   :: dtOld,
!!                         real(in)   :: dtNew,
!!                         integer(in):: ind)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  
!!  This version implements the improved Euler method, also called Heun's method,
!!  for integration in time. The improved Euler methos is one of the family of
!!  2-stage, second-order Runge-Kutta methods.  (It can probably also be framed
!!  as a simple Predictor-Corrector method, with a first-order predictor and
!!  corrector.)
!!
!!  In detail:
!!      x*(t+dtNew) = x(t)  + dtNew *        v(x(t),t)
!!      x(t+dtNew)  = x(t)  + dtNew * 1/2* [ v(x(t),t) + v(x*(t+dtNew),t+dtNew) ]
!!      v(t+dtNew)  =  v(x(t+dtNew),t+dtNew)
!!  where x* is ephemeral (a predicted position that then is corrected after
!!  another evaluation there).
!!
!!  Implementation detail:
!!  This can be rewritten to save on memory for intermediate result storage:
!!      x*(t+dtNew) = x (t)        + dtNew *          v(x(t),t)
!!      x(t+dtNew)  = x*(t+dtNew)  + dtNew * 1/2* [ - v(x(t),t) + v(x*(t+dtNew),t+dtNew) ]
!!      v(t+dtNew)  =  v(x(t+dtNew),t+dtNew)
!!  where x* can be stored in the same location as the previous and final x.
!!
!! ARGUMENTS
!!
!!   dtOld -- not used in this first-order scheme
!!   dtNew -- current time increment
!!   ind -- index into pt_typeInfo and into pt_containers[] array for this type of particles
!!
!! SIDE EFFECTS
!!
!!  Updates the POS{X,Y,Z} and VEL{X,Y,Z} properties of particles in the particles structure.
!!  Sorts particles in the particles structure by calling Grid_sortParticles.
!!
!! NOTES
!!
!!  No special handling is done for the first call - it is assumed that particle
!!  initialization fills in initial velocity components properly.
!!
!!***

!===============================================================================

subroutine pt_advanceRK(dtOld,dtNew, p_beg, p_end, ind)
    
  use Particles_data, ONLY: pt_numLocal, pt_maxPerProc, &
       useParticles, pt_typeInfo, &
       pt_gcMaskForAdvance, pt_gcMaskSizeForAdvance, pt_meshMe, &
       pt_posAttrib, pt_velNumAttrib,pt_velAttrib
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator, Grid_mapMeshToParticles
  use Particles_data, ONLY : pt_containers
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,        ONLY : Grid_tile_t
  use amrex_particlecontainer_module, ONLY : amrex_particle, amrex_particlecontainer,&
        amrex_particlecontainer_build, amrex_particlecontainer_destroy
  use amrex_amr_module, ONLY : amrex_get_amrcore
  
  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Particles.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  integer, INTENT(in) :: ind, p_beg, p_end

  integer       :: i,particleTypes, p_count
  
  integer,dimension(MAXBLOCKS, 1) :: perBlk
  real          :: jumpx,jumpy,jumpz
  real, dimension(NDIM)          :: jump
  real,allocatable :: origVel(:,:)
  integer :: part_props=NPART_PROPS

  integer :: mapType 
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)    :: tileDesc
  type(amrex_particle), pointer :: particles(:)
  type(amrex_particle) :: oneParticle
  type(amrex_particle), pointer :: particlesSaveVelocity(:)
  type(amrex_particlecontainer) :: pcSaveVelocity
  integer :: j
!!------------------------------------------------------------------------------
  
  ! Don't do anything if runtime parameter isn't set
  if (.not.useParticles ) return

  mapType=pt_typeInfo(PART_MAPMETHOD,ind)

  call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
  do while(itor%isValid())
    call itor%currentTile(tileDesc)
    particles => pt_containers(ind)%get_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index)
    p_count = size(particles)

  ! Update the particle positions to temporary ("predicted") values
    do i = 1, p_count
        do j=1,NDIM
            jump(j) = dtNew * particles(i)%vel(j)
            particles(i)%pos(j) = particles(i)%pos(j) + jump(j)
        end do
    enddo
    call itor%next()
  enddo             ! leaf itor enddo
  call Grid_releaseTileIterator(itor)


  ! Now save the original velocity values in a temporary amrex_particlecontainer
  
  call amrex_particlecontainer_build(pcSaveVelocity, amrex_get_amrcore())
    call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
  do while(itor%isValid())
    call itor%currentTile(tileDesc)
    particles => pt_containers(ind)%get_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index)
    p_count = size(particles)
  ! Copy values of particle velocities
    do i = 1, p_count
        oneParticle = particles(i)
        call pcSaveVelocity%add_particle(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index, oneParticle)
    enddo
    call itor%next()
  enddo             ! leaf itor enddo
  call Grid_releaseTileIterator(itor)


  ! Map the updated gas velocity field at the temporary positions to
  ! obtain a second estimate of velocities;

  call Grid_mapMeshToParticles(ind,&
       part_props, BLK_PART_PROP, &
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

  ! Adjust particle positions, using the second point velocities
  call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
  do while(itor%isValid())
    call itor%currentTile(tileDesc)
    particles => pt_containers(ind)%get_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index)
    particlesSaveVelocity => pcSaveVelocity%get_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index)
    p_count = size(particles)
    do i = 1, p_count
        do j=1,NDIM
            particles(i)%pos(j) =  particles(i)%pos(j) + &
                dtNew * 0.5*(particles(i)%vel(j) - particlesSaveVelocity(i)%vel(j) )
        end do
    enddo
    call itor%next()
  enddo             ! leaf itor enddo
  call Grid_releaseTileIterator(itor)

  ! done with this temporary amrex_particlecontainer
  call amrex_particlecontainer_destroy(pcSaveVelocity)

  ! Map the updated gas velocity field onto the current particle positions to
  ! obtain the updated particle velocities - for the next integration step
  ! as well as for particle plot files etc.

  call Grid_mapMeshToParticles(ind,&
       part_props, BLK_PART_PROP,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)
  
  return
!!------------------------------------------------------------------------------
  
end subroutine pt_advanceRK



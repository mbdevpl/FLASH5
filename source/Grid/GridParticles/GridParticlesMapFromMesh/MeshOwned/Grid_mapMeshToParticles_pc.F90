!!****f* source/Grid/GridParticles/GridParticlesMapFromMesh/Amrex
!!
!! NAME
!!
!!  Grid_mapMeshToParticles
!!
!! SYNOPSIS
!!
!!  Grid_mapMeshToParticles(int(in) :: ptContainerPos,
!!                          integer(in) :: part_props,
!!                          integer(in) :: part_blkID,
!!                          integer(in) :: posAttrib(POS),
!!                          integer(in) :: numAttrib,
!!                          integer(in) :: attrib(:,numAttrib),
!!                          integer(in) :: mapType,
!!               optional,  integer(in) :: gridDataStruct)
!!
!! DESCRIPTION
!!
!!  Routine to map a quantity defined on the mesh onto the particle positions.
!!
!! ARGUMENTS
!!
!!     ptContainerPos:     Position of particles type in AMReX's ParticleContainer
!!                                    array pt_containers of Particle_data module
!!     part_props : number of particle attributes
!!     part_blkID : the index of particle attributes that carries the block number (only for paramesh)
!!     posAttrib           : particles data structure indices pointing
!!                           to particle positions
!!     numAttrib           : number of attributes that need to be mapped
!!     attrib              : list of attributes and their corresponding
!!                           mesh data structure indices
!!                           processor
!!     mapType : method for mapping grid quantities to particles
!!     gridDataStruct: The Grid data structure that varGrid refers to; one of
!!                    CENTER, FACEX, FACEY, FACEZ, GRIDVAR.  
!!                    If this argument is not present, the default is CENTER, 
!!                    so that attrib indicates a list of variable in UNK.
!! NOTES: 
!! Arguments part_blkID is dummy for amrex implementation
!!  
!!***
!*******************************************************************************

#include "Flash.h"
#include "constants.h"
#include "Particles.h"

subroutine Grid_mapMeshToParticles_pc (ptContainerPos, part_props,part_blkID,&
                                    posAttrib,&
                                    numAttrib, attrib,&
                                    mapType,gridDataStruct)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use Particles_interface, ONLY : Particles_mapFromMesh
  use Particles_data, ONLY : pt_containers

  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,        ONLY : Grid_tile_t
  use amrex_particlecontainer_module, ONLY : amrex_particle


  implicit none

  integer, INTENT(in) :: part_props, numAttrib, part_blkID
  integer, INTENT(IN) :: ptContainerPos
  integer,dimension(MDIM), intent(IN) :: posAttrib
  integer, dimension(PART_ATTR_DS_SIZE,numAttrib),INTENT(in) :: attrib
  integer, INTENT(IN) :: mapType
  integer, optional, intent(IN) :: gridDataStruct
  !-----------------------------------------------------------------------------------------------!
  integer :: gDataStruct
  real, contiguous, pointer :: solnData(:,:,:,:)
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)    :: tileDesc
  type(amrex_particle), pointer :: particles(:)
  integer :: numParticlesOnBlock, i, j
  real, dimension(LOW:HIGH,MDIM) :: bndBox
  real,dimension(MDIM) :: delta, pos
  real, dimension(numAttrib) :: partAttribVec

  nullify(solnData)
  if(present(gridDataStruct)) then
     gDataStruct=gridDataStruct
  else
     gDataStruct=CENTER
  end if
!! Using grid iterator will return all blocks/tiles.  Even the ones that do not have 
!! particles on this container. Better way is to use the "ParIter'' which is not
!! yet exposed from AMReX in Fortran. There is now wrapper in FLASH for the same.
  call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
     do while(itor%isValid())
        call itor%currentTile(tileDesc)
        call tileDesc%getDataPtr(solnData, CENTER)
!!TODO:: Particles_getFromTile(pt_containers(ptContainerPos), tileDesc, particles)   !!Correct usage once a wrapper get function is introduced
        particles => pt_containers(ptContainerPos)%get_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index)
        numParticlesOnBlock = size(particles)
            if(numParticlesOnBlock>0) then
                call tileDesc%boundBox(bndBox)
                call tileDesc%deltas(delta)
                do i = 1, numParticlesOnBlock
                    do j = 1,NDIM
                        pos(j) = particles(i)%pos(j)
                    end do
                    call Particles_mapFromMesh (mapType, numAttrib, attrib,&
                    pos, bndBox,delta, tileDesc%limits, solnData, partAttribVec)
!                   Assign values to particles(i)%vel from output partAttribVec
!                   Assuming that velocities are in partAttribVec are in indices 1 to 3
!                   To change this, keep the old particle array and then 1) assign partAttribVec(j)
!                   to particles(attrib(PART_DS_IND,j),i) ; 2) particles(i)%vel(j) = particle(pt_velAttrib(PART_DS_IND?? , j), i)
                    do j = 1,NDIM
                        particles(i)%vel(j) = partAttribVec(j)
                    end do
                end do
            end if
        call tileDesc%releaseDataPtr(solnData, CENTER)
        nullify(solnData)
        call itor%next()
     enddo
  call Grid_releaseTileIterator(itor)

  return
end subroutine Grid_mapMeshToParticles_pc
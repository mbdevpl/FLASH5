!!****f* source/Grid/Amrex
!!
!! NAME
!!
!!  Grid_mapMeshToParticles
!!
!! SYNOPSIS
!!
!!  Grid_mapMeshToParticles(int(in) :: pt_containerPos,
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
!!     pt_containerPos:     Position of particles type in AMReX's ParticleContainer
!!                                    array pt_containers of Particle_data module
!!     part_props : number of particle attributes
!!     part_blkID : the index of particle attributes that carries the block number
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
!!
!!  
!!***
!*******************************************************************************

subroutine Grid_mapMeshToParticles_pc (pt_containerPos, part_props,part_blkID,&
                                    posAttrib,&
                                    numAttrib, attrib,&
                                    mapType,gridDataStruct)

#include "Flash.h"
  implicit none

  integer, INTENT(in) :: part_props, numAttrib, part_blkID
  integer, INTENT(IN) :: pt_containerPos
  integer,dimension(MDIM), intent(IN) :: posAttrib
  integer, dimension(PART_ATTR_DS_SIZE,numAttrib),INTENT(in) :: attrib
  integer, INTENT(IN) :: mapType
  integer, optional, intent(IN) :: gridDataStruct

  return
end subroutine Grid_mapMeshToParticles_pc

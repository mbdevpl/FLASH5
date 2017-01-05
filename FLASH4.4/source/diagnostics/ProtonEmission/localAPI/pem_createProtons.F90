!!****if* source/diagnostics/ProtonEmission/localAPI/pem_createProtons
!!
!! NAME
!!
!!  pem_createProtons
!!
!! SYNOPSIS
!!
!!  call pem_createProtons (integer, intent (in)    :: blockCount, 
!!                          integer, intent (in)    :: blockList (:),
!!                          real,    intent (in)    :: timeStep,
!!                          logical, intent (in)    :: countOnly,
!!                          logical, intent (inout) :: emissionInitiate,
!!                          logical, intent (inout) :: emissionIncomplete)
!!
!! DESCRIPTION
!!
!!  Generates emission protons inside the domain and places them in their initial blocks.
!!  This routine calls the appropriate subroutines according to the domain grid geometry
!!  specified. Currently only 3D cartesian geometries are allowed.
!!
!! ARGUMENTS
!!
!!  blockCount         : Number of blocks on current processor
!!  blockList          : All block ID numbers
!!  timeStep           : The current time step duration
!!  countOnly          : Should the protons only be counted and not created?
!!  emissionInitiate   : Is true, if the emission process starts
!!  emissionIncomplete : If true, the emission process is not done yet
!!
!!***

subroutine pem_createProtons (blockCount,                  &
                              blockList,                   &
                              timeStep,                    &
                              countOnly,                   &
                                         emissionInitiate, &
                                         emissionIncomplete)

  implicit none

  integer, intent (in)    :: blockCount
  integer, intent (in)    :: blockList (1:blockCount)
  real,    intent (in)    :: timeStep
  logical, intent (in)    :: countOnly
  logical, intent (inout) :: emissionInitiate
  logical, intent (inout) :: emissionIncomplete

  emissionInitiate = .false.
  emissionIncomplete = .false.

  return
end subroutine pem_createProtons

!!****if* source/diagnostics/ProtonEmission/localAPI/pem_createProtons3DRec
!!
!! NAME
!!
!!  pem_createProtons3DRec
!!
!! SYNOPSIS
!!
!!  call pem_createProtons3DRec (integer, intent (in)    :: blockCount,
!!                               integer, intent (in)    :: blockList (:),
!!                               real,    intent (in)    :: timeStep,
!!                               logical, intent (in)    :: countOnly,
!!                               logical, intent (inout) :: emissionInitiate,
!!                               logical, intent (inout) :: emissionIncomplete)
!!
!! DESCRIPTION
!!
!!  Generates emission protons inside the cells of each 3D cartesian block on the current processor.
!!  The location of the protons is at the corresponding cell centers with their velocity arrows
!!  pointing radially outward.
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
!! NOTES
!!
!!***

subroutine pem_createProtons3DRec (blockCount,                  &
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
end subroutine pem_createProtons3DRec

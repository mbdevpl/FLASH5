!!****if* source/diagnostics/ProtonEmission/localAPI/pem_createBlockProtons3DRec
!!
!! NAME
!!
!!  pem_createBlockProtons3DRec
!!
!! SYNOPSIS
!!
!!  call pem_createBlockProtons3DRec (integer, intent (in)    :: blockID,
!!                                    integer, intent (in)    :: iminBlock,
!!                                    integer, intent (in)    :: imaxBlock,
!!                                    integer, intent (in)    :: jminBlock,
!!                                    integer, intent (in)    :: jmaxBlock,
!!                                    integer, intent (in)    :: kminBlock,
!!                                    integer, intent (in)    :: kmaxBlock,
!!                                    real   , intent (in)    :: timeStep,
!!                                    real   , intent (in)    :: cellVolume,
!!                                    logical, intent (in)    :: countOnly,
!!                                    real   , intent (in)    :: blockData (:,:,:,:),
!!                                    logical, intent (inout) :: blockIncomplete)
!!
!! DESCRIPTION
!!
!!  Creates emission protons for each cell inside the block. On exit, a set of emmision
!!  protons has been added to the overall set of emission protons on the current processor.
!!  Each emission proton will have a position (center of a cell) and velocity components
!!  in the direction of emission.
!!
!!  The space of emission is defined as a spherical cone with origin each cell center and
!!  a specific direction position. The spherical cone is further characterized by its
!!  half apex angle and can have a range from 0 degrees to 180 degrees. A half apex angle
!!  of 180 is identical with a complete 3D sphere. In this case the protons are emitted
!!  statistically around the entire 3D space.
!!
!! ARGUMENTS
!!
!!  blockID         : The block ID number
!!  iminBlock       : Minimum cell i-index limit defining the interior block
!!  imaxBlock       : Maximum cell i-index limit defining the interior block
!!  jminBlock       : Minimum cell j-index limit defining the interior block
!!  jmaxBlock       : Maximum cell j-index limit defining the interior block
!!  kminBlock       : Minimum cell k-index limit defining the interior block
!!  kmaxBlock       : Maximum cell k-index limit defining the interior block
!!  timeStep        : The current time step duration
!!  cellVolume      : The volume of each cell in the block
!!  countOnly       : Should the protons only be counted and not created?
!!  blockData       : Four-dimensional array containing the block data
!!  blockIncomplete : If true, some cells of the block have not been processed
!!
!!***

subroutine pem_createBlockProtons3DRec (blockID,              &
                                        iminBlock, imaxBlock, &
                                        jminBlock, jmaxBlock, &
                                        kminBlock, kmaxBlock, &
                                        timeStep,             &
                                        cellVolume,           &
                                        countOnly,            &
                                        blockData,            &
                                        blockIncomplete       )

  implicit none

  integer, intent (in)    :: blockID
  integer, intent (in)    :: iminBlock, imaxBlock
  integer, intent (in)    :: jminBlock, jmaxBlock
  integer, intent (in)    :: kminBlock, kmaxBlock
  real,    intent (in)    :: timeStep
  real,    intent (in)    :: cellVolume
  logical, intent (in)    :: countOnly
  real,    intent (in)    :: blockData (:,:,:,:)
  logical, intent (inout) :: blockIncomplete

  blockIncomplete = .false.

  return
end subroutine pem_createBlockProtons3DRec

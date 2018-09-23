!!****f* source/Simulation/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
!!                            integer(IN)  :: blockDesc  )
!!
!!
!!
!! DESCRIPTION
!!  This routine applies initial conditions of a specific simulation
!!  to the specified block.
!!
!! 
!! ARGUMENTS
!!
!!  solnData  -        pointer to solution data
!!  blockDesc -        describes the block to initialize
!!
!! PARAMETERS
!!
!!  eosModeInit -     after this routine sets up initial conditions,
!!                    the grid package calls Eos to insure the
!!                    values are thermodynamically consistent.  This
!!                    parameter controls the mode of application of
!!                    the Eos.  Its default is "dens_ie", and it can
!!                    be one of "dens_ie", "dens_pres", "dens_temp".
!!                    Setting this value to dens_temp, for instance,
!!                    would make it possible to leave this routine
!!                    having just initialized density and temperature,
!!                    leaving Eos to calculate the rest.
!!
!! SEE ALSO
!!
!!  Eos_wrapped
!!***

subroutine Simulation_initBlock(solnData,blockDesc)

  use block_metadata, ONLY : block_metadata_t
  
  implicit none
  
  real,dimension(:,:,:,:),pointer :: solnData
  type(block_metadata_t), intent(in) :: blockDesc
  
  return
end subroutine Simulation_initBlock

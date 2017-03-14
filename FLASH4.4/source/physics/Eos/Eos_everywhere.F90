!!****f* source/physics/Eos/Eos_everywhere
!! NAME
!!
!!  Eos_everywhere
!!
!! SYNOPSIS
!!
!!  call Eos_everywhere(  integer(IN) :: mode,
!!               optional,integer(IN) :: gridDataStruct )
!!
!! DESCRIPTION
!!
!! Apply Eos() to all blocks.
!!
!!  ARGUMENTS
!!
!!   mode : determines which variables are used as Eos input.
!!          The valid values are MODE_DENS_EI (where density and internal
!!          energy are inputs), MODE_DENS_PRES (density and pressure as inputs)
!!          MODE_DENS_TEMP (density and temperature are inputs).
!!          These quantities are defined in constants.h, the argument is
!!          forwarded unchanged to the Eos function call.
!!
!!   gridDataStruct : the grid data structure on whose data Eos is to be applied
!!
!!
!!  EXAMPLE
!!
!!  NOTES
!!
!!  SEE ALSO
!!
!!     Eos
!!     Eos_wrapped
!!     Grid_makeVector
!!
!!***


subroutine Eos_everywhere(mode,gridDataStruct)
  implicit none

  integer, intent(in) :: mode
  integer, optional, intent(IN) :: gridDataStruct

  return
end subroutine Eos_everywhere

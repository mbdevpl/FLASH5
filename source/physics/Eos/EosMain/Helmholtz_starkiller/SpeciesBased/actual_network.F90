! An automatically generated "null" network.  This provides the properties
! of a set of non-reacting species.
!
! nspec            -- the number of species
! naux             -- the number of auxiliary variables
!
! aion             -- atomic number
! zion             -- proton number
!
! spec_names       -- the name of the isotope
! short_spec_names -- an abbreviated form of the species name
!
! aux_names        -- the name of the auxiliary variable
! short_aux_names  -- an abbreviated form of the auxiliary variable


module actual_network

  use bl_types

  implicit none

#include "Flash.h"

  integer, parameter :: nspec = NSPECIES
  integer, parameter :: nspec_evolve = NSPECIES
  integer, parameter :: naux =  NMASS_SCALARS

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len=16), save :: aux_names(naux)
  character (len= 5), save :: short_aux_names(naux)

  double precision  , save :: aion(nspec), zion(nspec)

  !$acc declare create(aion, zion)

  integer, parameter :: nrates = 0
  integer, parameter :: num_rate_groups = 0

contains
  
  subroutine actual_network_init

     use Multispecies_interface, only: Multispecies_getPropertyVector

     implicit none

#include "constants.h"
#include "Multispecies.h"

     integer :: i

     do i = 1, nspec
        call Simulation_mapIntToStr(SPECIES_BEGIN+i-1,short_spec_names(i),MAPBLOCK_UNK)
        spec_names(i) = short_spec_names(i)
     end do

     call Multispecies_getPropertyVector(A,aion)
     call Multispecies_getPropertyVector(Z,zion)

     do i = 1, naux
        call Simulation_mapIntToStr(MASS_SCALARS_BEGIN+i-1,short_aux_names(i),MAPBLOCK_UNK)
        aux_names(i) = short_aux_names(i)
     end do

    !$acc update device(aion, zion)

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    ! Nothing to do here.

  end subroutine actual_network_finalize

end module actual_network

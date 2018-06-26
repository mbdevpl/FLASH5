! Function to compute fluid integrals on volume occupied by
! solid bodies.
! Used on Integral form of nodal force computations (Uhlman 2005).


subroutine sm_assemble_FluidVolIntegrals()

  use sm_assemble_interface, only : sm_assemble_FluidVolInt_rigid

  implicit none

  ! Call assemble fluid integrals on all rigid bodies:
  call sm_assemble_FluidVolInt_rigid()

  ! Call assemble fluid integrals on all flexible bodies:
  ! something like call sm_assemble_FluidIntVol_flex()     

  return
end subroutine sm_assemble_FluidVolIntegrals

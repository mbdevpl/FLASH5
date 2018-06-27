!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Testing/DryRun/sm_surf_currentFluidForce
!!
!! NAME
!! 
!!
!! SYNOPSIS
!!    stub function
!!  
!! DESCRIPTION 
!! 
!!
!! ARGUMENTS 
!!
!!***

subroutine sm_surf_currentFluidForce(ibd, ndofs, Hs_pres, Hs_visc)
  implicit none
  integer, intent(in)  :: ibd, ndofs
  real,    intent(out) :: Hs_pres(ndofs), Hs_visc(ndofs)

  return

end subroutine sm_surf_currentFluidForce

!!****if* source/Grid/localAPI/gr_amrexGet/Geom
!!
!!  NAME 
!!
!! gr_amrexGetGeom
!!
!!  SYNOPSIS
!!
!!  call gr_amrexGetGeom()
!!
!!
!!  DESCRIPTION 
!! This routine returns the amrex_geometry associated with the multifab unk
!! 
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  Solver settings used ?? DESCRIPTION??
!!
!!***
subroutine gr_amrexGetGeom(geom)
  use amrex_geometry_module,     ONLY : amrex_geometry
  implicit none
  type(amrex_geometry), intent(inout) :: geom(:)
  return
end subroutine gr_amrexGetGeom

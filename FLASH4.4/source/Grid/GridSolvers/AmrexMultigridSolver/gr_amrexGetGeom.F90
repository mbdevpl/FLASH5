!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexGetGeom
!!
!!  NAME 
!!
!! gr_amrexGetBa
!!
!!  SYNOPSIS
!!
!!  call gr_amrexGetBa()
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
      use amrex_amr_module, ONLY : amrex_geom, amrex_get_finest_level
      use amrex_geometry_module, ONLY : amrex_geometry
      implicit none
      type(amrex_geometry), allocatable, intent(inout) :: geom(:)
      geom = amrex_geom(0:amrex_get_finest_level())
  return
end subroutine gr_amrexGetGeom

!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexGetDm
!!
!!  NAME 
!!
!! gr_amrexGetDm
!!
!!  SYNOPSIS
!!
!!  call gr_amrexGetDm()
!!
!!
!!  DESCRIPTION 
!! This routine returns the amrex_distromap associated with the multifab unk
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

subroutine gr_amrexGetDm(dm)
      use amrex_distromap_module,     ONLY : amrex_distromap
      use gr_physicalMultifabs, ONLY : unk
      implicit none
      type(amrex_distromap ), intent(inout) :: dm(:)
      dm=unk%dm
      return
end subroutine gr_amrexGetDm

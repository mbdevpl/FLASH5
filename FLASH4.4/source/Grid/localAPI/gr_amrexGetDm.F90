!!****if* source/Grid/localAPI/gr_amrexGetDm
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
      implicit none
      type(amrex_distromap), intent(inout) :: dm(:)
      return
end subroutine gr_amrexGetDm

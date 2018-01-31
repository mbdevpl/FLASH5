!!****if* source/Grid/localAPI/gr_amrexLsSolvePoissonUnk
!!
!!  NAME 
!!
!!  gr_amrexLsSolvePoissonUnk
!!
!!  SYNOPSIS
!!
!!  call gr_amrexLsSolvePoissonUnk()
!!
!!  DESCRIPTION 
!!      This routine solves AX=B using Linear Solvers provided in AMReX on Unk multifab
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!  Stub implementation.
!!
!!***


subroutine gr_amrexLsSolvePoissonUnk (iSoln, iSrc, bcTypes, bcValues, poisfact)
  
  implicit none
    integer, intent(in)    :: iSoln, iSrc
    integer, intent(in)    :: bcTypes(6)
    real, intent(in)       :: bcValues(2,6)
    real, intent(inout)    :: poisfact
  return
  
end subroutine gr_amrexLsSolvePoissonUnk

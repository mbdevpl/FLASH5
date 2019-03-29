!!****f* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Burn_finalize
!!
!! NAME
!!  
!!  Burn_finalize
!!
!!
!! SYNOPSIS
!! 
!!  call Burn_finalize()
!!
!!  
!! DESCRIPTION
!!
!!  Finalizes the Burn module.
!!
!! NOTES
!!  
!!  There is no implementation that does anything.
!!
!!  The NSE arrays used with parametricBurn are deallocated by
!!  NSE_finalize, which should be called directly frm Driver_finalizeFlash.
!!
!!***


subroutine Burn_finalize()

  use bnNetwork_interface, ONLY : bn_finalizeNetwork

  implicit none

  call bn_finalizeNetwork

  return

end subroutine Burn_finalize

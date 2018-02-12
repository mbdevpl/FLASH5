!!****if* source/Grid/localAPI/gr_amrexLsInit
!!
!! NAME
!!
!!  gr_amrexLsInit
!!
!!
!! SYNOPSIS
!!
!!  call gr_amrexLsInit()
!!
!! Description
!!
!!  Initializes local data for amrex linear solvers.
!!  All the variables here are initialized by calling the
!!  RuntimeParameters_get subroutine. These data variables are for
!!  Unit Scope ->  Diffuse.
!!
!! ARGUMENTS
!!
!!  none  
!!
!! PARAMETERS
!!
!!***

subroutine gr_amrexLsInit() 
  
  implicit none  
  
  
  return
  
end subroutine gr_amrexLsInit

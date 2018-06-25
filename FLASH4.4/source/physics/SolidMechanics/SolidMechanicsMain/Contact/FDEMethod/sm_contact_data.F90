!!****if* source/physics/SolidMechanics/SolidMechanicsMain/Contact/FDEMethod/sm_contact_data
!!
!! NAME
!!  gr_sbDistributedForces
!!
!! SYNOPSIS
!!
!!  gr_sbDistributedForces()
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!***
module sm_contact_data

  implicit none

  type contact_structure
     
     integer, allocatable,dimension(:,:) :: intersectedNodes
     real, allocatable,dimension(:,:) :: intersectForce
     integer, allocatable,dimension(:) :: intersect_count
     real, allocatable,dimension(:) :: intersect_vol
     !integer :: intersect_count
  end type contact_structure

  type(contact_structure),save,dimension(:),pointer :: sm_ContactInfo 
  

end module sm_contact_data

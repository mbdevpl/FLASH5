!!****if* source/diagnostics/ProtonEmission/localAPI/pem_setupEmissionBoxes
!!
!! NAME
!!
!!  pem_setupEmissionBoxes
!!
!! SYNOPSIS
!!
!!  call pem_setupEmissionBoxes ()
!!
!! DESCRIPTION
!!
!!  Sets up the emission boxes. These are defined through their bounding box coordinates (lower
!!  and upper corners). A check is done on each emission box to see if it has a positive volume.
!!  In case the user has not defined any emission box, the routine does nothing and exits right
!!  away.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  The emission boxes need not be entirely inside the domain. They can also overlap in volume,
!!  i.e. an active emission region can be defined by more than one emission box.
!!
!!***

subroutine pem_setupEmissionBoxes ()

  implicit none

  return
end subroutine pem_setupEmissionBoxes

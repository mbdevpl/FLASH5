!!****if* source/diagnostics/ProtonEmission/localAPI/pem_openDetectorFiles
!!
!! NAME
!!
!!  pem_openDetectorFiles
!!
!! SYNOPSIS
!!
!!  call pem_openDetectorFiles ()
!!
!! DESCRIPTION
!!
!!  Opens the detector files for recording emission screen protons. This is done at
!!  initialization of a run, since the emission protons will be recorded on the same
!!  screen during the entire simulation. To each detector there corresponds one file,
!!  characterized by the detector number. Only the master processor opens the detector
!!  files. The names of the detector files is as follows:
!!
!!              <basenm> EmissionProtonDetector <detectorID>
!!
!!  where <basenm> is the simulation base name and <detectorID> contains the detector
!!  number.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!          
!!  The detector number range is currently limited to 01-99.
!!
!!***

subroutine pem_openDetectorFiles ()

  implicit none

  return
end subroutine pem_openDetectorFiles

!!****if* source/diagnostics/ProtonEmission/localAPI/pem_write2DetectorFile
!!
!! NAME
!!
!!  pem_write2DetectorFile
!!
!! SYNOPSIS
!!
!!  call pem_write2DetectorFile (integer, intent (in) :: detector,
!!                               integer, intent (in) :: recordCount,
!!                               integer, intent (in) :: bucketCount,
!!                               real,    intent (in) :: bucket (:,:))
!!
!! DESCRIPTION
!!
!!  Writes a bucket of emission screen protons to the appropriate detector file on disk.
!!  The file is already open and the screen protons data is appended. Only the master processor
!!  writes out the bucket.
!!
!! ARGUMENTS
!!
!!  detector      : The detector number
!!  recordCount   : Number of records in bucket (2 only: screen [x,y] pairs)
!!  bucketCount   : Number of elements in bucket
!!  bucket        : The bucket containing all the elements.
!!
!!***

subroutine pem_write2DetectorFile (detector, recordCount, bucketCount, bucket)

  implicit none

  integer, intent (in) :: detector
  integer, intent (in) :: recordCount
  integer, intent (in) :: bucketCount
  real,    intent (in) :: bucket (1:recordCount,1:bucketCount)

  return
end subroutine pem_write2DetectorFile

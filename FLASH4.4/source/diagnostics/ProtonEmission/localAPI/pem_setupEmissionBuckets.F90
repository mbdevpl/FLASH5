!!****if* source/diagnostics/ProtonEmission/localAPI/pem_setupEmissionBuckets
!!
!! NAME
!!
!!  pem_setupEmissionBuckets
!!
!! SYNOPSIS
!!
!!  call pem_setupEmissionBuckets ()
!!
!! DESCRIPTION
!!
!!  Sets up the proton emission buckets. A bucket is defined as a set of emission protons that
!!  fit safely into the given proton memory. Each bucket is given a consecutive integer number
!!  starting from 1. At any given time only one bucket of protons is being actively created and
!!  transported through the domain.
!!
!!  Procedure:
!!
!!  The predetermined count of emission protons on each processor is stored and broadcast to
!!  all processors. All processors set up the emission buckets. There are two situations that
!!  can happen:
!!
!!    1) One/several processors can be treated completely in one bucket.
!!
!!    2) The number of emitted protons on one processor is larger than
!!       the available memory -> the code handles partial bucket processing
!!       as well.
!!
!!  Only one bucket can be handled partially at a time, so in 1) only buckets that can be
!!  treated completely are included. As an example, suppose we have 8 processors with the
!!  following emission proton count:
!!
!!              proc | proton count
!!            ----------------------
!!               0   |   1,000
!!               1   |  20,000
!!               2   |     500
!!               3   |   7,000
!!               4   |      50
!!               5   | 140,000
!!               6   |  40,000
!!               7   |   6,000
!!
!!  and say we have 30,000 proton memory available. Then the first bucket consists of
!!  processors [0-4], the second bucket (partially handled) of processor [5], the third
!!  bucket of (partially handled) of processor [6] and the last bucket of (completely
!!  handled) processor [7]. For a memory of 200,000 we would just have two buckets, which
!!  are handled completely: first bucket [0-5] and second bucket [6,7].
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  none
!!
!!***

subroutine pem_setupEmissionBuckets ()

  implicit none

  return
end subroutine pem_setupEmissionBuckets

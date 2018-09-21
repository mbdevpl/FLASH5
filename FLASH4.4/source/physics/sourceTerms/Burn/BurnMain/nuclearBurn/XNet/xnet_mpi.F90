!***************************************************************************************************
! xnet_mpi.f90 06/27/18
! These wrappers are used so that a non-MPI version can coexist with the MPI version.
! Ported from https://github.com/AMReX-Codes/amrex/blob/master/Src/F_BaseLib/parallel.f90 
!***************************************************************************************************

module xnet_mpi
  implicit none
  include 'mpif.h'
  !use mpi
end module xnet_mpi

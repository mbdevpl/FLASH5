!     
! File:   sm_pk_data.F90
! Author: tim
!
!

#include "SolidMechanics.h"
 
module sm_pk_data

      implicit none

      type sm_pk_dataset
           integer :: flag      ! which kinematics is this
           integer :: NumParams ! Size of params
           real, allocatable, dimension(:) :: params
      end type sm_pk_dataset

      integer :: sm_pk_NumKinematics
      type(sm_pk_dataset), save, dimension(:), pointer :: sm_pk_info
      real :: sm_pk_timedelay
      
end module sm_pk_data


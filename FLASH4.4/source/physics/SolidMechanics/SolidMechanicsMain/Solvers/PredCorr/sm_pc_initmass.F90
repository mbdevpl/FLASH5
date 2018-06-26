#include "SolidMechanics.h"

subroutine sm_pc_initmass(ibd)
      use sm_assemble_interface, only: sm_assemble_mass
      use SolidMechanics_data, only :  sm_BodyInfo, sm_structure
      use sm_Misc_interface, only: sm_c2fortran_dgssv
      implicit none

      ! IO
      integer, intent(in) :: ibd ! body number

      ! Internal variables
      type(sm_structure), pointer :: body
      integer :: iopt, ldb, info, nrhs

      ! Get the body
      body => sm_BodyInfo(ibd)

      ! initalize the mass, body%M and body%Mqv
      call sm_assemble_mass(ibd)

      ! init dyn_rhs variable
      body%dyn_rhs = 0.

      ! Build LU = M, and store in LU_factors
      iopt = 1
      nrhs = 1
      ldb = body%neq
      write(*,*) 'init mass'
      call sm_c2fortran_dgssv( iopt, body%neq, body%qq_nnz, nrhs, body%M, &
                            body%qq_ia, body%qq_ja, body%dyn_rhs, ldb, &
                            body%lu_factors, info )


end subroutine sm_pc_initmass


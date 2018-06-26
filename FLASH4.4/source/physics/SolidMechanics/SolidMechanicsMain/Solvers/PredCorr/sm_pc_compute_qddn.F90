#include "SolidMechanics.h"

subroutine sm_pc_compute_qddn(ibd)
      use SolidMechanics_data, only :  sm_BodyInfo, sm_structure
      use sm_Misc_interface, only: sm_c2fortran_dgssv,  DGEMV_sparse, sm_ludcmp, sm_lubksb
      implicit none

      ! IO
      integer, intent(in) :: ibd ! body number
 
      ! Internal variables
      type(sm_structure), pointer :: body
      integer :: iopt, ldb, info, nrhs, neq,j, nstep
      real, allocatable, dimension(:) :: rhs
      integer, allocatable, dimension(:) :: indx
      real :: d

      ! Get the body
      body => sm_BodyInfo(ibd)

      neq = body%neq

      !
      ! Clear the dyn_rhs container
      !
      body%dyn_rhs(1:neq) = 0.


      if (body%bodyType .eq. BODYTYPE_RIGID ) then

         allocate(indx(neq))

         ! For now only works for all rigid body dofs free. 
         call sm_ludcmp(body%M_rigid(1:neq,1:neq),neq,neq,indx(1:neq),d)

         ! Add force vectors : Here Qs already have their rhs sign.
         ! Hack Set Hs to zero:         
         body%dyn_rhs(1:neq) = body%dyn_rhs(1:neq) + body%Qs(1:neq) + body%Hs(1:neq) 

         call sm_lubksb(body%M_rigid(1:neq,1:neq),neq,neq,indx(1:neq),body%dyn_rhs(1:neq))

         deallocate(indx)

      else


      ! if there are applied kinematics, then 
      ! $rhs = rhs - M_{qv}*\ddot{v}$
      if( allocated( body%restraints_surf ) ) then
         call DGEMV_sparse(body%ng,body%qv_nnz,body%Mqv,  &
                           body%qv_IA,body%qv_JA,         &
                           -1.0, body%qddn(neq+1), 1.0, &  !note that qddn must contain the vdd values for time n
                           body%dyn_rhs )

         if( body%damping_flag == SM_TRUE ) then
            call DGEMV_sparse(body%ng,body%qv_nnz,body%Dampqv,  &
                              body%qv_IA,body%qv_JA,            &
                              -1.0, body%qdn(neq+1), 1.0,     & 
                              body%dyn_rhs )
            
         end if

      end if

      !
      ! If there is Prop. Damping
      !
      if( body%damping_flag == SM_TRUE ) then
         call DGEMV_sparse(body%neq,body%qq_nnz,body%damp, &
              body%qq_IA, body%qq_JA,                      &
              -1.0, body%qdn, 1.0,                       & 
              body%dyn_rhs )    
      end if

      ! Reduce 
      ! rhs = rhs + Fext - Fint
      ! Fext: body%Hs: external forces
      !                body forces (like gravity), and surface forces from the fluid
      ! Fint: body%Qs: internal forces
      !                elastic forces, etc
      body%dyn_rhs(1:neq) = body%dyn_rhs(1:neq) + body%Hs(1:neq) - body%Qs(1:neq)

      ! Build LU = M, and store in LU_factors
      iopt = 2
      nrhs = 1
      ldb  = neq
      call sm_c2fortran_dgssv( iopt, body%neq, body%qq_nnz, nrhs, body%M, &
                               body%qq_ia, body%qq_ja, body%dyn_rhs, ldb, &
                               body%lu_factors, info )

      endif

      ! update qddn with dyn_rhs
      body%qddn(1:neq) = body%dyn_rhs(1:neq)

      return

end subroutine sm_pc_compute_qddn


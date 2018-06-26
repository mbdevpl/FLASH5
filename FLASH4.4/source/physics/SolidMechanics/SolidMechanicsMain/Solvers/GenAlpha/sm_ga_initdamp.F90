#include "SolidMechanics.h"

subroutine sm_ga_initdamp(ibd)
  use sm_assemble_interface, only: sm_assemble_stiff
  use SolidMechanics_data, only :  sm_BodyInfo, sm_structure
  implicit none

  ! IO
  integer, intent(in) :: ibd ! body number

  ! Internal variables
  type(sm_structure), pointer :: body
  integer :: nnz,ndofs, neq
  real    :: kappa_M, kappa_K0
  real, allocatable, dimension(:) :: qi, qn, Qs

  ! Get the body
  body => sm_BodyInfo(ibd)

  write(*,*) ' body%damping_flag = ', body%damping_flag
  if( body%damping_flag /= SM_TRUE ) then
     return
  end if

  write(*,*) 'init damping'

  kappa_M = Body%damping_kappa_M
  kappa_K0= Body%damping_kappa_K0

  ! Check to see if building K is needed
  if( abs(kappa_K0) > 1.e-14 ) then

     ! initalize the stiffness matrix in the reference configuration, body%K and body%Kqv
     ! *** sm_assemble_stiff uses qi, loaded IC's are in qn and qi at this point
     ndofs = body%ndofs
     neq   = body%neq
     allocate( qi(ndofs), qn(ndofs), Qs(neq) )
     ! save a copy of qi and qn from body
     qi(1:ndofs) = body%qi(1:ndofs)
     qn(1:ndofs) = body%qn(1:ndofs)
     Qs(1:neq) = body%Qs(1:neq)
     ! zero qi and qn out for reference configuration
     body%qi(1:ndofs) = 0.
     body%qn(1:ndofs) = 0.
     ! compute stiffness K at reference (q = 0)
     call sm_assemble_stiff(ibd)
     ! replace the values of qi and qn
     body%qi(1:ndofs) = qi(1:ndofs)
     body%qn(1:ndofs) = qn(1:ndofs)
     ! Get rid of qi and qn
     deallocate( qi, qn )

     ! Compute Damping matrix
     ! Damp = k_M*[M] + k_K0*[K]
     nnz = body%qq_nnz
     body%Damp(1:nnz) = kappa_M*body%M(1:nnz) + kappa_K0*body%K(1:nnz)

     nnz = body%qv_nnz
     body%Dampqv(1:nnz) = kappa_M*body%Mqv(1:nnz) + kappa_K0*body%Kqv(1:nnz)

     ! Rebuild the stiffness matrix with the real qn
     !call sm_assemble_stiff(ibd)
     body%Qs(1:neq) = Qs(1:neq)

  else
     ! K is not needed

     ! Compute Damping matrix
     ! Damp = k_M*[M]
     nnz = body%qq_nnz
     body%Damp(1:nnz) = kappa_M*body%M(1:nnz)

     nnz = body%qv_nnz
     body%Dampqv(1:nnz) = kappa_M*body%Mqv(1:nnz)
  end if

  return

end subroutine sm_ga_initdamp


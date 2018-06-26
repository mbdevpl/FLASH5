!     Stub Function
! File:  
! Author: tim

subroutine sm_assemble_mass_rigid(ibd, flag)

  use Driver_interface, only : Driver_abortFlash
  use SolidMechanics_Data, only : sm_nen, sm_BodyInfo, sm_structure
  use sm_element_interface, only : el15_mass

  implicit none
#include "Flash.h"    
#include "constants.h"
#include "SolidMechanics.h"

  !IO Variables
  integer, intent(in)    :: ibd, flag ! body number
  

  ! Local Variables:
  type(sm_structure), pointer :: body
  real :: mass
  real :: I_newton(NDIM,NDIM)  
  real, allocatable, dimension(:,:) :: M
  integer, allocatable, dimension(:) :: dof_glob,dof_loc,rest_glob,rest_loc
  integer :: neq,ng,nmaxdofs,ndof,nrest,i,imaster

  body => sm_BodyInfo(ibd)

  ! Max dofs in nodes of body
  nmaxdofs = body%max_dofs_per_node
  allocate(M(nmaxdofs,nmaxdofs))
  allocate(dof_glob(nmaxdofs),dof_loc(nmaxdofs),rest_glob(nmaxdofs),rest_loc(nmaxdofs))

  ! Unrestrained coordinates:
  neq = body%neq
  body%M_rigid = 0.
  ! Restrained Coordinates:
  ng  = body % ng  
  body%Mqv_rigid = 0.

  ! Body Mass:
  mass = body%mass

  ! Body Inertia Moments in Global axes
#if NDIM == 2
  I_newton(1,1) = body%I_newton(1,1)
#elif NDIM == MDIM
  I_newton(1:NDIM,1:NDIM) = body%I_newton(1:NDIM,1:NDIM)
#endif

  imaster = body%borigin_node

  select case(body%eltype(imaster))
  case (ONE_NODE_POINT)

     ! Build mass matrix for all gen coords:
     call el15_mass(nmaxdofs,ibd,M)

     ! Figure out which of the BODYFRAME_NODE gen coords are dofs
     ! and which are restraints
     ndof = 0; nrest = 0;
     do i=1,nmaxdofs
        if (body%ID(i,BODYFRAME_NODE) .le. neq) then ! Dof

          ndof = ndof + 1
          dof_glob(ndof)=body%ID(i,BODYFRAME_NODE)
          dof_loc(ndof) = i

        else

          nrest = nrest + 1
          rest_glob(nrest)=body%ID(i,BODYFRAME_NODE)
          rest_loc(nrest) = i

        endif
     enddo
    
     rest_glob = rest_glob - neq
 
     ! Now add to free-free and free-restrained parts of matrix.
     if (ndof .eq. neq) then
         body%M_rigid(dof_glob(1:ndof),dof_glob(1:ndof)) =    &
              M(dof_loc(1:ndof),dof_loc(1:ndof))
     else
        write(*,*) 'Body=',ibd,'ndof node 1=',ndof,', neq=',neq
        call Driver_abortFlash("sm_assemble_mass_rigid : ndof not equal to neq.")
     endif
     if ((nrest .gt. 0) .and. (nrest .le. ng)) then
         body%Mqv_rigid(dof_glob(1:ndof),rest_glob(1:nrest)) = &
              M(dof_loc(1:ndof),rest_loc(1:nrest)) 
     endif

     !call Driver_abortFlash('sm_assemble_mass_rigid: mass assembly not completed')
     !body % M_rigid(1:neq,1:neq) = M(1:neq,1:neq)     ! All coords unrestrained for now.
     !body % Mqv_rigid(1:ng,1:ng) = Mqv(1:ng,1:ng) 

  case default

     call Driver_abortFlash('sm_assemble_mass_rigid: mass assembly not defined')

  end select

  deallocate(M,dof_glob,dof_loc,rest_glob,rest_loc)
  
  return
    
end subroutine sm_assemble_mass_rigid


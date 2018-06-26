! For 3D Flexible Body
! Get the maximum eigenvalue from an element

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_3DFlexible_getElement_EvalMax(lambda_max, body, e)

  USE SolidMechanics_data, only: sm_nen, sm_structure
  USE sm_element_interface, only: el05_mass, el05_stiff_biot, el05_stiff_kirch, &
                                  el12_mass, el12_stiff_biot, el12_stiff_kirch
  USE Driver_interface, only: Driver_abortFlash

  implicit none
    
  ! IO Variables
  real,   intent(out) :: lambda_max
  type(sm_structure), intent(in)  :: body
  integer,            intent(in)  :: e

  ! Define internal variables
  integer :: nee, info, nen_e
  real, allocatable, dimension(:,:) :: Me, ke, XYZ_e, qn_e
  real, allocatable, dimension(:)   :: qes, lambda, work
  integer :: lwork
  integer, parameter :: liwork = 1
  integer, dimension(liwork) :: iwork
  character*1 :: UPLO = 'L'    

  nen_e = sm_nen( body%eltype(e) )
  nee = 3*nen_e

  allocate( me(nee, nee), ke(nee, nee), qes(nee), lambda(nee) )   

  ! Get XYZ for the element
  allocate(XYZ_e(nen_e,3), qn_e(nen_e,3))
  call get_Nodal_XYZ(e, nen_e, XYZ_e, body)
  call get_Nodal_UVW_qn(e, nen_e, qn_e, body) 

  lwork = 2*nee+1
  allocate( work(lwork) )
  work = 0.

  select case( body%eltype(e) )

  case( EIGHT_NODE_HEXAHEDRON )

     call el05_mass(Me, XYZ_e, body%MatDensity(e))

     select case( body%MatType(e) )

     case( MATERIAL_KIRCHHOFF )

        call el05_stiff_Kirch( ke, qes, XYZ_e, &
             body%YoungsModulus(e),            &
             body%PoissonsRatio(e), qn_e)

     case( MATERIAL_BIOT )

        call el05_stiff_Biot( ke, qes, XYZ_e, &
             body%YoungsModulus(e),           &
             body%PoissonsRatio(e), qn_e)

     end select

  case( TWSEVEN_NODE_HEXAHEDRON )

     call el12_mass(Me, XYZ_e, body%MatDensity(e))

     select case( body%MatType(e) )

     case( MATERIAL_KIRCHHOFF )

        call el12_stiff_Kirch( ke, qes, XYZ_e, &
             body%YoungsModulus(e),            &
             body%PoissonsRatio(e), qn_e)

     case( MATERIAL_BIOT )

        call el12_stiff_Biot( ke, qes, XYZ_e, &
             body%YoungsModulus(e),           &
             body%PoissonsRatio(e), qn_e)

     end select

  case default

     call Driver_abortFlash("unknown element type")

  end select

  ! Compute the eigenvalues of this [Ke]{q} = lambda [Me]
  call DSYGVD( 1, 'N', UPLO, nee, ke, nee, Me, nee, lambda, work, lwork, iwork, liwork, info )

  ! get max eigenvalue
  lambda_max = lambda(nee)

  ! Cleanup the workspace
  deallocate(me, ke, qes)
  deallocate(lambda, work)

  return

end subroutine sm_3DFlexible_getElement_EvalMax

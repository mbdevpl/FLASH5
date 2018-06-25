Module sm_Misc_interface

#include "Flash.h"
#include "SolidMechanics.h"

    interface
        subroutine get_Nodal_XYZ(e, nen_e, Xe, body)
            use SolidMechanics_data, only: sm_structure
            implicit none
            integer, intent(in) :: e, nen_e
            real, intent(out), dimension(nen_e,3) :: Xe
            type(sm_structure), intent(in) :: body
        end subroutine
    end interface

    interface
        subroutine get_Nodal_UVW_qi(e, nen_e, Qe, body)
            use SolidMechanics_data, only: sm_structure
            implicit none
            integer, intent(in) :: e, nen_e
            real, intent(out), dimension(nen_e,3) :: Qe
            type(sm_structure) :: body
        end subroutine
    end interface


    interface
        subroutine get_Nodal_UVW_qn(e, nen_e, Qe, body)
            use SolidMechanics_data, only: sm_structure
            implicit none
            integer, intent(in) :: e, nen_e
            real, intent(out), dimension(nen_e,3) :: Qe
            type(sm_structure) :: body
        end subroutine
    end interface

    interface
        subroutine get_Nodal_UVW_qdn(e, nen_e, qde, body)
            use SolidMechanics_data, only: sm_structure
            implicit none
            integer, intent(in) :: e, nen_e
            real, intent(out), dimension(nen_e,3) :: qde
            type(sm_structure) :: body
        end subroutine
    end interface

    ! get_Nodal_UVW_qdn.F90 get_Nodal_UVW_qi.F90  get_Nodal_qms.F90     get_Nodal_qn.F90

    interface
        subroutine get_Nodal_UVW_qdm(e, nen_e, qde, body)
            use SolidMechanics_data, only: sm_structure
            implicit none
            integer, intent(in) :: e, nen_e
            real, intent(out), dimension(nen_e,3) :: qde
            type(sm_structure) :: body
        end subroutine
    end interface

    interface
        subroutine get_Nodal_qms(e, nen_e, Xe, body)
            use SolidMechanics_data, only: sm_structure
            implicit none
            integer, intent(in) :: e, nen_e
            real, intent(out), dimension(nen_e,3) :: Xe
            type(sm_structure), intent(in) :: body
        end subroutine
    end interface

    interface
        subroutine get_Nodal_qn(e, nen_e, Xe, body)
            use SolidMechanics_data, only: sm_structure
            implicit none
            integer, intent(in) :: e, nen_e
            real, intent(out), dimension(nen_e,3) :: Xe
            type(sm_structure), intent(in) :: body
        end subroutine
    end interface

    interface
        subroutine get_Nodal_UVW_qdi(e, nen_e, qde, body)
            use SolidMechanics_data, only: sm_structure
            implicit none
            integer, intent(in) :: e, nen_e
            real, intent(out), dimension(nen_e,3) :: qde
            type(sm_structure) :: body
        end subroutine
    end interface
    interface
        subroutine compute_Jacobian3d(nen_e, Xe, NNxi, NNeta, NNzeta, Jinv, detJ)
            implicit none
            integer, intent(in)                  :: nen_e
            real, intent(in), dimension(nen_e)   :: NNxi, NNeta, NNzeta
            real, intent(in), dimension(nen_e,3) :: Xe
            real, intent(out), dimension(3,3)    :: Jinv
            real, intent(out)                    :: detJ
        end subroutine
    end interface

    interface
        subroutine DGEMV_sparse(n,nnz,A,IA,JA,alpha,x,beta,y)
            ! compute y = alpha*A*x + beta*y
            implicit none
            integer, intent(in) :: n, nnz
            real, intent(in)    :: A(nnz)
#if FEM_MATFORMAT == FEM_CSR
            integer, intent(in) :: IA(n+1), JA(nnz)
#elif FEM_MATFORMAT == FEM_CSC
            integer, intent(in) :: IA(nnz), JA(n+1)
#endif
            real, intent(in)    :: alpha
            real, intent(in)    :: x(*)
            real, intent(in)    :: beta
            real, intent(inout) :: y(*)
        end subroutine
    end interface

    interface
        subroutine sm_deallocateBody ( ibd )
            implicit none
            integer, intent(in) :: ibd
        end subroutine
    end interface
    
    interface
        subroutine sm_c2fortran_dgssv( iopt, n, nnz, nrhs, values, rowind, colptr, b, ldb, f_factors, info ) BIND(C,NAME='sm_c2fortran_dgssv_')
                    implicit none
                    integer, intent(in)    :: iopt, n, nnz, nrhs, rowind(*), colptr(*), ldb
                    integer, intent(inout) :: info
                    real, intent(inout) :: values(*), b(*)
                    integer*8 f_factors  ! /* a handle containing the address pointing to the factored matrices */
        end subroutine
    end interface

    interface
        subroutine get_Nodal_XYZ_ws(e, nen_e, Xe, body)
            use SolidMechanics_data, only: sm_structure
            implicit none
            integer, intent(in) :: e, nen_e
            real, intent(out), dimension(nen_e,NDIM) :: Xe
            type(sm_structure), intent(in) :: body
        end subroutine
    end interface

    interface
        subroutine get_Nodal_UVW_qi_ws(e, nen_e, Qe, body)
            use SolidMechanics_data, only: sm_structure
            implicit none
            integer, intent(in) :: e, nen_e
            real, intent(out), dimension(nen_e,NDIM) :: Qe
            type(sm_structure) :: body
        end subroutine
    end interface

     interface
        subroutine get_Nodal_UVW_qn_ws(e, nen_e, Qe, body)
            use SolidMechanics_data, only: sm_structure
            implicit none
            integer, intent(in) :: e, nen_e
            real, intent(out), dimension(nen_e,NDIM) :: Qe
            type(sm_structure) :: body
        end subroutine
    end interface

    interface
       function sm_crossProd(a,b)
         implicit none
         real, intent(in) :: a(NDIM),b(NDIM)
#if NDIM == 2
         real :: sm_crossProd
#elif NDIM == 3
         real :: sm_crossProd(NDIM)
#endif
       end function sm_crossProd
    end interface

    interface
       function sm_crossProd_mat(a)
         implicit none
         real, intent(in) :: a(NDIM)
         real :: sm_crossProd_mat(NDIM,NDIM)
       end function sm_crossProd_mat
    end interface

    interface
       subroutine sm_get_NumBodies(numBodies)
         implicit none
         integer, intent(out) :: numBodies
       end subroutine sm_get_NumBodies
    end interface
      interface
       subroutine  machineepsilon (eps)
         implicit none
         double precision,intent(OUT) :: eps
         double precision :: MACHEPS=1.D0 
       end subroutine machineepsilon
    end interface
    
    interface
       Subroutine int2char(i,strng)
         implicit none
         integer, intent(IN):: i
         character(6),intent(OUT):: strng
       end subroutine int2char
    end interface

    interface
       subroutine sm_ludcmp(a,n,np,indx,d)
         INTEGER n,np,indx(n)
         REAL d,a(np,np)
       end subroutine sm_ludcmp
    end interface

    interface
       subroutine sm_lubksb(a,n,np,indx,b)
         INTEGER n,np,indx(n)
         REAL a(np,np),b(n)
       end subroutine sm_lubksb
    end interface

end Module sm_Misc_interface

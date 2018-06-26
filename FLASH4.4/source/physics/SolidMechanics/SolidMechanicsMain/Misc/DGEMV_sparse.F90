#include "SolidMechanics.h" 

subroutine DGEMV_sparse(n,nnz,A,IA,JA,alpha,x,beta,y)
    implicit none
        
#if FEM_MATFORMAT == FEM_CSR
    ! for a CSR matrix A, dense x,y
    ! compute y = alpha*A*x + beta*y
    
    ! I/O:
    integer, intent(in)             :: n, nnz
    real, intent(in)    :: A(nnz)
    integer, intent(in)             :: IA(n+1), JA(nnz)
    real, intent(in)    :: alpha
    real, intent(in)    :: x(*)
    real, intent(in)    :: beta
    real, intent(inout) :: y(*)
    
    ! Internal variables
    integer :: i, k1, k2
    real, external :: DDOT
    
!$omp parallel do private(i,k1,k2)    
    do i = 1,n
        k1 = IA(i)
        k2 = IA(i+1)-1
        y(i) = alpha*ddot(n, A(k1:k2), 1, x(JA(k1:k2)), 1) + beta*y(i)
    enddo
    
#elif FEM_MATFORMAT == FEM_CSC
    ! for a CSC matrix A, dense x,y
    ! compute y = alpha*A*x + beta*y
    
    ! I/O:
    integer, intent(in)             :: n, nnz
    real, intent(in)    :: A(nnz)
    integer, intent(in)             :: IA(nnz), JA(n+1)
    real, intent(in)    :: alpha
    real, intent(in)    :: x(*)
    real, intent(in)    :: beta
    real, intent(inout) :: y(*)
    
    ! Internal variables
    integer :: i,j, k1, k2
    real, external :: DDOT
    
    ! Internal variables
    call dscal(n,beta,y,1)    
    do j = 1,n
        k1 = JA(j)
        k2 = JA(j+1)-1
        do i = k1,k2
            y(IA(i)) = y(IA(i)) + x(j)*A(i)*alpha
        enddo
    enddo
    
#endif
    
    
end
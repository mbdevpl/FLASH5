subroutine det33( A, det )
    implicit none
    real, intent(in), dimension(3,3) :: A
    real, intent(out) :: det
        
    det = A(1,3)*(-(A(2,2)*A(3,1)) + A(2,1)*A(3,2)) &
        + A(1,2)*(A(2,3)*A(3,1) - A(2,1)*A(3,3))    &
        + A(1,1)*(-(A(2,3)*A(3,2)) + A(2,2)*A(3,3))
    
end
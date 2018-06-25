subroutine compute_Jacobian3d(nen_e, Xe, NNxi, NNeta, NNzeta, Jinv, detJ)
        
    implicit none
    
    ! input:
    integer, intent(in)                              :: nen_e
    real, intent(in), dimension(nen_e)   :: NNxi, NNeta, NNzeta
    real, intent(in), dimension(nen_e,3) :: Xe
    
    ! output:
    real, intent(out), dimension(3,3) :: Jinv
    real, intent(out)                 :: detJ
    
    ! Internal:
    real, dimension(3)   :: J1, J2, J3
    
    
    ! Compute Jacobian
    !J = [Nxi*Xe, Neta*Xe, Nzeta*Xe]
    !detJ = det(J)
    call DGEMV('T', nen_e, 3, 1.d0, Xe, nen_e, NNxi,   1, 0.d0, J1, 1)
    call DGEMV('T', nen_e, 3, 1.d0, Xe, nen_e, NNeta,  1, 0.d0, J2, 1)
    call DGEMV('T', nen_e, 3, 1.d0, Xe, nen_e, NNzeta, 1, 0.d0, J3, 1)

    detJ = +J1(1)*(J2(2)*J3(3) - J3(2)*J2(3)) &
            -J2(1)*(J1(2)*J3(3) - J3(2)*J1(3)) &
            +J3(1)*(J1(2)*J2(3) - J2(2)*J1(3))
    Jinv(1:3,1) = (/ -J2(3)*J3(2) + J2(2)*J3(3),  J1(3)*J3(2) - J1(2)*J3(3), -J1(3)*J2(2) + J1(2)*J2(3) /)
    Jinv(1:3,2) = (/  J2(3)*J3(1) - J2(1)*J3(3), -J1(3)*J3(1) + J1(1)*J3(3),  J1(3)*J2(1) - J1(1)*J2(3) /)
    Jinv(1:3,3) = (/ -J2(2)*J3(1) + J2(1)*J3(2),  J1(2)*J3(1) - J1(1)*J3(2), -J1(2)*J2(1) + J1(1)*J2(2) /)
end
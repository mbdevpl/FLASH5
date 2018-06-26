subroutine get_Nodal_UVW_qdn(e, nen_e, qde, body)
    use SolidMechanics_data, only: sm_structure
    implicit none
    
    integer, intent(in) :: e, nen_e
    real, intent(out), dimension(nen_e,3) :: qde
    type(sm_structure) :: body
    
    integer :: a, idx1, idx2
    
    ! get qde = [ue, ve, we]
    do a = 1,nen_e
        
        idx1 = body%IEN(a,e)
        
        idx2 = body%ID(1,idx1)
        qde(a,1) = body%qdn(idx2)
        
        idx2 = body%ID(2,idx1)
        qde(a,2) = body%qdn(idx2)
        
        idx2 = body%ID(3,idx1)
        qde(a,3) = body%qdn(idx2)
        
    enddo
    
end
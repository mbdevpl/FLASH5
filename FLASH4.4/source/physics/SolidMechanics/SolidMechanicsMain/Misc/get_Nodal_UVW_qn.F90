subroutine get_Nodal_UVW_qn(e, nen_e, Qe, body)
    use SolidMechanics_data, only: sm_structure
    implicit none
    
    integer, intent(in) :: e, nen_e
    real, intent(out), dimension(nen_e,3) :: Qe
    type(sm_structure) :: body
    
    integer :: a, idx1, idx2
    
    ! get qe = [ue, ve, we]
    do a = 1,nen_e
        
        idx1 = body%IEN(a,e)
        
        idx2 = body%ID(1,idx1)
        Qe(a,1) = body%qn(idx2)
        
        idx2 = body%ID(2,idx1)
        Qe(a,2) = body%qn(idx2)
        
        idx2 = body%ID(3,idx1)
        Qe(a,3) = body%qn(idx2)
        
    enddo

end subroutine get_Nodal_UVW_qn

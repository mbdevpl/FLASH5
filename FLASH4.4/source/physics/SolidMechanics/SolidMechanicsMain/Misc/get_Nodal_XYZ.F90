subroutine get_Nodal_XYZ(e, nen_e, Xe, body)
    use SolidMechanics_data, only: sm_structure
    implicit none
    
    integer, intent(in) :: e, nen_e
    real, intent(out), dimension(nen_e,3) :: Xe
    type(sm_structure), intent(in) :: body
   
    integer :: a, idx
    
    ! get Xe = [xe, ye, ze]
    do a = 1,nen_e
        idx = body%IEN(a,e)
        Xe(a,1) = body%x(idx)
        Xe(a,2) = body%y(idx)
        Xe(a,3) = body%z(idx)
    enddo
    
end

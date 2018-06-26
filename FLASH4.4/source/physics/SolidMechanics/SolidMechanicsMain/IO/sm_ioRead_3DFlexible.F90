!     
! File:   stub function
! Author: tim
! Edit by Hussein Ezzeldin hmezz@gwu.edu Jan 2013 

subroutine sm_ioRead_3DFlexible(ibd,file)

  use Driver_interface, only: Driver_abortFlash

  USE HDF5
    
  implicit none
    
  integer, intent(IN) :: ibd
  INTEGER(HID_T), intent(in) :: file

  call Driver_abortFlash('3DFlexible IO not configured')
  
end subroutine sm_ioRead_3DFlexible




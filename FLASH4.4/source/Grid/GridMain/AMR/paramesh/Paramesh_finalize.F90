!!****if* source/Grid/GridMain/paramesh/Paramesh_finalize
!!
!! NAME
!!  Paramesh_finalize
!! 
!! SYNOPSIS
!!  call Paramesh_finalize()
!!
!! DESCRIPTION
!!  Simple interface to finalize the mesh.
!!  Each mesh package will have a different copy of this function.
!!
!! USED BY
!!  source/driver/end_flash
!!***
subroutine Paramesh_finalize()
  use paramesh_interfaces, ONLY : amr_close
  implicit none
  call amr_close()
end subroutine Paramesh_finalize


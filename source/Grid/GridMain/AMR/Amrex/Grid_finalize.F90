!!****if* source/Grid/GridMain/Chombo/Grid_finalize
!!
!! NAME
!!
!!  Grid_finalize
!!
!!
!! SYNOPSIS
!!
!!  Grid_finalize()
!!
!!
!! DESCRIPTION
!!
!!  Deallocates memory that has been allocated in the Grid Unit
!!
!!***


subroutine Grid_finalize()
    use gr_bcInterface,    ONLY : gr_bcFinalize
    use gr_ptInterface,    ONLY : gr_ptFinalize
    use gr_amrexInterface, ONLY : gr_amrexFinalize

    implicit none

    call gr_solversFinalize()
    call gr_ptFinalize()
    call gr_bcFinalize()
    call gr_amrexFinalize()
end subroutine Grid_finalize

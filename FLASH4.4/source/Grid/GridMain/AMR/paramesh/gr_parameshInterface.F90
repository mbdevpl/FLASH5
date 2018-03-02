!!****ih* source/Grid/localAPI/gr_parameshInterface
!!
!! NAME
!!
!!  gr_parameshInterface
!!
!! SYNOPSIS
!!
!!  use gr_parameshInterface
!!
!! DESCRIPTION
!!
!!  Interfaces for some subprograms private to the paramesh subunit
!!  implementation.
!!
!!***

module gr_parameshInterface
  implicit none

  interface
     function gr_blockMatch(blkID, ntype, refinementLevel) result(match)
       implicit none
       integer, intent(IN)           :: blkID
       integer, intent(IN)           :: ntype
       integer, intent(IN), OPTIONAL :: refinementLevel
       logical                       :: match
     end function gr_blockMatch
  end interface

end module gr_parameshInterface


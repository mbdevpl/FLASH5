!!****ih* source/Grid/localAPI/gr_amrextInterface
!!
!! NAME
!!
!!  gr_amrextInterface
!!
!! SYNOPSIS
!!
!!  use gr_amrextInterface
!!
!! DESCRIPTION
!!
!!  Interfaces for some subprograms related to AMReX transition and private to the GridMain subunit.
!!
!!
!!***


module gr_amrextInterface
#include "constants.h"
#include "Flash.h"
  implicit none

#ifdef FLASH_GRID_ANYAMREX
  interface
     subroutine gr_amrextBuildMultiFabsFromF4Grid(phi_mf, maxLev, nodetype)
       use amrex_multifab_module, ONLY : amrex_multifab
       implicit none
       type(amrex_multifab),intent(OUT) :: phi_mf(:)
       integer,intent(IN) :: maxLev
       integer,intent(IN),OPTIONAL :: nodetype
     end subroutine gr_amrextBuildMultiFabsFromF4Grid
  end interface
#endif

  interface
     subroutine gr_fillMetaData(blockID, blockDesc)
       use block_metadata, ONLY : block_metadata_t
       implicit none
       integer,intent(IN) :: blockID
       type(block_metadata_t),intent(OUT) :: blockDesc
     end subroutine gr_fillMetaData
  end interface

end module gr_amrextInterface

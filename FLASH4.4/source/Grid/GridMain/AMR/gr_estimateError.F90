!!****if* source/Grid/GridMain/AMR/gr_estimateError
!!
!! NAME
!!  gr_estimateError
!!
!! SYNOPSIS
!!
!!  gr_estimateError(real(INOUT) :: error(MAXBLOCKS),
!!                   integer(IN) :: iref,
!!                   real(IN)    :: refine_filter)
!!  
!!  DESCRIPTION
!!  
!!  For each block, estimate the error associated with the given variable to
!!  help determine if the block needs refinement or derefinement.  Update the
!!  corresponding value in the error array to be the maximum of the incoming
!!  value and the value calculated here.
!!
!!  ARGUMENTS 
!!
!!    error - indexed by block IDs.
!!
!!    iref - index of the refinement variable in data structure "unk"
!!
!!    refine_filter - makes sure that error calculations to determine refinement
!!                    don't diverge numerically 
!! 
!!  NOTES
!!  
!!    See Grid_markRefineDerefine
!!
!!  SEE ALSO
!!  
!!    Grid_markRefineDerefine
!!
!!***

subroutine gr_estimateError(error, iref, refine_filter)

  use Grid_data, ONLY: gr_geometry, &
       gr_meshComm, gr_meshMe,gr_delta, gr_domainBC
  use Grid_interface, ONLY : Grid_getBlkBC, &
                             Grid_getBlkPtr, Grid_releaseBlkPtr
  use gr_interface,   ONLY : gr_estimateBlkError
  use gr_specificData, ONLY : gr_oneBlock
  use block_iterator, ONLY : block_iterator_t, destroy_iterator
  use block_metadata, ONLY : block_metadata_t

  implicit none

#include "Flash_mpi.h"
#include "Flash.h"
#include "constants.h"  
#ifdef INDEXREORDER
  integer,parameter::IX=1,IY=2,IZ=3
#else
  integer,parameter::IX=2,IY=3,IZ=4
#endif  
  integer, intent(IN) :: iref
  real, intent(IN) ::  refine_filter
  real,intent(INOUT) :: error(MAXBLOCKS)
  integer, parameter :: SQNDIM = NDIM*NDIM
  
  real,dimension(MDIM) ::  del, del_f, psize
  integer,dimension(MDIM) :: ncell
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC,face,bdry
  real,allocatable,dimension(:,:,:,:)::delu,delua

  real delu2(SQNDIM), delu3(SQNDIM), delu4(SQNDIM)

  real num,denom

  integer i,j,k
  integer ierr,grd
  integer,dimension(MDIM)::bstart,bend 
  integer nsend,nrecv

  integer :: kk

  integer :: idest, iopt, nlayers, icoord
  logical :: lcc, lfc, lec, lnc, l_srl_only, ldiag
  type(block_iterator_t) :: itor
  type(block_metadata_t) :: blockDesc
  integer :: blkLevel, blkID

!==============================================================================

! A non-directional guardcell fill for CENTER (and also EOS calls for
! all block cells, including guardcells, if any refinement variables
! refine_var_# require this to be current) must have been performed
! when this routine is invoked. Moreover, there must not be any
! intervening calls that would modify the solution data in unk (at
! least, for the variables to be used for refinement criteria).
! Finally, this should be true for both LEAF and PARENT blocks
! (node types 1 and 2).
! Normally the caller of this routine, Grid_markRefineDerefine, takes care
! of all that.
!
! If this routine must be used in a situation where the conditions above
! are not true, the simplest (but probably inefficient) way of adapting
! this code to that situation would be uncommenting the following line:
!!$  call Grid_fillGuardCells(CENTER_FACES,ALLDIR)


! We are using more cell layers, including guardcells, from unk.

     
  !==============================================================================


  itor = block_iterator_t(ACTIVE_BLKS)
  do while(itor%is_valid())
     call itor%blkMetaData(blockDesc)

     blkID       = blockDesc%id
     call gr_estimateBlkError(error(blkID), blockDesc, iref, refine_filter)

     call itor%next()
  end do
#if defined(__GFORTRAN__) && (__GNUC__ <= 4)
  call destroy_iterator(itor)
#endif
end subroutine gr_estimateError


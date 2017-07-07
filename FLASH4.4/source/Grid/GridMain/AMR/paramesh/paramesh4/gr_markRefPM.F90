!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_markRefineDerefine
!!
!! NAME
!!  gr_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  gr_markRefineDerefine(integer(IN) :: iref,
!!                        real(IN) :: refine_cutoff,
!!                        real(IN) :: derefine_cutoff,
!!                        real(IN) :: refine_filter)
!!  
!!  DESCRIPTION
!!  
!!    Blocks are marked for refining or derefining.
!!    This version uses the second derivative calculations on the specified variable to 
!!    determine if the block needs more resoultion (refine) or less resolution (derefine).
!!    The arguments de/refine_cutoff are the thresholds for triggering the corresponding action.
!!
!!    After blocks have been marked, control is meant to be passed to Paramesh for actually
!!    updating the refinement of the grid.
!!
!!  ARGUMENTS 
!!
!!
!!    iref - index of the refinement variable in data structure "unk"
!!
!!    refine_cutoff - the threshold value for triggering refinement 
!!
!!    derefine_cutoff - the threshold for triggereing derefinement
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

!!REORDER(5): unk, unk1


subroutine gr_markRefPM(error, refine_cutoff,derefine_cutoff)

  use tree

  use Grid_data, ONLY: gr_meshComm, gr_meshMe, gr_maxRefine

  implicit none


#include "Flash_mpi.h"
#include "Flash.h"
#include "constants.h"  
  real, intent(IN) :: refine_cutoff, derefine_cutoff
  real,intent(IN) :: error(MAXBLOCKS)
  
  integer ierr
  integer nsend,nrecv
  integer reqr(MAXBLOCKS),reqs(MAXBLOCKS*nchild)
!
  integer :: lb,i,j,k
  integer :: statr(MPI_STATUS_SIZE,MAXBLOCKS)
  integer :: stats(MPI_STATUS_SIZE,MAXBLOCKS*nchild)
  real :: error_par(MAXBLOCKS)
  logical :: gcell_on_cc_backup(NUNK_VARS)


  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.


  
! MARK FOR REFINEMENT OR DEREFINEMENT

! first communicate error of parent to its children
! Children collect messages from parents.

  error_par(1:lnblocks) = 0.
  nrecv = 0
  do lb = 1,lnblocks
     if(parent(1,lb).gt.-1) then
        if (parent(2,lb).ne.gr_meshMe) then
           nrecv = nrecv + 1
           call MPI_IRecv(error_par(lb),1, &
                MPI_DOUBLE_PRECISION, &
                parent(2,lb), &
                lb, &
                gr_meshComm, &
                reqr(nrecv), &
                ierr)
        else
           error_par(lb) = error(parent(1,lb))
        end if
     end if
  end do
 
  ! parents send error to children

  nsend = 0
  do lb = 1,lnblocks
     do j = 1,nchild
        if(child(1,j,lb).gt.-1) then
           if (child(2,j,lb).ne.gr_meshMe) then
              nsend = nsend + 1
              call MPI_ISend(error(lb), &
                   1, &
                   MPI_DOUBLE_PRECISION, &
                   child(2,j,lb), &  ! PE TO SEND TO
                   child(1,j,lb), &  ! THIS IS THE TAG
                   gr_meshComm, &
                   reqs(nsend), &
                   ierr)
           end if
        end if
     end do
  end do

  if (nsend.gt.0) then
     call MPI_Waitall (nsend, reqs, stats, ierr)
  end if
  if (nrecv.gt.0) then
     call MPI_Waitall (nrecv, reqr, statr, ierr)
  end if

  do lb = 1,lnblocks

     if (nodetype(lb).eq.1 .or. nodetype(lb).eq.2) then
        
        ! test for derefinement
        
        if (nodetype(lb).eq.1 .and. .not.refine(lb).and..not.stay(lb) &
             &          .and.error(lb).le.derefine_cutoff &
             &          .and.error_par(lb).le.derefine_cutoff) then
           derefine(lb) = .TRUE.
        else
           derefine(lb) = .FALSE.
        end if
        
        ! test for refinement
        if (error(lb).gt.refine_cutoff) then
           derefine(lb) = .FALSE.
           refine(lb) = .TRUE.
        end if

        if (nodetype(lb).eq.1 .and. &
            error(lb).gt.derefine_cutoff.or.error_par(lb).gt.derefine_cutoff)  &
             &           stay(lb) = .TRUE.

        if (lrefine(lb).ge.gr_maxRefine)  &
             &           refine(lb) = .FALSE.


     end if
     
  end do

  !restore to the state when we came in
  ! When the flag arrays are passed to Paramesh for processing, only leaf
  ! blocks should be marked. - KW
  where (nodetype(:) .NE. LEAF)
     refine(:)   = .false.
     derefine(:) = .false.
  end where
  
  !=========================================================================
  return
end subroutine gr_markRefPM















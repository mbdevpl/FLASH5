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


subroutine gr_markRefineDerefine(&
                              iref,refine_cutoff,derefine_cutoff,refine_filter)

  use physicaldata, ONLY : gcell_on_cc, unk, unk1, no_permanent_guardcells
  use tree
  use Grid_data, ONLY: gr_geometry,  gr_maxRefine, &
       gr_meshComm, gr_meshMe,gr_delta
  use Grid_interface, ONLY : Grid_getblkIndexLimits
  use gr_specificData, ONLY : gr_oneBlock

  implicit none

#include "Flash_mpi.h"
#include "Flash.h"
#include "constants.h"  

  integer, intent(IN) :: iref
  real, intent(IN) :: refine_cutoff, derefine_cutoff, refine_filter
  integer, parameter :: SQNDIM = NDIM*NDIM
  
  real,dimension(MDIM) ::  del, del_f, psize
  integer,dimension(MDIM) :: ncell
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real,allocatable,dimension(:,:,:,:)::delu,delua

  real delu2(SQNDIM), delu3(SQNDIM), delu4(SQNDIM)

  real num,denom,error(MAXBLOCKS),error_par(MAXBLOCKS)

  integer lb,i,j,k
  integer ierr
  integer kstart,kend,jstart,jend,istart,iend
  integer nsend,nrecv
  integer reqr(MAXBLOCKS),reqs(MAXBLOCKS*nchild)
!
  integer :: kk
  integer :: statr(MPI_STATUS_SIZE,MAXBLOCKS)
  integer :: stats(MPI_STATUS_SIZE,MAXBLOCKS*nchild)

  real, pointer :: solnData(:,:,:)
  logical :: gcell_on_cc_backup(NUNK_VARS)
  integer :: idest, iopt, nlayers, icoord
  logical :: lcc, lfc, lec, lnc, l_srl_only, ldiag

!==============================================================================

  if (no_permanent_guardcells) gcell_on_cc_backup = gcell_on_cc

  if (.not. no_permanent_guardcells) then
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

     
!     work(:,:,:,:,1)=unk(iref,:,:,:,:)

  end if
  !==============================================================================



#define XCOORD(I) (gr_oneBlock(lb)%firstAxisCoords(CENTER,I))
#define YCOORD(I) (gr_oneBlock(lb)%secondAxisCoords(CENTER,I))
  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.


  do lb = 1,lnblocks


     error(lb) = 0.

     if (nodetype(lb).eq.1.or.nodetype(lb).eq.2) then

        solnData => unk(iref,:,:,:,lb)
        call Grid_getBlkIndexLimits(lb,blkLimits,blkLimitsGC)

        del=0.0
        ncell(:)=blkLimits(HIGH,:)-blkLimits(LOW,:)+1
        psize(:)=ncell(:)*gr_delta(:,lrefine(lb))
        del(IAXIS:NDIM) = 0.5e0*float(ncell(IAXIS:NDIM))/psize(IAXIS:NDIM)
        del_f(JAXIS:NDIM) = del(JAXIS:NDIM)
        allocate(delu(MDIM,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
 
        allocate(delua(MDIM,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
        ! Compute first derivatives
        
        do k = blkLimitsGC(LOW,KAXIS)+K3D*1,blkLimitsGC(HIGH,KAXIS)-K3D*1
           do j = blkLimitsGC(LOW,JAXIS)+K2D*1,blkLimitsGC(HIGH,JAXIS)-K2D*1
              do i = blkLimitsGC(LOW,IAXIS)+1,blkLimitsGC(HIGH,IAXIS)-1
                 
                 if (gr_geometry == SPHERICAL) &
                      del(IAXIS) = 1.0/(XCOORD(i+1) - XCOORD(i-1))
                 
                 ! d/dx
                 delu(1,i,j,k) = solnData(i+1,j,k) - solnData(i-1,j,k)
                 delu(1,i,j,k) = delu(1,i,j,k)*del(IAXIS)
                 
                 delua(1,i,j,k) = abs(solnData(i+1,j,k)) + &
                      abs(solnData(i-1,j,k))
                 delua(1,i,j,k) = delua(1,i,j,k)*del(IAXIS)
                 
#if N_DIM >= 2
                 ! d/dy
                 if ((gr_geometry == SPHERICAL) .or. (gr_geometry == POLAR)) then
                    del_f(JAXIS) = del(JAXIS)/XCOORD(i)
                 end if
                 
                 delu(2,i,j,k) = solnData(i,j+1,k) - solnData(i,j-1,k)
                 delu(2,i,j,k) = delu(2,i,j,k)*del_f(JAXIS)
                 
                 delua(2,i,j,k) = abs(solnData(i,j+1,k)) + &
                      abs(solnData(i,j-1,k))
                 delua(2,i,j,k) = delua(2,i,j,k)*del_f(JAXIS)
#endif
                 
#if N_DIM == 3
                 ! d/dz
                 if (gr_geometry == SPHERICAL) then
                    del_f(KAXIS) = del(KAXIS)/(  XCOORD(i) &
                              &    * sin(YCOORD(j))  )
                 else if (gr_geometry == CYLINDRICAL) then
                    del_f(KAXIS) = del(KAXIS)/XCOORD(i)
                 end if
                 delu(3,i,j,k) = solnData(i,j,k+1) -  solnData(i,j,k-1)
                 delu(3,i,j,k) = delu(3,i,j,k)*del_f(KAXIS)
                 
                 delua(3,i,j,k) = abs(solnData(i,j,k+1)) + &
                      abs(solnData(i,j,k-1))
                 delua(3,i,j,k) = delua(3,i,j,k)*del_f(KAXIS)
#endif
                 
              end do
           end do
        end do
        
        ! Compute second derivatives
        
        ! Test if at a block boundary
        
        ! Two guardcells
        kstart = 2*K3D+1
        kend   = ncell(KAXIS)+(K3D*((2*NGUARD)-2))
        jstart = 2*K2D+1
        jend   = ncell(JAXIS)+(K2D*((2*NGUARD)-2))
        istart = 3
        iend   = ncell(IAXIS)+(2*NGUARD)-2
        ! One guardcell
        !            kstart = 2*K3D+1+K3D
        !            kend   = ncell(KAXIS)+(K3D*((2*NGUARD)-2))-K3D
        !            jstart = 2*K2D+1+K2D
        !            jend   = ncell(JAXIS)+(K2D*((2*NGUARD)-2))-K2D
        !            istart = NGUARD
        !            iend   = ncell(IAXIS)+(2*NGUARD)-3
        ! No guardcells
        !            kstart = K3D*NGUARD+1
        !            kend   = ncell(KAXIS)+K3D*NGUARD
        !            jstart = K2D*NGUARD+1
        !            jend   = ncell(JAXIS)+K2D*NGUARD
        !            istart = NGUARD+1
        !            iend   = ncell(IAXIS)+NGUARD
        
        if (neigh(1,1,lb).le.-20) istart = NGUARD+1
        if (neigh(1,2,lb).le.-20) iend   = NGUARD+ncell(IAXIS)
        
#if N_DIM >= 2
        if (neigh(1,3,lb).le.-20) jstart = NGUARD*K2D+1
        if (neigh(1,4,lb).le.-20) jend   = NGUARD*K2D+ncell(JAXIS)
#endif
#if N_DIM == 3
        if (neigh(1,5,lb).le.-20) kstart = NGUARD*K3D+1
        if (neigh(1,6,lb).le.-20) kend   = NGUARD*K3D+ncell(KAXIS)
#endif
        
        do k = kstart,kend
           do j = jstart,jend
              do i = istart,iend
                 
                 if (gr_geometry == SPHERICAL) &
                      del(IAXIS) = 1.0/(XCOORD(i+1) - XCOORD(i-1))

                 ! d/dxdx
                 delu2(1) = delu(1,i+1,j,k) - delu(1,i-1,j,k)
                 delu2(1) = delu2(1)*del(IAXIS)
                 
                 delu3(1) = abs(delu(1,i+1,j,k)) + abs(delu(1,i-1,j,k))
                 delu3(1) = delu3(1)*del(IAXIS)
                 
                 delu4(1) = delua(1,i+1,j,k) + delua(1,i-1,j,k)
                 delu4(1) = delu4(1)*del(IAXIS)
                 
#if N_DIM >= 2
                 if ((gr_geometry == SPHERICAL) .or. (gr_geometry == POLAR)) then
                    del_f(JAXIS) = del(JAXIS)/XCOORD(i)
                 end if

                 ! d/dydx
                 delu2(2) = delu(1,i,j+1,k) - delu(1,i,j-1,k)
                 delu2(2) = delu2(2)*del_f(JAXIS)
                 
                 delu3(2) = abs(delu(1,i,j+1,k)) + abs(delu(1,i,j-1,k))
                 delu3(2) = delu3(2)*del_f(JAXIS)
                 
                 delu4(2) = delua(1,i,j+1,k) + delua(1,i,j-1,k)
                 delu4(2) = delu4(2)*del_f(JAXIS)
                 
                 ! d/dxdy
                 delu2(3) = delu(2,i+1,j,k) - delu(2,i-1,j,k)
                 delu2(3) = delu2(3)*del(IAXIS)
                 
                 delu3(3) = abs(delu(2,i+1,j,k)) + abs(delu(2,i-1,j,k))
                 delu3(3) = delu3(3)*del(IAXIS)
                 
                 delu4(3) = delua(2,i+1,j,k) + delua(2,i-1,j,k)
                 delu4(3) = delu4(3)*del(IAXIS)
                 
                 ! d/dydy
                 delu2(4) = delu(2,i,j+1,k) - delu(2,i,j-1,k)
                 delu2(4) = delu2(4)*del_f(JAXIS)
                 
                 delu3(4) = abs(delu(2,i,j+1,k)) +  &
                      &                          abs(delu(2,i,j-1,k))
                 delu3(4) = delu3(4)*del_f(JAXIS)
                 
                 delu4(4) = delua(2,i,j+1,k) + delua(2,i,j-1,k)
                 delu4(4) = delu4(4)*del_f(JAXIS)
#endif
                 
#if N_DIM == 3
                 if (gr_geometry == SPHERICAL) then
                    del_f(KAXIS) = del(KAXIS)/(  XCOORD(i) &
                              &    * sin(YCOORD(j))  )
                 else if (gr_geometry == CYLINDRICAL) then
                    del_f(KAXIS) = del(KAXIS)/XCOORD(i)
                 end if

                 ! d/dzdx
                 delu2(5) = delu(1,i,j,k+1) - delu(1,i,j,k-1)
                 delu2(5) = delu2(5)*del_f(KAXIS)
                 
                 delu3(5) = abs(delu(1,i,j,k+1)) + abs(delu(1,i,j,k-1))
                 delu3(5) = delu3(5)*del_f(KAXIS)
                 
                 delu4(5) = delua(1,i,j,k+1) + delua(1,i,j,k-1)
                 delu4(5) = delu4(5)*del_f(KAXIS)
                 
                 ! d/dzdy
                 delu2(6) = delu(2,i,j,k+1) - delu(2,i,j,k-1)
                 delu2(6) = delu2(6)*del_f(KAXIS)
                 
                 delu3(6) = abs(delu(2,i,j,k+1)) + abs(delu(2,i,j,k-1))
                 delu3(6) = delu3(6)*del_f(KAXIS)
                 
                 delu4(6) = delua(2,i,j,k+1) + delua(2,i,j,k-1)
                 delu4(6) = delu4(6)*del_f(KAXIS)
                 
                 ! d/dxdz
                 delu2(7) = delu(3,i+1,j,k) - delu(3,i-1,j,k)
                 delu2(7) = delu2(7)*del(IAXIS)
                 
                 delu3(7) = abs(delu(3,i+1,j,k)) + abs(delu(3,i-1,j,k))
                 delu3(7) = delu3(7)*del(IAXIS)
                 
                 delu4(7) = delua(3,i+1,j,k) + delua(3,i-1,j,k)
                 delu4(7) = delu4(7)*del(IAXIS)
                 
                 ! d/dydz
                 delu2(8) = delu(3,i,j+1,k) - delu(3,i,j-1,k)
                 delu2(8) = delu2(8)*del_f(JAXIS)
                 
                 delu3(8) = abs(delu(3,i,j+1,k)) + abs(delu(3,i,j-1,k))
                 delu3(8) = delu3(8)*del_f(JAXIS)
                 
                 delu4(8) = delua(3,i,j+1,k) + delua(3,i,j-1,k)
                 delu4(8) = delu4(8)*del_f(JAXIS)
                 
                 ! d/dzdz
                 delu2(9) = delu(3,i,j,k+1) - delu(3,i,j,k-1)
                 delu2(9) = delu2(9)*del_f(KAXIS)
                 
                 delu3(9) = abs(delu(3,i,j,k+1)) + abs(delu(3,i,j,k-1))
                 delu3(9) = delu3(9)*del_f(KAXIS)
                 
                 delu4(9) = delua(3,i,j,k+1) + delua(3,i,j,k-1)
                 delu4(9) = delu4(9)*del_f(KAXIS)
#endif

! compute the error
                 num = 0.
                 denom = 0.
                 
                 do kk = 1, SQNDIM
                    num = num + delu2(kk)**2
                    denom = denom + (delu3(kk) + &
                         (refine_filter*delu4(kk)))**2
                 end do
                    
                    ! mz -- compare the square of the error
                 if (denom .eq. 0.0 .AND. num .ne. 0.0) then
                    error(lb) = HUGE(1.0)
                 else if (denom .ne. 0.0) then
                    error(lb) = max(error(lb),num/denom)
                 end if
                    
              end do
           end do
        end do
           
           ! store the maximum error for the current block
        error(lb) = sqrt(error(lb))
        deallocate(delu)
        deallocate(delua)
     end if
        
  end do
     
  
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
  
  if (no_permanent_guardcells) gcell_on_cc = gcell_on_cc_backup
  !=========================================================================
  return
end subroutine gr_markRefineDerefine















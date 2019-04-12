!!****if* source/Grid/GridMain/AMR/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub. The AMReX-based Grid
!!  currently does not understand this way of implementing
!!  refinement criteria either, and uses callbacks instead,
!!  which are implemented as private functions of the Grid unit.
!!
!!  With the PARAMESH-based Grid implementation,
!!  this routine is normally called by the implementation of
!!  Grid_updateRefinement. It may also get called repeatedly
!!  during the initial construction of the Grid from
!!  Grid_initDomain.
!!
!! ARGUMENTS
!!
!!  none
!!
!! SEE ALSO
!!
!!  Grid_updateRefinement
!!  Grid_initDomain
!!  gr_expandDomain
!!
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

!!REORDER(4): solnData

subroutine Grid_markRefineDerefine()

  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Grid_iterator,  ONLY : Grid_iterator_t
  use Grid_tile,      ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_fillGuardCells,  &
                             Grid_getTileIterator, &
                             Grid_releaseTileIterator

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount,&
                        gr_enforceMaxRefinement, gr_maxRefine,&
                        gr_lrefineMaxByTime,&
                        gr_lrefineMaxRedDoByTime,&
                        gr_lrefineMaxRedDoByLogR,&
                        gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK,&
                        gr_eosModeNow, gr_eosMode, &
                        gr_meshMe, gr_meshComm
  use Simulation_data, ONLY : sim_initDens, sim_ictr,sim_jctr,&
                              sim_kctr, sim_initRad

  use gr_interface, ONLY : gr_markRefineDerefine
  use gr_parameshInterface, ONLY : gr_pmGetListOfBlocks
  use tree, ONLY : newchild, refine, derefine, stay, nodetype,&
                   lrefine,lrefine_max, parent, nchild,child


  implicit none

  include 'Flash_mpi.h'
#include "constants.h"
#include "Flash.h"
  
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc
  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref,blkCount,lb,j
  logical,save :: gcMaskArgsLogged = .FALSE.
  integer,save :: eosModeLast = 0
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask
  real, dimension(MAXBLOCKS) :: err
  integer, dimension(MAXBLOCKS) :: blkList

  real, dimension(:,:,:,:), pointer :: solnData
  real :: maxdens(MAXBLOCKS),maxdens_parent(MAXBLOCKS)
  integer :: nsend,nrecv,ierr

  integer, dimension(MAXBLOCKS) :: reqr
  integer, dimension(MAXBLOCKS*nchild) :: reqs
  integer, dimension(MPI_STATUS_SIZE,MAXBLOCKS) :: statr
  integer, dimension(MPI_STATUS_SIZE,MAXBLOCKS*nchild) :: stats

  if(gr_lrefineMaxRedDoByTime) then
     call gr_markDerefineByTime()
  end if
  
  if(gr_lrefineMaxByTime) then
     call gr_setMaxRefineByTime()
  end if

  if (gr_eosModeNow .NE. eosModeLast) then
     gcMaskArgsLogged = .FALSE.
     eosModeLast = gr_eosModeNow
  end if

  ! that are implemented in this file need values in guardcells

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  gcMask(NUNK_VARS+1:min(maskSize,NUNK_VARS+NDIM*NFACE_VARS)) = .TRUE.
!!$  gcMask(NUNK_VARS+1:maskSize) = .TRUE.


  if (.NOT.gcMaskArgsLogged) then
     call Logfile_stampVarMask(gcMask, .true., '[Grid_markRefineDerefine]', 'gcArgs')
  end if

!!$  force_consistency = .FALSE.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.,doLogMask=.NOT.gcMaskArgsLogged,&
       selectBlockType=ACTIVE_BLKS)
  gcMaskArgsLogged = .TRUE.
!!$  force_consistency = .TRUE.

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     err(:)      = 0.0
     call gr_estimateError(err, iref, ref_filter)
     call gr_markRefineDerefine(err, ref_cut, deref_cut)
  end do


!------------------------------------------------------------------------------
!
! Apply problem-specific refinement criteria.

! Dust collapse problem:  refine center of cloud.  _Don't_ refine blocks
! that are in the "fluff" (max density less than 0.5*starting density of cloud).

  call gr_markInRadius(sim_ictr, sim_jctr, sim_kctr, sim_initRad, lrefine_max)

  nullify(solnData)
  call Grid_getTileIterator(itor, ACTIVE_BLKS, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc % getDataPtr(solnData,CENTER)
     lb = tileDesc % id
     maxdens(lb) = maxval(solnData(DENS_VAR,:,:,:))
     call tileDesc % releaseDataPtr(solnData,CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call gr_pmGetListOfBlocks(ACTIVE_BLKS, blkList, blkCount)

! Communicate maxdens of parents to their leaf children.
! Maximally refined children collect messages from parents.

  maxdens_parent(:) = 0.0
  nrecv = 0
  do i = 1, blkCount
     lb = blkList(i)
     if (nodetype(lb) == LEAF .AND. lrefine(lb) == lrefine_max) then
        if(parent(1,lb).gt.-1) then
           if (parent(2,lb).ne.gr_meshMe) then
              nrecv = nrecv + 1
              call MPI_IRecv(maxdens_parent(lb),1, &
                   FLASH_REAL, &
                   parent(2,lb), &
                   lb, &
                   gr_meshComm, &
                   reqr(nrecv), &
                   ierr)
           else
              maxdens_parent(lb) = maxdens(parent(1,lb))
           end if
        end if
     end if
  end do

  ! parents send maxdens to children

  nsend = 0
  do i = 1, blkCount
     lb = blkList(i)
     if (nodetype(lb) == PARENT_BLK .AND. lrefine(lb) == lrefine_max-1) then
        do j = 1,nchild
           if(child(1,j,lb).gt.-1) then
              if (child(2,j,lb).ne.gr_meshMe) then
                 nsend = nsend + 1
                 call MPI_ISend(maxdens(lb), &
                      1, &
                      FLASH_REAL, &
                      child(2,j,lb), &  ! PE TO SEND TO
                      child(1,j,lb), &  ! THIS IS THE TAG
                      gr_meshComm, &
                      reqs(nsend), &
                      ierr)
              end if
           end if
        end do
     end if
  end do

  if (nsend.gt.0) then
     call MPI_Waitall (nsend, reqs, stats, ierr)
  end if
  if (nrecv.gt.0) then
     call MPI_Waitall (nrecv, reqr, statr, ierr)
  end if

!!  maxdens_parent(:) = 0.0  ! <-- uncomment line for previous behavior
  do i = 1, blkCount
     lb = blkList(i)
     if (nodetype(lb) == LEAF) then
        if (maxdens(lb) < 0.5*sim_initDens) then
           refine(lb)   = .false.
!           if (maxdens_parent(lb) < 0.5*sim_initDens .AND. .NOT. stay(lb)) derefine(lb)   = .true.
           if (maxdens_parent(lb) < 0.5*sim_initDens) derefine(lb)   = .true.
        endif
     else if (nodetype(lb) == PARENT_BLK) then
        if (maxdens(lb) < 0.5*sim_initDens)  refine(lb)   = .false.
     end if
  enddo

  if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()

  if(gr_enforceMaxRefinement) call gr_enforceMaxRefine(gr_maxRefine)

  if(gr_lrefineMaxRedDoByLogR) &
       call gr_unmarkRefineByLogRadius(gr_lrefineCenterI,&
       gr_lrefineCenterJ,gr_lrefineCenterK)
  

  ! When the flag arrays are passed to Paramesh for processing, only leaf
  ! blocks should be marked. - KW
  where (nodetype(:) .NE. LEAF)
     refine(:)   = .false.
     derefine(:) = .false.
  end where
  
  return
end subroutine Grid_markRefineDerefine


!!****if* source/Simulation/SimulationMain/unitTest/Multigrid_Amrex
!!
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  call Grid_unitTest(integer(in):: fileUnit,
!!                     logical(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This unit test is to test the multigrid solver of AMReX library.
!!  It uses known analytical function as well as harmonic
!!
!! ARGUMENTS
!!
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!!
!! NOTES
!!
!!***

!!------!! Do not REORDER(4): solnData

subroutine Grid_unitTest(fileUnit,perfect)

  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, GRID_PDE_BND_NEUMANN,&
!       Grid_getBlkIndexLimits, &
       Grid_solvePoisson, &
       Grid_getBlkPtr,Grid_releaseBlkPtr, &
       Grid_getDeltas, Grid_fillGuardCells, &
       Grid_getLeafIterator, Grid_releaseLeafIterator
  use gr_interface ,ONLY : gr_findMean
  use Grid_data, ONLY : gr_meshMe, gr_meshComm
  use leaf_iterator, ONLY : leaf_iterator_t
  use block_metadata, ONLY : block_metadata_t
use amrex_amr_module,     ONLY : amrex_init_from_scratch, &
                                   amrex_max_level
!  use gr_amrexLsInterface, ONLY : Grid_amrexLsSolvePoissonUnk

#include "Flash.h"
#include "constants.h"


  implicit none
#include "Flash_mpi.h"

  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors

  real,dimension(:,:,:,:),pointer :: solnData
  type(leaf_iterator_t) :: itor
  type(block_metadata_t) :: block

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(2*MDIM) :: bcTypes
  real,dimension(2,2*MDIM) :: bcValues

  logical :: gcMask(NUNK_VARS), evaluateData

  real :: poisfact
  integer,dimension(MAXBLOCKS) :: blkList
  integer :: blkCount=0,lb,i,j,k
  real:: del(MDIM)
  real meanASOL,meanPFFT
  integer nx,ny,nz
  real, allocatable, dimension(:,:,:) :: fsrc, usol
  integer blkpoints, blkpointsaux,blkCountaux
  real L2_err, L2_erraux, Linf_err, Linf_erraux, Tvol, Tvolaux, vcell
  integer :: refinelevel
  real, parameter :: tol = 1.e-6
  real, parameter :: tolInf = 4.1e-2
  real, parameter :: tol2  = 2.24e-2

  integer TA(2),count_rate,ierr
  real :: ET

  ! -------------------------------------------------------------------
  bcTypes(:)=GRID_PDE_BND_NEUMANN  
  bcTypes(3:4)=GRID_PDE_BND_DIRICHLET
  bcValues(:,:)=0.

  call mpi_barrier(gr_meshComm,ierr)
  if (gr_meshMe .eq. 0) CALL SYSTEM_CLOCK(TA(1),count_rate)  
  poisfact=1.
   call Grid_solvePoisson(NSOL_VAR, RHS_VAR, bcTypes, bcValues, poisfact)
  call mpi_barrier(gr_meshComm,ierr)
  if (gr_meshMe .eq. 0) then
     CALL SYSTEM_CLOCK(TA(2),count_rate)
     ET=REAL(TA(2)-TA(1))/count_rate
     write(*,*) ' ' 
     write(*,*) '3 PERIODIC Poisson Solver time = ',ET,' sec.'
  endif



  ! Check error in the solution:
  L2_err = 0.
  blkpoints = 0.
  Linf_err = 0.
  Tvol = 0.
  ! Get Block iterator
!   itor = block_iterator_t(LEAF)
  call Grid_getLeafIterator(itor)
  do while (itor%is_valid())
     call itor%blkMetaData(block)
     !get the index limits of the block
     blkLimits   = block%limits
     blkLimitsGC = block%limitsGC

     ! get a pointer to the current block of data
!     call Grid_getBlkPtr(block, solnData)
     call Grid_getBlkPtr(block,solnData,CENTER)

     call Grid_getDeltas(block%level,del)

     select case (NDIM)
     case(1)
        vcell = del(IAXIS)
     case(2)
        vcell = del(IAXIS)*del(JAXIS)
     case(3)
        vcell = del(IAXIS)*del(JAXIS)*del(KAXIS)
     end select

     blkCount = blkCount + 1
     blkpoints = blkpoints + &
          (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) * &
          (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) * &
          (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)

     Tvol = Tvol + vcell*real( (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) * &
          (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) * &
          (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1))

     solnData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS),DIFF_VAR)   =       &
          abs(  solnData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS),NSOL_VAR)   - &
          solnData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS),ASOL_VAR)   )

     ! L2 norm of error:
     L2_err = L2_err + sum( vcell*solnData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS),DIFF_VAR)**2.)

     ! Linf norm of error:
     Linf_err = max(Linf_err,maxval(solnData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS),DIFF_VAR) ))

     call Grid_releaseBlkPtr(block,solnData,CENTER)
     call itor%next()
  enddo
 call Grid_releaseLeafIterator(itor)
! #if defined(__GFORTRAN__) && (__GNUC__ <= 4)
!   call destroy_iterator(itor)
! #endif


  ! Sum processors Volumes
  Tvolaux = Tvol
  call MPI_Allreduce(Tvolaux, Tvol, 1, FLASH_REAL,&
       MPI_SUM, gr_meshComm, ierr)

  ! Sum processors points
  blkpointsaux = blkpoints
  call MPI_Allreduce(blkpointsaux, blkpoints, 1, FLASH_INTEGER,&
       MPI_SUM, gr_meshComm, ierr)

  ! Sum processors L2 norm of error squared
  L2_erraux = L2_err
  call MPI_Allreduce(L2_erraux, L2_err, 1, FLASH_REAL,&
       MPI_SUM, gr_meshComm, ierr)


  ! Compute L2 norm for whole domain
  L2_err = sqrt(L2_err/Tvol)

  ! Sum processors Linf norm of error
  Linf_erraux = Linf_err
  call MPI_Allreduce(Linf_erraux, Linf_err, 1, FLASH_REAL,&
       MPI_MAX, gr_meshComm, ierr)


  ! Fill GuardCells:
  gcMask = .TRUE.                         
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS,mask=gcMask)          

  ! Export to Tecplot:
  !call outtotecplot(gr_meshMe,0.0,1.,1,0,0.0,blkList,blkCount,0)


   call gr_findMean(ASOL_VAR,2,.false.,meanASOL)
   call gr_findMean(NSOL_VAR,2,.false.,meanPFFT)

  !Unit test gives a mean analytical solution of zero.  Ensure the absolute
  !value of the mean numerical solution is between 0 and the tolerance value.
  perfect = abs(meanPFFT) < MAX(TINY(1.), tol)

  if (perfect) perfect = (abs(Linf_err) .LE. tolInf)
  if (perfect) perfect = (L2_err .LE. tol2)

  ! Sum processors points
  call MPI_Allreduce(blkCount, blkCountaux, 1, FLASH_INTEGER,&
       MPI_SUM, gr_meshComm, ierr)


  if (gr_meshMe .eq. 0) then
     write(*,*) ' ' 
     write(*,'(A,2g16.8)') ' Mean Analytical, Numerical Sol=',meanASOL,meanPFFT
     write(*,'(A,1g16.8)') " ||Phi - PhiAnalytical||inf =" ,Linf_err
     write(*,'(A,1g16.8)') " ||Phi - PhiAnalytical||2   =" ,L2_err
     write(*,*) " Total Volume =",Tvol
     write(*,*) " Total Number of Leaf Blocks=", blkCountaux
     write(*,*) ' ' 
  endif
  return

end subroutine Grid_unitTest

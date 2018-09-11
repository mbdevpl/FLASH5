!!****f* source/physics/sourceTerms/Burn/Burn_finalize
!!
!! NAME
!!  
!!  Burn_finalize
!!
!!
!! SYNOPSIS
!! 
!!  call Burn_finalize()
!!
!!  
!! DESCRIPTION
!!
!!  Finalizes the Burn module.
!!
!! NOTES
!!  
!!  There is no implementation that does anything.
!!
!!  The NSE arrays used with parametricBurn are deallocated by
!!  NSE_finalize, which should be called directly frm Driver_finalizeFlash.
!!
!!***


subroutine Burn_finalize()
  use timers, ONLY : timer_burner, timer_xnet, timer_tstep, timer_nraph, timer_deriv, &
     timer_jacob, timer_solve, timer_csect, timer_scrn, timer_eoscrn
  use bn_xnetData, ONLY : xnet_myid, xnet_nproc, xnet_mythread, xnet_nthread, xnet_writeTimers
  use xnet_interface, ONLY : jacobian_finalize

  use mpi
  !$ use omp_lib

  implicit none

  integer, parameter :: ntimers = 10

  real(8), allocatable, dimension(:) :: &
     t_MPI_burner, t_MPI_xnet, t_MPI_tstep, &
     t_MPI_nraph, t_MPI_deriv, t_MPI_jacob, &
     t_MPI_solve, t_MPI_csect, t_MPI_scrn, &
     t_MPI_eoscrn

  real(8), allocatable, dimension(:,:) :: &
     t_ALL_burner, t_ALL_xnet, t_ALL_tstep, &
     t_ALL_nraph, t_ALL_deriv, t_ALL_jacob, &
     t_ALL_solve, t_ALL_csect, t_ALL_scrn, &
     t_ALL_eoscrn

  real(8), allocatable, dimension(:,:,:) :: t_ALL

  real(8), allocatable, dimension(:,:) :: &
     t_ALL_omp_min, t_ALL_omp_max, t_ALL_omp_avg

  real(8), allocatable, dimension(:) :: &
     t_ALL_min, t_ALL_max, t_ALL_avg

  integer :: ierr, i, j, k

  character(2) :: cstrLen, cntimers, cstrLenM9
  integer, parameter :: strLen = 14
  character(strLen), parameter :: tHeader(ntimers+2) = [ character(strLen) :: &
    "           MPI", &
    "           OMP", &
    "      t_burner", &
    "        t_xnet", &
    "       t_tstep", &
    "       t_nraph", &
    "       t_deriv", &
    "       t_jacob", &
    "       t_solve", &
    "       t_csect", &
    "        t_scrn", &
    "      t_eoscrn" ]

  if (xnet_writeTimers) then

     write(cstrLen,'(i2)') strLen
     write(cntimers,'(i2)') ntimers
     write(cstrLenM9,'(i2)') max(strLen-9,1)

     allocate(t_MPI_burner(0:xnet_nthread-1))
     allocate(t_MPI_xnet(0:xnet_nthread-1))
     allocate(t_MPI_tstep(0:xnet_nthread-1))
     allocate(t_MPI_nraph(0:xnet_nthread-1))
     allocate(t_MPI_deriv(0:xnet_nthread-1))
     allocate(t_MPI_jacob(0:xnet_nthread-1))
     allocate(t_MPI_solve(0:xnet_nthread-1))
     allocate(t_MPI_csect(0:xnet_nthread-1))
     allocate(t_MPI_scrn(0:xnet_nthread-1))
     allocate(t_MPI_eoscrn(0:xnet_nthread-1))

     allocate(t_ALL_burner(0:xnet_nthread-1,0:xnet_nproc-1))
     allocate(t_ALL_xnet(0:xnet_nthread-1,0:xnet_nproc-1))
     allocate(t_ALL_tstep(0:xnet_nthread-1,0:xnet_nproc-1))
     allocate(t_ALL_nraph(0:xnet_nthread-1,0:xnet_nproc-1))
     allocate(t_ALL_deriv(0:xnet_nthread-1,0:xnet_nproc-1))
     allocate(t_ALL_jacob(0:xnet_nthread-1,0:xnet_nproc-1))
     allocate(t_ALL_solve(0:xnet_nthread-1,0:xnet_nproc-1))
     allocate(t_ALL_csect(0:xnet_nthread-1,0:xnet_nproc-1))
     allocate(t_ALL_scrn(0:xnet_nthread-1,0:xnet_nproc-1))
     allocate(t_ALL_eoscrn(0:xnet_nthread-1,0:xnet_nproc-1))

     allocate(t_ALL(ntimers,0:xnet_nthread-1,0:xnet_nproc-1))
     allocate(t_ALL_omp_min(ntimers,0:xnet_nproc-1))
     allocate(t_ALL_omp_max(ntimers,0:xnet_nproc-1))
     allocate(t_ALL_omp_avg(ntimers,0:xnet_nproc-1))
     allocate(t_ALL_min(ntimers))
     allocate(t_ALL_max(ntimers))
     allocate(t_ALL_avg(ntimers))

     !$omp parallel default(shared)
!     t_MPI_burner(xnet_mythread) = timer_burner
!     t_MPI_xnet(xnet_mythread)   = timer_xnet
!     t_MPI_tstep(xnet_mythread)  = timer_tstep
!     t_MPI_nraph(xnet_mythread)  = timer_nraph
!     t_MPI_deriv(xnet_mythread)  = timer_deriv
!     t_MPI_jacob(xnet_mythread)  = timer_jacob
!     t_MPI_solve(xnet_mythread)  = timer_solve
!     t_MPI_csect(xnet_mythread)  = timer_csect
!     t_MPI_scrn(xnet_mythread)   = timer_scrn
!     t_MPI_eoscrn(xnet_mythread) = timer_eoscrn
     !$omp end parallel

     call MPI_GATHER(t_MPI_burner, xnet_nthread, MPI_REAL8, t_ALL_burner, xnet_nthread, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(t_MPI_xnet, xnet_nthread, MPI_REAL8, t_ALL_xnet, xnet_nthread, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(t_MPI_tstep, xnet_nthread, MPI_REAL8, t_ALL_tstep, xnet_nthread, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(t_MPI_nraph, xnet_nthread, MPI_REAL8, t_ALL_nraph, xnet_nthread, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(t_MPI_deriv, xnet_nthread, MPI_REAL8, t_ALL_deriv, xnet_nthread, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(t_MPI_jacob, xnet_nthread, MPI_REAL8, t_ALL_jacob, xnet_nthread, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(t_MPI_solve, xnet_nthread, MPI_REAL8, t_ALL_solve, xnet_nthread, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(t_MPI_csect, xnet_nthread, MPI_REAL8, t_ALL_csect, xnet_nthread, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(t_MPI_scrn, xnet_nthread, MPI_REAL8, t_ALL_scrn, xnet_nthread, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(t_MPI_eoscrn, xnet_nthread, MPI_REAL8, t_ALL_eoscrn, xnet_nthread, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

     if (xnet_myid == 0) then
        t_ALL(1,:,:) = t_ALL_burner
        t_ALL(2,:,:) = t_ALL_xnet
        t_ALL(3,:,:) = t_ALL_tstep
        t_ALL(4,:,:) = t_ALL_nraph
        t_ALL(5,:,:) = t_ALL_deriv
        t_ALL(6,:,:) = t_ALL_jacob
        t_ALL(7,:,:) = t_ALL_solve
        t_ALL(8,:,:) = t_ALL_csect
        t_ALL(9,:,:) = t_ALL_scrn
        t_ALL(10,:,:) = t_ALL_eoscrn
        t_ALL_omp_min(:,:) = minval(t_ALL,dim=2)
        t_ALL_omp_max(:,:) = maxval(t_ALL,dim=2)
        t_ALL_omp_avg(:,:) = sum(t_ALL,dim=2) / xnet_nthread
        t_ALL_min(:) = minval(t_ALL_omp_min,dim=2)
        t_ALL_max(:) = maxval(t_ALL_omp_max,dim=2)
        t_ALL_avg(:) = sum(t_ALL_omp_max,dim=2) / xnet_nproc
        write(*,*)
        write(*,'(2A'//trim(cstrLen)//','//trim(cntimers)//'A'//trim(cstrLen)//')') (tHeader(k),k=1,ntimers+2)
        do i = 0, xnet_nproc-1
           write(*,'(A)') repeat("-",strLen*(ntimers+2))
           do j = 0, xnet_nthread-1
              write(*,'(2I'//trim(cstrLen)//','//trim(cntimers)//'ES'//trim(cstrLen)//'.'//trim(cstrLenM9)//')') &
                i, j, (t_ALL(k,j,i),k=1,ntimers)
           end do
           write(*,'(2('//trim(cstrLen)//'x),A)') repeat("-",strLen*ntimers)
           write(*,'(I'//trim(cstrLen)//'A'//trim(cstrLen)//','//trim(cntimers)//'ES'//trim(cstrLen)//'.'//trim(cstrLenM9)//')') &
             i, 'MIN:', (t_ALL_omp_min(k,i),k=1,ntimers)
           write(*,'(I'//trim(cstrLen)//'A'//trim(cstrLen)//','//trim(cntimers)//'ES'//trim(cstrLen)//'.'//trim(cstrLenM9)//')') &
             i, 'MAX:', (t_ALL_omp_max(k,i),k=1,ntimers)
           write(*,'(I'//trim(cstrLen)//'A'//trim(cstrLen)//','//trim(cntimers)//'ES'//trim(cstrLen)//'.'//trim(cstrLenM9)//')') &
             i, 'AVG:', (t_ALL_omp_avg(k,i),k=1,ntimers)
        end do
        write(*,'('//trim(cstrLen)//'x,A)') repeat("-",strLen*(ntimers+1))
        write(*,'(2A'//trim(cstrLen)//','//trim(cntimers)//'A'//trim(cstrLen)//')') (tHeader(k),k=1,ntimers+2)
        write(*,'('//trim(cstrLen)//'x,A)') repeat("-",strLen*(ntimers+1))
        write(*,'('//trim(cstrLen)//'x,A'//trim(cstrLen)//','//trim(cntimers)//'ES'//trim(cstrLen)//'.'//trim(cstrLenM9)//')') &
          'MIN:', (t_ALL_min(k),k=1,ntimers)
        write(*,'('//trim(cstrLen)//'x,A'//trim(cstrLen)//','//trim(cntimers)//'ES'//trim(cstrLen)//'.'//trim(cstrLenM9)//')') &
          'MAX:', (t_ALL_max(k),k=1,ntimers)
        write(*,'('//trim(cstrLen)//'x,A'//trim(cstrLen)//','//trim(cntimers)//'ES'//trim(cstrLen)//'.'//trim(cstrLenM9)//')') &
          'AVG:', (t_ALL_avg(k),k=1,ntimers)
     end if

     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     deallocate(t_MPI_burner)
     deallocate(t_MPI_xnet)
     deallocate(t_MPI_tstep)
     deallocate(t_MPI_nraph)
     deallocate(t_MPI_deriv)
     deallocate(t_MPI_jacob)
     deallocate(t_MPI_solve)
     deallocate(t_MPI_csect)
     deallocate(t_MPI_scrn)
     deallocate(t_MPI_eoscrn)

     deallocate(t_ALL_burner)
     deallocate(t_ALL_xnet)
     deallocate(t_ALL_tstep)
     deallocate(t_ALL_nraph)
     deallocate(t_ALL_deriv)
     deallocate(t_ALL_jacob)
     deallocate(t_ALL_solve)
     deallocate(t_ALL_csect)
     deallocate(t_ALL_scrn)
     deallocate(t_ALL_eoscrn)

     deallocate(t_ALL)
     deallocate(t_ALL_omp_min)
     deallocate(t_ALL_omp_max)
     deallocate(t_ALL_omp_avg)
     deallocate(t_ALL_min)
     deallocate(t_ALL_max)
     deallocate(t_ALL_avg)

  end if

  call jacobian_finalize

  return

end subroutine Burn_finalize

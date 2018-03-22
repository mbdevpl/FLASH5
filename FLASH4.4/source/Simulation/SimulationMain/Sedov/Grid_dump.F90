!!****if* source/Simulation/SimulationMain/Sedov/Grid_dump
!!
!! NAME
!!  Grid_dump
!!
!! SYNOPSIS
!!
!!  call Grid_dump(integer(IN) :: var(num),
!!                 integer(IN) :: num,
!!                 integer(IN) :: blockID,
!!                 logical(IN) :: gcell)
!!
!! DESCRIPTION 
!!  
!! Dumps the variables specified in "var" to a file. Can be done from 
!! anywhere in the code, and is useful for diagnostic purposes
!! With paramesh this function doesn not work in parallel, but works
!! only with a single block
!!  
!! ARGUMENTS 
!!
!!  var :: array containing the indices of the variables to be dumped
!!  num :: number of variables being dumped.
!!  blockID :: local number of block to be dumped
!!  gcell :: indicates whether to include guardcells in the dump.
!!             
!! EXAMPLE
!!  
!!  num = 3  !dumping 3 variables
!!  var(1) = DENS_VAR
!!  var(2) = PRES_VAR
!!  var(3) = TEMP_VAR
!!  blockID = 1  ! local block number
!!  gcell = .false.
!!
!!  call Grid_dump(var, num, blockID, gcell)
!!  
!!  will dump the interior cell values of density, pressure and temperature
!!  for local block number 1.
!!
!! NOTES
!!  DENS_VAR, PRES_VAR, TEMP_VAR etc are #defined values in Flash.h
!!  indicating the index in the physical data array.
!!  The routine calling Grid_dump will need to include Flash.h 
!! 
!!***

#include "Flash.h"

subroutine Grid_dump(var,num, blockID, gcell)

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Driver_interface, ONLY : Driver_getNStep, Driver_getSimTime, Driver_getDt

#ifdef FIXEDBLOCKSIZE
  use Grid_data, ONLY : gr_ilo, gr_ihi, gr_jlo, gr_jhi, &
       gr_klo, gr_khi, gr_iloGC, gr_ihiGC, gr_jloGC, gr_jhiGC, &
       gr_kloGC, gr_khiGC
#endif

  use Simulation_data , ONLY : sim_fileUnitOutNum, sim_fileUnitOutAna
  use Simulation_data , ONLY : sim_ffNum, sim_ffAna
  use Simulation_data , ONLY : sim_globalNumProcs, sim_globalMe

#include "constants.h"

  implicit none

  integer, intent(IN) :: num, blockID
  integer, dimension(num), intent(IN) :: var
  logical, intent(IN) :: gcell
  
  integer,target :: blkLimitsGC(LOW:HIGH,MDIM)
  integer,target :: blkLimits(LOW:HIGH,MDIM)
  integer,pointer :: dumpLimits(:,:)

!!$  character(len=80) :: ffTmp
  integer,dimension(4), save :: filecount = 0
  integer :: i,j,k,count,nstep
  integer,save :: iglobalCell = 0
  real,pointer :: blkPtr(:,:,:,:)
  real :: rho
  real,allocatable,dimension(:) :: x,y,z

#ifdef FIXEDBLOCKSIZE
  integer,parameter :: bxn=GRID_IHI_GC-GRID_ILO_GC+1
  integer,parameter :: byn=GRID_JHI_GC-GRID_JLO_GC+1 
  integer,parameter :: bzn=GRID_KHI_GC-GRID_KLO_GC+1 
#else
  integer :: bxn,byn,bzn
#endif
  real :: stime, dt
  integer,save :: numCalls = 0, lastBlockID = -1

  if (lastBlockID < 0 .OR. blockID .LE. lastBlockID) then
     sim_ffNum = "sedSol-num"//char(48+filecount(4))//char(48+filecount(3))//&
          char(48+filecount(2))//char(48+filecount(1))
     sim_ffAna = "sedSol-ana"//char(48+filecount(4))//char(48+filecount(3))//&
          char(48+filecount(2))//char(48+filecount(1))
     if (sim_globalNumProcs > 1) then
99      format("p",I6.6,".out")
        write(sim_ffNum(15:),99) sim_globalMe
        write(sim_ffAna(15:),99) sim_globalMe
     end if
     numCalls = 0
     !print*,'filecount',filecount
     filecount(1) = filecount(1) + 1
     do i = 1,3
        if(filecount(i)==10)then
           filecount(i) = 0
           filecount(i+1)=filecount(i+1)+1
        end if
     end do
     close(sim_fileUnitOutNum)
     close(sim_fileUnitOutAna)
     print*,'Opening new file:',sim_ffNum
     open(sim_fileUnitOutNum,file=sim_ffNum,form='formatted')
     write(sim_fileUnitOutNum,'(a10,7a15)') &
          '#   cellno', 'x       ', 'y     ', 'dens    ', 'pres    ', 'velx    ', 'eint (spec.)'
     print*,'Opening new file:',sim_ffAna
     open(sim_fileUnitOutAna,file=sim_ffAna,form='formatted')
     write(sim_fileUnitOutAna,'(a10,7a15)') &
          '#   cellno', 'x       ', 'y     ', 'dens    ', 'pres    ', 'velx    ', 'eint density'
  end if

  if (numCalls == 0) then
     call Driver_getNStep(nstep)
     call Driver_getSimTime(stime)
     call Driver_getDt(dt)
     write(sim_fileUnitOutNum, '(1x,"#nstep, stime, dt =",i8,7(1pe15.6))') &    ! ,ADVANCE='NO') &
          nstep, stime, dt
     write(sim_fileUnitOutAna, '(1x,"#nstep, stime, dt =",i8,8(1pe15.6))') &
          nstep, stime, dt
  end if
       
  call Grid_getBlkPtr(blockID,blkPtr)
  call Grid_getBlkIndexLimits(blockID, blkLimits,blkLimitsGC)
#ifndef FIXEDBLOCKSIZE
  bxn = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  byn = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  bzn = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
#endif
  count = bxn*byn*bzn


  allocate(x(bxn))
  allocate(y(byn))
  allocate(z(bzn))
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.TRUE.,x,bxn)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.TRUE.,y,byn)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.TRUE.,z,bzn)


  if(.not. gcell) then
     dumpLimits => blkLimits
#ifdef DEBUG_SIM
#ifdef FIXEDBLOCKSIZE
     print '(8F11.6)', blkPtr(var(1), gr_ilo:gr_ihi, gr_jlo:gr_jhi, gr_klo:gr_khi)
#endif
#endif
  else
     dumpLimits = blkLimitsGC
#ifdef DEBUG_SIM
#ifdef FIXEDBLOCKSIZE
     print '(16F7.3)', blkPtr(var(1), gr_iloGc:gr_ihiGc, gr_jloGc:gr_jhiGc, gr_kloGc:gr_khiGc)
#endif
#endif
  end if

  iglobalCell = 0
     do k = dumpLimits(LOW,KAXIS), dumpLimits(HIGH,KAXIS)
        do j = dumpLimits(LOW,JAXIS), dumpLimits(HIGH,JAXIS)
           do i = dumpLimits(LOW,IAXIS), dumpLimits(HIGH,IAXIS)

              iglobalCell = iglobalCell + 1
              ! Write temperatures (and whatever else) to a file. If
              ! you add additional quantities here, you might want to
              ! also change the file header in Simulation_init.F90:
              rho = blkPtr(DENS_VAR,i,j,k)
              write(sim_fileUnitOutNum, '(i10,7(1pe15.6))') &
                   iglobalCell, x(i), y(j), &
                   blkPtr(DENS_VAR,i,j,k), &
                   blkPtr(PRES_VAR,i,j,k), &
                   blkPtr(VELX_VAR,i,j,k), &
                   blkPtr(EINT_VAR,i,j,k)
              write(sim_fileUnitOutAna, '(i10,8(1pe15.6))') &
                   iglobalCell, x(i), y(j), &
                   blkPtr(DENA_VAR,i,j,k), &
                   blkPtr(PRSA_VAR,i,j,k), &
                   blkPtr(VLXA_VAR,i,j,k), &
                   blkPtr(EINT_VAR,i,j,k)*rho
              if ((blkPtr(VLXA_VAR,i,j,k)==0.0) .AND. (blkPtr(VELX_VAR,i,j,k)==0.0)) go to 100
!!$              exit
           enddo
        end do
     end do
     100 continue



  !print '(8F11.6)', unk(var(1), gr_ilo:gr_ihi, gr_jlo:gr_jhi, gr_klo:gr_khi, blockID) 
 

  deallocate(x,y,z)
  lastBlockID = blockID
  call Grid_releaseBlkPtr(blockID,blkPtr)

  numCalls = numCalls + 1
  return
end subroutine Grid_dump

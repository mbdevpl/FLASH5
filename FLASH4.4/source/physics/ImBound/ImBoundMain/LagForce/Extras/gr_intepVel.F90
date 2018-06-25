!!****if* source/physics/ImBound/ImBoundMain/LagForce/Extras/gr_intepVel
!!
!! NAME
!!  gr_intepVel
!!
!! SYNOPSIS
!!
!!  gr_intepVel()
!!
!! DESCRIPTION
!!  Modified based on gr_sbUpdateForces by Shizhao Wang
!!  Sep 10, 2014
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

subroutine gr_intepVel(file_name, ind)
  
  use Grid_interface, ONLY : Grid_getListOfBlocks

  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, totalPart, &
                        gr_sbParticleCount, solid_body

  use gr_ptData, ONLY :  gr_ptBlkList, gr_ptBlkCount

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use gr_ptInterface, ONLY : gr_ptMove, gr_ptSetIndices, gr_ptResetIndices

  use Particles_data, ONLY : pt_indexList, pt_indexCount,pt_maxPerProc,pt_posinitialized

#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY :  gr_maxParticlesPerProc, gr_meshMe, gr_meshComm
#else
  use Grid_data, ONLY :  gr_meshMe, gr_meshComm
#endif

  implicit none
#include "Flash_mpi.h"
  character(len=*) file_name
  integer :: ind

  type(solid_body), pointer :: bodyInfo
  integer :: i,j,k,b,p, gettingFrom, blockID, recvCount
  real :: particleData(NPART_PROPS)
  integer, save, dimension(MAXBLOCKS) :: listOfBlocks
  integer, save :: count
  integer :: ct,localCount,globalCount,localPart
  logical :: moveDone,coords_in_blks=.true.


#ifdef FLASH_GRID_PARAMESH
  real :: particles2(NPART_PROPS,gr_maxParticlesPerProc)
#else
  real :: particles2(NPART_PROPS,pt_maxPerProc)
#endif


#ifdef FLASH_GRID_PARAMESH
  call Grid_getListOfBlocks(LEAF, listOfBlocks, count)
  call Grid_getListOfBlocks(LEAF, gr_ptBlkList, gr_ptBlkCount)
#else
  count=CONSTANT_ONE
  listOfBlocks(CONSTANT_ONE)=CONSTANT_ONE
#endif

  call Timers_start("InterpolateLagrangianVel")

  ! Loop over all bodies
  do b = 1, gr_sbNumBodies
     bodyInfo => gr_sbBodyInfo(b)

     call Timers_start("Create_Particles2")
     particles2(:,:) = 0.
#ifdef FLASH_GRID_PARAMESH
     globalCount = gr_maxParticlesPerProc  
     moveDone=.false.
#else
     globalCount=pt_maxPerProc 
#endif

     !! Now get all the indices into the data structure setup right
     !! for the rest of the unit
     call gr_ptSetIndices(pt_indexList, pt_indexCount)

     if (bodyInfo % myPE == bodyInfo % bodyMaster) then

        totalPart = bodyInfo % totalPart

        ct = 0 
        do i = 1, totalPart
           if (int(bodyInfo % particles(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then
              ct= ct+1
              particles2(:,ct) = bodyInfo % particles(:,i)
           endif
        enddo
        
        localCount  = ct

        ! Generate virtual Particles: This call to gr_ptMove is calles with the
        ! only purpose of generating virtual particles for surface markers on the
        ! framework of the Immersed Boundary method. This might have to change if
        ! Different types of particles are used.
#ifdef FLASH_GRID_PARAMESH
        call gr_ptMove(particles2,NPART_PROPS, localCount,globalCount, moveDone)
#else
        pt_posinitialized=.true.
        call Grid_moveParticles(particles2,NPART_PROPS, globalCount, localCount, &
             pt_indexList, pt_indexCount, coords_in_blks)
#endif
     else

        localCount = gr_sbParticleCount(b)
  
        if (localCount .gt. 0) &
        particles2(:,1:localCount) = bodyInfo%particles(:,1:localCount)

        ! Generate virtual Particles:
#ifdef FLASH_GRID_PARAMESH     
        call gr_ptMove(particles2,NPART_PROPS, localCount,globalCount, moveDone)
#else
        pt_posinitialized=.true.
        call Grid_moveParticles(particles2,NPART_PROPS, globalCount, localCount, &
             pt_indexList, pt_indexCount, coords_in_blks)
#endif
     endif
     call Timers_stop("Create_Particles2")


     call gr_ptResetIndices(pt_indexList, pt_indexCount)     

     recvCount = 0
     if (bodyInfo % myPE == bodyInfo % bodyMaster) then

        ! Find number of particles in the Master Processor:
        localPart = 0
        do j=1,size(bodyInfo%particlesPerProc,DIM=2)
           if(bodyInfo%particlesPerProc(1,j) .eq. bodyInfo%bodyMaster) then
              localPart = bodyInfo%particlesPerProc(2,j)
           endif
        enddo

        recvCount = localCount

        do i = 1, localCount

           if (int(particles2(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then
#ifdef FLASH_GRID_PARAMESH
              blockID = int(particles2(BLK_PART_PROP,i))
#else
              blockID=CONSTANT_ONE
#endif
              particleData = particles2(1:NPART_PROPS,i)

              !call Timers_start("Grid_updateSolidBodyForces")      
              call ib_getVel(b,int(particleData(GLOB_PART_PROP)),blockID, particleData)
              !call Timers_stop("Grid_updateSolidBodyForces")

              if (i .le. localPart) &
              bodyInfo % particles(1:NPART_PROPS,i) = particleData   

           end if
        end do

     else ! if (bodyInfo % myPE == bodyInfo % bodyMaster) then 

           gettingFrom = gr_sbParticleCount(b)
 
           recvCount = localCount           

           do p = 1, recvCount

              blockID = int(particles2(BLK_PART_PROP,p))
              particleData =  particles2(1:NPART_PROPS,p)

              !call Timers_start("Grid_updateSolidBodyForces") 
              call ib_getVel(b,p,blockID, particleData)
              !call Timers_stop("Grid_updateSolidBodyForces") 

              if (p .le. gettingFrom) &
              bodyInfo % particles(1:NPART_PROPS,p) = particleData 
           enddo

     end if ! if (bodyInfo % myPE == bodyInfo % bodyMaster) then

  end do  ! End Loop over Bodies

  call Timers_stop("InterpolateLagrangianVel")


  ! Output Lagragian velocity
  call sm_uvTecplot(file_name, ind) 

  return
end subroutine gr_intepVel

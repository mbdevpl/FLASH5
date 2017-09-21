!!****if* source/Simulation/SimulationMain/unitTest/Eos/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer (in) :: blockId, 
!!                   
!!
!!
!!
!! DESCRIPTION
!!   initializes a 3-d domain by setting with a density 
!!   gradient in x, temperature gradient in y, and composition gradient
!!   in z. An Eos call is made with mode = MODE_DENSITY_TEMP, and the
!!   resultant pressure and energy are saved. 
!!  
!!
!! 
!! ARGUMENTS
!!
!!  blockId -        the ID of block to put initial condtions into
!!  
!!
!!
!!***

subroutine Simulation_initBlock(solnData,block)


  use Simulation_data, ONLY : sim_xmin,sim_xmax,sim_ymin,sim_ymax,&
                              sim_zmin,sim_zmax,sim_smallx,&
                              sim_densMin, sim_tempMin, sim_xnMin, sim_presMin,&
                              sim_densMax, sim_tempMax, sim_xnMax, sim_presMax, &
                              sim_initialMass
  use Grid_interface, ONLY : Grid_getCellCoords, Grid_getGlobalIndexLimits
    
  use block_metadata, ONLY : block_metadata_t

  implicit none

#include "constants.h"
#include "Flash.h"
  
  real,dimension(:,:,:,:),pointer :: solnData
  type(block_metadata_t), intent(in) :: block

  integer :: iMax, jMax, kMax
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  
  integer :: i, j, k, n, even, put

  integer ::  npts_x, npts_y, npts_z
  real :: xpos, ypos, zpos
  real :: dx, dy, dz

  real ldensMax, ldensMin, dlogrho
  real ltempMax, ltempMin, dlogT
  real lpresMax, lpresMin, dlogP
  real lxnMax, lxnMin, dlogxn
  real :: dens,temp,pres
  real,dimension(NSPECIES) :: xn

  
  real,allocatable,dimension(:) :: xCenter,yCenter,zCenter
  logical,parameter :: gcell = .true.
  integer :: sizeX, sizeY, sizeZ
  integer,dimension(MDIM) :: pos,globalInd
  

  call Grid_getGlobalIndexLimits(globalInd)

  blkLimits = block%limits
  blkLimitsGC = block%limitsGC
  
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
  
  allocate(xCenter( blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS) ))
  allocate(yCenter( blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS) ))
  allocate(zCenter( blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) ))
  
  zCenter=0.0
  yCenter=0.0
  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, block, CENTER, gcell, zCenter, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, block, CENTER, gcell, yCenter, sizeY)

  call Grid_getCellCoords(IAXIS, block, CENTER, gcell, xCenter, sizeX)


  npts_x = globalInd(IAXIS)
  npts_y = globalInd(JAXIS)
  npts_z = globalInd(KAXIS)
  
! compute the log interval for dens, temp, and xn
  dlogrho = (log10(sim_densMax) - log10(sim_densMin)) / float(npts_x)
  dlogT   = (log10(sim_tempMax) - log10(sim_tempMin)) / float(npts_y)
  dlogP   = (log10(sim_presMax) - log10(sim_presMin)) / float(npts_y)

  if(NSPECIES > 1) then
     dlogxn  = (log10(sim_xnMax)   - log10(sim_xnMin)  ) / float(npts_z)
  end if
! compute the distance between coordinatess -- ideally xmin=0, xmax=1
  dx = (sim_xmax - sim_xmin) / float(npts_x)
  dy = (sim_ymax - sim_ymin) / float(npts_y)
  dz = (sim_zmax - sim_zmin) / float(npts_z)

! since we will be making the density, etc. equally spaced in log, compute
! the log of the extrema
  ldensMax = log10(sim_densMax)
  ldensMin = log10(sim_densMin)
  
  ltempMax = log10(sim_tempMax)
  ltempMin = log10(sim_tempMin)
  
  lpresMax = log10(sim_presMax)
  lpresMin = log10(sim_presMin)
  
  lxnMax   = log10(sim_xnMax)
  lxnMin   = log10(sim_xnMin)

! now fill the master arrays
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)

     zpos = (zCenter(k) - sim_zmin) / dz
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        
        ypos = (yCenter(j) - sim_ymin) / dy
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           xpos = (xCenter(i) - sim_xmin) / dx
           ! compute the density, temperature, and composition
           dens = 10.e0**(ldensMin + xpos*dlogrho)
           temp = 10.e0**(ltempMin + ypos*dlogT)
           pres = 10.e0**(lpresMin + ypos*dlogP)
           solnData(DENS_VAR,i,j,k)=dens
           solnData(CTMP_VAR,i,j,k)=temp
           solnData(TEMP_VAR,i,j,k)=temp
           solnData(CPRS_VAR,i,j,k)=pres
           
#if NSPECIES > 0
           ! for the composition, divide the material according to the input parameter
           !  sim_initialMass = -1 to have a gradient divided between species 1 and NSPECIES
           !  sim_initialMass = 0 to divide evenly
           !  sim_initialMass = i to put all mass in the ith species
           
           if(NSPECIES==1) then
              xn = 1.0
           else
              select case (sim_initialMass)
              case (0)              !! divide evenly
                 do even=1, NSPECIES
                    xn(even) = 1.0 / NSPECIES
                 enddo
              case (-1)            !! anshu's original
                 xn(:) = sim_smallx
                 xn(1) =  10.e0**(lxnMin + zpos*dlogxn)
                 xn(NSPECIES) = 1.e0 - xn(1)
              case default          !!  put in ith variable
                 xn(:) = sim_smallx
                 xn(sim_initialMass) = 1.0
              end select
              
              
           end if
           !! Now put the data in UNK
           do put=1,NSPECIES
              solnData(SPECIES_BEGIN+put-1,i,j,k)=xn(put)
           enddo
#endif
           
        enddo
     enddo
  enddo
  
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
  
  return
end subroutine Simulation_initBlock





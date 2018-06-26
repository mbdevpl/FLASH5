!
!
!
!**

Subroutine ib_mapParticles(ib)

  use gr_sbData, ONLY : gr_sbBodyInfo, NodesPerElem

  use ib_interface, ONLY : ib_mapParticlesElem

  use ImBound_data, only : omg,freq_t,Ro,ao,ib_dt

  use Driver_data, only : dr_simTime

  implicit none
#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: ib


  integer :: i,ii,p,NumPart,nXi,nEta,ptelem,maxptelem,vert_elem

  integer :: lb,blockID,aelem(NodesPerElem+CONSTANT_ONE)

  real, allocatable, dimension(:) :: xpos,ypos,zpos,xvel,yvel,zvel,areai
  real, allocatable, dimension(:) :: xacc,yacc,zacc,xnrm,ynrm,znrm

  real, dimension(NodesPerElem) :: xi,yi,zi,ui,vi,wi,udi,vdi,wdi,nxLi,nyLi,nzLi

  real :: tita

  real :: nrm

  ! Here do surface Kinematics, for aero nodes:
  ! Do surf kinematics, for vertex points and internal markers:
  ! The general version of this part (vertex points) will be defined in Create Particles.

#if NDIM == 2
  tita  = omg*dr_simTime
  do i=1,gr_sbBodyInfo(ib)%NumVertices


    ! Positions of markers, normals and arclength var
    gr_sbBodyInfo(ib)%xb(i) =  gr_sbBodyInfo(ib)%xo +  gr_sbBodyInfo(ib)%xbus(i)*cos(tita) &
                              -gr_sbBodyInfo(ib)%ybus(i)*sin(tita)    
    gr_sbBodyInfo(ib)%yb(i) =  gr_sbBodyInfo(ib)%yo +  gr_sbBodyInfo(ib)%xbus(i)*sin(tita) &
                             + gr_sbBodyInfo(ib)%ybus(i)*cos(tita) +ao*sin(2.*PI*freq_t*dr_simTime)

    
    gr_sbBodyInfo(ib)%nxL(i) =  (gr_sbBodyInfo(ib)%xb(i)-gr_sbBodyInfo(ib)%xo)/Ro
    gr_sbBodyInfo(ib)%nyL(i) =  (gr_sbBodyInfo(ib)%yb(i)-gr_sbBodyInfo(ib)%yo)/Ro

    ! Set velocities and accelerations
    gr_sbBodyInfo(ib)%ubd(i) = -omg*gr_sbBodyInfo(ib)%xbus(i)*sin(tita) - &
                                omg*gr_sbBodyInfo(ib)%ybus(i)*cos(tita) 
    gr_sbBodyInfo(ib)%vbd(i) =  omg*gr_sbBodyInfo(ib)%xbus(i)*cos(tita) - &
                                omg*gr_sbBodyInfo(ib)%ybus(i)*sin(tita) + &
                                (2.*PI*freq_t)*ao*cos(2.*PI*freq_t*dr_simTime)
    gr_sbBodyInfo(ib)%ubdd(i)= -omg**2.*gr_sbBodyInfo(ib)%xbus(i)*cos(tita) + &
                                omg**2.*gr_sbBodyInfo(ib)%ybus(i)*sin(tita)
    gr_sbBodyInfo(ib)%vbdd(i)= -omg**2.*gr_sbBodyInfo(ib)%xbus(i)*sin(tita) - &
                                omg**2.*gr_sbBodyInfo(ib)%ybus(i)*cos(tita) - &
                                (2.*PI*freq_t)**2.*ao*sin(2.*PI*freq_t*dr_simTime)

  enddo 

#elif NDIM == 3

  do i=1,gr_sbBodyInfo(ib)%NumVertices

    ! Positions of markers, normals and arclength var
    gr_sbBodyInfo(ib)%xb(i) =  gr_sbBodyInfo(ib)%xo + gr_sbBodyInfo(ib)%xbus(i)   
    gr_sbBodyInfo(ib)%yb(i) =  gr_sbBodyInfo(ib)%yo + gr_sbBodyInfo(ib)%ybus(i)
    gr_sbBodyInfo(ib)%zb(i) =  gr_sbBodyInfo(ib)%zo + gr_sbBodyInfo(ib)%zbus(i)
    
    gr_sbBodyInfo(ib)%nxL(i) =  (gr_sbBodyInfo(ib)%xb(i)-gr_sbBodyInfo(ib)%xo)
    gr_sbBodyInfo(ib)%nyL(i) =  (gr_sbBodyInfo(ib)%yb(i)-gr_sbBodyInfo(ib)%yo)
    gr_sbBodyInfo(ib)%nzL(i) =  (gr_sbBodyInfo(ib)%zb(i)-gr_sbBodyInfo(ib)%zo)

    nrm = sqrt(gr_sbBodyInfo(ib)%nxL(i)**2 + gr_sbBodyInfo(ib)%nyL(i)**2 + gr_sbBodyInfo(ib)%nzL(i)**2)
    
    gr_sbBodyInfo(ib)%nxL(i) = gr_sbBodyInfo(ib)%nxL(i)/nrm
    gr_sbBodyInfo(ib)%nyL(i) = gr_sbBodyInfo(ib)%nyL(i)/nrm
    gr_sbBodyInfo(ib)%nzL(i) = gr_sbBodyInfo(ib)%nzL(i)/nrm

    ! Set velocities and accelerations
    gr_sbBodyInfo(ib)%ubd(i) = 0.
    gr_sbBodyInfo(ib)%vbd(i) = 0.
    gr_sbBodyInfo(ib)%wbd(i) = 0.
    gr_sbBodyInfo(ib)%ubdd(i)= 0.
    gr_sbBodyInfo(ib)%vbdd(i)= 0.
    gr_sbBodyInfo(ib)%wbdd(i)= 0.

  enddo 

#endif


  ! Maximum number of particles per element:
  maxptelem = maxval(gr_sbBodyInfo(ib)%sbPtNumElem(:))

  ! Allocate vars:
  allocate(xpos(maxptelem),ypos(maxptelem),zpos(maxptelem))
  allocate(xvel(maxptelem),yvel(maxptelem),zvel(maxptelem),areai(maxptelem))
  allocate(xacc(maxptelem),yacc(maxptelem),zacc(maxptelem))
  allocate(xnrm(maxptelem),ynrm(maxptelem),znrm(maxptelem))
 
  numPart = 0
  do i=1,gr_sbBodyInfo(ib)%NumAelem

     aelem(:) = gr_sbBodyInfo(ib)%AELEM(:,i)

     vert_elem = aelem(1)

     ! Elements nodes Positions, Velocities, Accelerations and Normals
     xi(1:vert_elem)   = gr_sbBodyInfo(ib)%xb(aelem(2:vert_elem+1))
     yi(1:vert_elem)   = gr_sbBodyInfo(ib)%yb(aelem(2:vert_elem+1))
     ui(1:vert_elem)   = gr_sbBodyInfo(ib)%ubd(aelem(2:vert_elem+1))
     vi(1:vert_elem)   = gr_sbBodyInfo(ib)%vbd(aelem(2:vert_elem+1))
     udi(1:vert_elem)  = gr_sbBodyInfo(ib)%ubdd(aelem(2:vert_elem+1))
     vdi(1:vert_elem)  = gr_sbBodyInfo(ib)%vbdd(aelem(2:vert_elem+1))
     nxLi(1:vert_elem) = gr_sbBodyInfo(ib)%nxL(aelem(2:vert_elem+1))
     nyLi(1:vert_elem) = gr_sbBodyInfo(ib)%nyL(aelem(2:vert_elem+1))

#if NDIM == 3
     zi(1:vert_elem)   = gr_sbBodyInfo(ib)%zb(aelem(2:vert_elem+1))
     wi(1:vert_elem)   = gr_sbBodyInfo(ib)%wbd(aelem(2:vert_elem+1))
     wdi(1:vert_elem)  = gr_sbBodyInfo(ib)%wbdd(aelem(2:vert_elem+1))
     nzLi(1:vert_elem) = gr_sbBodyInfo(ib)%nzL(aelem(2:vert_elem+1))
#endif 

     nXi    = gr_sbBodyInfo(ib)%sbPtNumXi(i) 
     nEta   = gr_sbBodyInfo(ib)%sbPtNumeta(i)
     ptelem = gr_sbBodyInfo(ib)%sbPtNumElem(i)
 
     ! Get elements Particle data:
     call ib_mapParticlesElem(gr_sbBodyInfo(ib)%AELTYPE(i),vert_elem,nXi,nEta,ptelem, &
                              xi(1:vert_elem),yi(1:vert_elem),zi(1:vert_elem), &
                              ui(1:vert_elem),vi(1:vert_elem),wi(1:vert_elem), &
                              udi(1:vert_elem),vdi(1:vert_elem),wdi(1:vert_elem), &
                              nxLi(1:vert_elem),nyLi(1:vert_elem),nzLi(1:vert_elem), &
                              xpos(1:ptelem),ypos(1:ptelem),zpos(1:ptelem), &
                              xvel(1:ptelem),yvel(1:ptelem),zvel(1:ptelem), &
                              xacc(1:ptelem),yacc(1:ptelem),zacc(1:ptelem), &
                              xnrm(1:ptelem),ynrm(1:ptelem),znrm(1:ptelem), &
                              areai(1:ptelem))


     ! Load elements Particles data into particles array:
     gr_sbBodyInfo(ib) % particles(BLK_PART_PROP,numPart+1:numPart+ptelem)  = UNKNOWN
     gr_sbBodyInfo(ib) % particles(POSX_PART_PROP,numPart+1:numPart+ptelem) = xpos(1:ptelem)
     gr_sbBodyInfo(ib) % particles(POSY_PART_PROP,numPart+1:numPart+ptelem) = ypos(1:ptelem)
     gr_sbBodyInfo(ib) % particles(POSZ_PART_PROP,numPart+1:numPart+ptelem) = zpos(1:ptelem)
     gr_sbBodyInfo(ib) % particles(VELX_PART_PROP,numPart+1:numPart+ptelem) = xvel(1:ptelem)
     gr_sbBodyInfo(ib) % particles(VELY_PART_PROP,numPart+1:numPart+ptelem) = yvel(1:ptelem)
     gr_sbBodyInfo(ib) % particles(VELZ_PART_PROP,numPart+1:numPart+ptelem) = zvel(1:ptelem)
     gr_sbBodyInfo(ib) % particles(ACCX_PART_PROP,numPart+1:numPart+ptelem) = xacc(1:ptelem)
     gr_sbBodyInfo(ib) % particles(ACCY_PART_PROP,numPart+1:numPart+ptelem) = yacc(1:ptelem)
     gr_sbBodyInfo(ib) % particles(ACCZ_PART_PROP,numPart+1:numPart+ptelem) = zacc(1:ptelem)
     gr_sbBodyInfo(ib) % particles(NMLX_PART_PROP,numPart+1:numPart+ptelem) = xnrm(1:ptelem)
     gr_sbBodyInfo(ib) % particles(NMLY_PART_PROP,numPart+1:numPart+ptelem) = ynrm(1:ptelem)
     gr_sbBodyInfo(ib) % particles(NMLZ_PART_PROP,numPart+1:numPart+ptelem) = znrm(1:ptelem)


     gr_sbBodyInfo(ib) % particles(FUL_PART_PROP,numPart+1:numPart+ptelem)  = 0.
     gr_sbBodyInfo(ib) % particles(FVL_PART_PROP,numPart+1:numPart+ptelem)  = 0. 
     gr_sbBodyInfo(ib) % particles(FWL_PART_PROP,numPart+1:numPart+ptelem)  = 0.
     gr_sbBodyInfo(ib) % particles(TAG_PART_PROP,numPart+1:numPart+ptelem)  = 1.
     gr_sbBodyInfo(ib) % particles(AREA_PART_PROP,numPart+1:numPart+ptelem) = areai(1:ptelem)

     do p = numPart+1,numPart+ptelem
        gr_sbBodyInfo(ib) % particles(GLOB_PART_PROP,p) = p !local particle counter in solid body.
     enddo

     numPart = numPart + ptelem

  enddo


  ! DeAllocate vars:
  deallocate(xpos,ypos,zpos)
  deallocate(xvel,yvel,zvel,areai)
  deallocate(xacc,yacc,zacc)
  deallocate(xnrm,ynrm,znrm)
 
  return

end Subroutine ib_mapParticles

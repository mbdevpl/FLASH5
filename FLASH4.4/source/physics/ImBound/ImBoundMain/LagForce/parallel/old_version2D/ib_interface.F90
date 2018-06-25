module ib_interface

#include "Flash.h"
#include "constants.h"

  interface 
     subroutine ib_forcing(blockID, particleData)
       
       use Grid_data, ONLY: gr_meshMe
       
       use ImBound_data
                   
       use Grid_interface, ONLY : Grid_getDeltas,          &
            Grid_getBlkCenterCoords, &
            Grid_getBlkPhysicalSize, &
            Grid_getBlkPtr,          &
            Grid_releaseBlkPtr
       
       use Driver_interface, ONLY : Driver_abortFlash
       
       implicit none
       
       !! ---- Argument List ----------------------------------
       integer, INTENT(IN) :: blockID
       real, INTENT(INOUT) :: particleData(NPART_PROPS)
       !! -----------------------------------------------------
     end subroutine ib_forcing
  end interface
  
  interface 
     subroutine ib_stencils(ng,nx,ny,xb,yb,dsx,dsy,   & 
          ielem,jelem,gridflag,                           &
          del,coord,bsize)
       
       use ImBound_data , ONLY : ib_stencil, ib_alphax, ib_alphay
       use gr_sbData, ONLY : gr_sbBodyInfo 
       
       !! ---- Argument List ----------------------------------
       integer, INTENT(IN) :: ng,nx,ny,gridflag
       real, INTENT(IN) :: xb,yb
       real, INTENT(INOUT) :: dsx,dsy
       integer, INTENT(INOUT) :: ielem(ib_stencil),jelem(ib_stencil)
       real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)
         !! -----------------------------------------------------
     end subroutine ib_stencils
  end interface
  
  interface 
     subroutine ib_interpLpoints(xb,yb,dsx,dsy,      &
          ielem,jelem,phile,zL,VarO,gridflag,            &
          nx1,ny1,nz1,ng,nx,ny,del,coord,bsize)
       
       use ImBound_data , only :ib_stencil,ib_interp,ib_npol
       implicit none

       integer, INTENT(IN) :: gridflag,nx1,ny1,nz1,ng,nx,ny
       real, INTENT(IN) :: xb,yb
       real, INTENT(IN) :: dsx ,dsy
       integer, INTENT(IN) :: ielem(ib_stencil),jelem(ib_stencil)
       real, INTENT(INOUT) :: phile(ib_stencil),zL
       real, INTENT(IN) :: VarO(nx1,ny1,nz1)
       real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)
     end subroutine ib_interpLpoints
  end interface
  
  
  interface
     subroutine ib_extrapEpoints(xb,yb,sb,           &
          ielem,jelem,phile,zL,VarO,      &
          nx1,ny1,nz1,del,coord,bsize)
       
       
       use ImBound_data , only : ib_stencil
       implicit none       
!!! Argument list
       integer, INTENT(IN) :: nx1,ny1,nz1
       real, INTENT(IN) :: xb,yb,sb
       integer, INTENT(IN) :: ielem(ib_stencil),jelem(ib_stencil)
       real, INTENT(INOUT) :: phile(ib_stencil)
       real, INTENT(IN) :: zL
       real, INTENT(INOUT) :: VarO(nx1,ny1,nz1)
       real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)
     end subroutine ib_extrapEpoints
  end interface



  interface

     Subroutine ib_countParticles(ib,numPart)
       implicit none
       integer, intent(IN)  :: ib
       integer, intent(OUT) :: numPart
     end Subroutine ib_countParticles

  end interface


  interface

     Subroutine ib_countParticlesElem(eltype,vert_elem,xi,yi,zi,Dmin,nXi,nEta,ptelem)
       implicit none
       integer, intent(IN)  :: eltype,vert_elem
       real, intent(IN)     :: xi(vert_elem),yi(vert_elem),zi(vert_elem),Dmin
       integer, intent(OUT) :: nXi,nEta,ptelem
     end Subroutine ib_countParticlesElem

  end interface


  interface

     Subroutine ib_mapParticles(ib)
       implicit none
       integer, intent(IN) :: ib
     end Subroutine ib_mapParticles

  end interface
  

  interface

     Subroutine ib_mapParticlesElem(eltype,vert_elem,nXi,nEta,ptelem,&
                              xi,yi,zi,ui,vi,wi,udi,vdi,wdi,     &
                              nxLi,nyLi,nzLi,                    &
                              xpos,ypos,zpos,xvel,yvel,zvel,     &
                              xacc,yacc,zacc,xnrm,ynrm,znrm,areai)
       implicit none
       integer, intent(IN)  :: eltype,vert_elem,nXi,nEta,ptelem
       real, intent(IN)     :: xi(vert_elem),yi(vert_elem),zi(vert_elem), &
                          ui(vert_elem),vi(vert_elem),wi(vert_elem), &
                          udi(vert_elem),vdi(vert_elem),wdi(vert_elem), &
                          nxLi(vert_elem),nyLi(vert_elem),nzLi(vert_elem)
       real, intent(OUT) :: xpos(ptelem),ypos(ptelem),zpos(ptelem), &
                          xvel(ptelem),yvel(ptelem),zvel(ptelem), &
                          xacc(ptelem),yacc(ptelem),zacc(ptelem), &
                          xnrm(ptelem),ynrm(ptelem),znrm(ptelem), &
                          areai(ptelem)
     end Subroutine ib_mapParticlesElem

  end interface

  interface

     subroutine ib_distributedForces(blockID, particleData, vortx, vorty, vortz)
     implicit none
     integer, intent(IN) :: blockID
     real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vortx 
     real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vorty
     real, intent(IN), dimension(GRID_IHI_GC*K1D+1,GRID_JHI_GC*K2D+1,GRID_KHI_GC*K3D+1) :: vortz
     real, intent(INOUT) :: particleData(NPART_PROPS)
     end subroutine

  end interface

  interface

     subroutine ib_buildABLan(stencil,npol,dsx,dsy,xp,yp, &
                              x,y,interp,A,B,buildflag)
       implicit none
       integer, INTENT(IN) :: stencil,npol,interp
       integer, INTENT(OUT) :: buildflag
       real, INTENT(IN) :: dsx,dsy,xp,yp
       real, INTENT(IN) :: x(stencil),y(stencil)
       real, INTENT(INOUT) :: A(npol,npol),B(npol,stencil)
     end subroutine ib_buildABLan

  end interface


  interface

     subroutine ib_weightfunc(dsx,dsy,xp,yp,x,y,wtype,wi,dwidx,dwidy)
       implicit none
       integer, INTENT(IN) ::  wtype
       real, INTENT(IN) :: dsx,dsy,xp,yp,x,y
       real, INTENT(OUT) :: wi,dwidx,dwidy
     end subroutine ib_weightfunc

  end interface


end Module ib_interface


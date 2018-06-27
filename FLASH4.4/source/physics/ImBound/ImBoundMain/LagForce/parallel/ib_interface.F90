module ib_interface

#include "Flash.h"
#include "constants.h"


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
     end subroutine ib_distributedForces

  end interface



  interface

     subroutine ib_stencils(xp,np,gridfl,del,coord,bsize,   & 
                                   ielem,hl,forcflag)
       use ImBound_data , ONLY : ib_stencil
       implicit none
       integer, INTENT(IN) :: gridfl(MDIM),forcflag
       real, INTENT(IN) :: xp(MDIM),np(MDIM)
       integer, INTENT(INOUT) :: ielem(ib_stencil,MDIM)
       real, INTENT(OUT) :: hl
       real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)

     end subroutine ib_stencils

  end interface


  interface

     subroutine ib_forcing(ibd,p,blockID, particleData)

       implicit none
       !! ---- Argument List ----------------------------------
       integer, INTENT(IN) :: blockID,ibd,p
       real, INTENT(INOUT) :: particleData(NPART_PROPS)
       !! -----------------------------------------------------

     end subroutine ib_forcing

  end interface

  interface

     subroutine ib_interpLpoints(xp,gridfl,          &
                    del,coord,bsize,ielem,phile,     &
                    zL,forcflag,blockID,faceind)
     use ImBound_data , only :ib_stencil
     implicit none
     integer, INTENT(IN) :: gridfl(MDIM),forcflag,blockID,faceind
     real, INTENT(IN)    :: xp(MDIM)
     integer, INTENT(IN) :: ielem(ib_stencil,MDIM)
     real, INTENT(INOUT) :: phile(ib_stencil,NDIM+1),zL
     real, INTENT(IN)    :: del(MDIM),coord(MDIM),bsize(MDIM)
   end subroutine ib_interpLpoints

  end interface

  interface

     subroutine ib_extrapEpoints(xp,sb,hl,del,ielem,phile,zL,blockID,faceind)
       use ImBound_data , only : ib_stencil
       implicit none  
       ! Argument list
       integer, INTENT(IN) :: blockID,faceind
       real, INTENT(IN) :: xp(MDIM),sb,hl
       integer, INTENT(IN) :: ielem(ib_stencil,MDIM)
       real, INTENT(IN) :: zL
       real, INTENT(IN) :: del(MDIM)
       real, INTENT(INOUT) :: phile(ib_stencil,NDIM+1)
     end subroutine ib_extrapEpoints

  end interface

  interface

     subroutine ib_getInterpFunc(xp,xyz_stencil,del,derivflag,phile)

       use ImBound_data , only :ib_stencil
       implicit none
       real, intent(IN) :: xp(MDIM),xyz_stencil(ib_stencil,MDIM),del(MDIM)
       integer, intent(IN) :: derivflag
       real, intent(OUT):: phile(ib_stencil,NDIM+1)

     end subroutine ib_getInterpFunc

  end interface

  interface

     subroutine ib_buildABLan(stencil,interp,npol,nderiv,dsx,dsy,dsz,xp,yp,zp, &
                              x,y,z,A,B,p,buildflag,derivflag)
       implicit none
       integer, INTENT(IN) :: stencil,interp,npol,nderiv,derivflag
       integer, INTENT(OUT) :: buildflag
       real, INTENT(IN) :: dsx,dsy,dsz,xp,yp,zp
       real, INTENT(IN) :: x(stencil),y(stencil),z(stencil)
       real, INTENT(INOUT) :: A(npol,npol,nderiv),B(npol,stencil,nderiv),p(npol,nderiv)
     end subroutine ib_buildABLan

  end interface



  interface

     subroutine ib_weightfunc(dsx,dsy,dsz,xp,yp,zp,x,y,z,wtype,wi,dwidx,dwidy,dwidz)
       implicit none
       integer, INTENT(IN) ::  wtype
       real, INTENT(IN) :: dsx,dsy,dsz,xp,yp,zp,x,y,z
       real, INTENT(OUT) :: wi,dwidx,dwidy,dwidz
     end subroutine ib_weightfunc

  end interface


  interface
    subroutine ib_ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n)
      REAL d,a(np,np)
    end subroutine ib_ludcmp
  end interface

  interface
     subroutine ib_lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
    end subroutine ib_lubksb
  end interface

  interface
     subroutine ib_forceInsideBody()
       implicit none
     end subroutine ib_forceInsideBody
  end interface

end Module ib_interface


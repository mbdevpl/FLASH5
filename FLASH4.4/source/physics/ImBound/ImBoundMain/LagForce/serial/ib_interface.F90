

Module ib_interface

#include "Flash.h"
#include "constants.h"

  interface

     SUBROUTINE ib_forcing(blockCount, blockList, dt)
       implicit none
       !! ---- Argument List ----------------------------------
       integer, INTENT(IN) :: blockCount
       integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
       real, INTENT(IN) :: dt
       !! -----------------------------------------------------
     end SUBROUTINE ib_forcing

  end interface


  interface

     subroutine ib_stencils(ng,nx,ny,nbd,xb,yb,sb,dsx,dsy, & 
          ielem,jelem,flago,flagi,nxp,nyp,lb,gridflag,&
          np,lpindex,nx1,ny1,nz1,del,coord,bsize)
       use ImBound_data , only : ib_nmaxa, ib_stencil
       implicit none
#include "Flash.h"
#include "constants.h"
       integer, INTENT(IN) :: ng,nx,ny,nbd,lb,gridflag,nx1,ny1,nz1
       integer, INTENT(INOUT) :: np
       real, INTENT(IN) :: xb(ib_nmaxa),yb(ib_nmaxa),sb(ib_nmaxa),nxp(ib_nmaxa),nyp(ib_nmaxa)
       real, INTENT(INOUT) :: dsx(ib_nmaxa),dsy(ib_nmaxa)
       integer, INTENT(INOUT) :: ielem(ib_stencil,ib_nmaxa),jelem(ib_stencil,ib_nmaxa)
       integer, INTENT(INOUT) :: flago(nx1,ny1,nz1),flagi(nx1,ny1,nz1),lpindex(ib_nmaxa)
       real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)
    end subroutine ib_stencils

  end interface



  interface

     subroutine ib_interpLpoints(nbd,xb,yb,sb,dsx,dsy, &
         ielem,jelem,phile,zL,VarO,lb,gridflag,np,lpindex, &
         nx1,ny1,nz1,ng,nx,ny,del,coord,bsize)
       use ImBound_data , only : ib_nmaxa,ib_stencil,ib_interp,ib_npol
       implicit none
#include "Flash.h"
#include "constants.h"
       integer, INTENT(IN) :: nbd,lb,gridflag,np,nx1,ny1,nz1,ng,nx,ny
       real, INTENT(IN) :: xb(ib_nmaxa),yb(ib_nmaxa),sb(ib_nmaxa)
       real, INTENT(IN) :: dsx(ib_nmaxa),dsy(ib_nmaxa)
       integer, INTENT(IN) :: ielem(ib_stencil,ib_nmaxa),jelem(ib_stencil,ib_nmaxa),lpindex(ib_nmaxa)
       real, INTENT(INOUT) :: phile(ib_stencil,ib_nmaxa),zL(ib_nmaxa)
       real, INTENT(IN) :: VarO(nx1,ny1,nz1)
       real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)
     end subroutine ib_interpLpoints

  end interface

  interface

      subroutine ib_extrapEpoints(nbd,xb,yb,sb,     &
                          ielem,jelem,phile,zL,VarO,&
                          lb,np,lpindex,nx1,ny1,nz1,&
                          del,coord,bsize)
      use ImBound_data , only : ib_ABODY,ib_nmaxa,ib_stencil
      implicit none
#include "Flash.h"
#include "constants.h"
      integer, INTENT(IN) :: nbd,lb,np,nx1,ny1,nz1
      real, INTENT(IN) :: xb(ib_nmaxa),yb(ib_nmaxa),sb(ib_nmaxa)
      integer, INTENT(IN) :: ielem(ib_stencil,ib_nmaxa),jelem(ib_stencil,ib_nmaxa),lpindex(ib_nmaxa)
      real, INTENT(INOUT) :: phile(ib_stencil,ib_nmaxa)
      real, INTENT(IN) :: zL(ib_nmaxa)
      real, INTENT(INOUT) :: VarO(nx1,ny1,nz1)
      real, INTENT(IN) :: del(MDIM),coord(MDIM),bsize(MDIM)
    end subroutine ib_extrapEpoints
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


  interface

     SUBROUTINE ib_CalcForce(blockCount, blockList, dt)
       implicit none
       !! ---- Argument List ----------------------------------
       integer, INTENT(IN) :: blockCount
       integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
       real, INTENT(IN) :: dt
       !! -----------------------------------------------------
     END SUBROUTINE ib_CalcForce

  end interface


  interface


     subroutine ib_forces2D(blockCount, blockList, dt)
       implicit none
       !! ---- Argument List ----------------------------------
       integer, INTENT(IN) :: blockCount
       integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
       real, INTENT(IN) :: dt
       !! -----------------------------------------------------
     end subroutine ib_forces2D
     
  end interface
  
  interface 

     subroutine ib_calcdistforce2D(nu,ng,nx,ny,nbd,nnoda,xb,yb,sb, &
                                   nxp,nyp,xbe,ybe,                &
                                   stencil,zL,VarO,presflag,       &
                                   ubdd,vbdd,lb,nx1,ny1)
       use imBound_data, ONLY   : ib_nmaxa
       implicit none
       ! --- ARGUMENTS -----------------------------------------
       integer, INTENT(IN) :: ng,lb
       integer, INTENT(IN) :: nx,ny,nbd,nnoda,stencil,presflag,nx1,ny1
       real  , INTENT(IN)  :: nu,xb(ib_nmaxa),yb(ib_nmaxa),sb(ib_nmaxa)
       real  , INTENT(IN)  ::  nxp(ib_nmaxa),nyp(ib_nmaxa)
       real  , INTENT(INOUT)  ::  xbe(ib_nmaxa),ybe(ib_nmaxa)
       real  , INTENT(INOUT)  ::  zL(ib_nmaxa)
       real  , INTENT(IN) :: ubdd(ib_nmaxa),vbdd(ib_nmaxa)
       real  ,INTENT(IN)   :: VarO(nx1,ny1)
       ! -------------------------------------------------------
     end subroutine ib_calcdistforce2D
     
  end interface

  end Module ib_interface
  

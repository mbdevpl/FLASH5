module sm_element_interface


#include "Flash.h"
#include "SolidMechanics.h"

  interface
     subroutine el01_ShapeFunc(NN,Nr,r)
       ! Compute shape functions for 2 node line
       implicit none
       real, intent(in)  :: r
       integer, parameter  :: n = TWO_NODES
       real, intent(out) :: NN(n), Nr(n)
     end subroutine el01_ShapeFunc
  end interface


  interface
     subroutine el02_ShapeFunc(NN,Nr,Ns,r,s)
       ! Compute shape functions for 3 node triangle
       implicit none
       real, intent(in)  :: r,s
       integer, parameter  :: n = THREE_NODES
       real, intent(out) :: NN(n), Nr(n), Ns(n)
     end subroutine el02_ShapeFunc
  end interface


  interface
     subroutine el05_ShapeFunc(NN,Nr,Ns,Nt,r,s,t)
       ! Compute shape functions for a linear 8 node hexahedron
       implicit none
       real, intent(in)  :: r,s,t
       integer, parameter  :: n = EIGHT_NODES
       real, intent(out) :: NN(n), Nr(n), Ns(n), Nt(n)
     end subroutine el05_ShapeFunc
  end interface


  interface
     subroutine el09_ShapeFunc(NN,Nr,Ns,r,s)
       ! Compute shape functions for 6 node triangle
       implicit none

       real, intent(in)  :: r,s
       integer, parameter  :: n = SIX_NODES
       real, intent(out) :: NN(n), Nr(n), Ns(n)
     end subroutine el09_ShapeFunc
  end interface


  interface
     subroutine el10_ShapeFunc(NN,Nr,Ns,r,s)
       ! Compute shape functions for 9 node quad
       implicit none
       real, intent(in)  :: r,s
       integer,parameter   :: n = NINE_NODES
       real, intent(out) :: NN(n), Nr(n), Ns(n)
     end subroutine el10_ShapeFunc
  end interface

  interface
     subroutine el11_ShapeFunc(NN,Nr,Ns,Nt,r,s,t)
       ! Compute shape functions for a 10-node second order tetrahedron 
       !(4 nodes associated with the vertices and 6 with the edges)
       implicit none
       real, intent(in)  :: r,s,t
       integer,parameter   :: n = TEN_NODES
       real, intent(out) :: NN(n), Nr(n), Ns(n), Nt(n)
     end subroutine el11_ShapeFunc
  end  interface

  interface
     subroutine el12_ShapeFunc(NN,Nr,Ns,Nt,r,s,t)
       ! Compute shape functions for a quadratic 27 node hexahedron
       implicit none
       real, intent(in)  :: r,s,t
       integer,parameter   :: n = TWENTYSEVEN_NODES
       real, intent(out) :: NN(n), Nr(n), Ns(n), Nt(n)
     end subroutine el12_ShapeFunc
  end interface


  interface
     subroutine ShapeFunc_Expand(NN, nen, N)
       ! expand [NN] -> [N]
       implicit none
       integer, intent(in) :: nen
       real, intent(in)  :: NN(nen)
       real, intent(out) :: N(NDIM, nen*NDIM)
     end subroutine ShapeFunc_Expand
  end interface

  interface
     subroutine el05_mass(me, XYZ, MatDensity)
       ! Compute element mass matrix for an 8 node hexahedron
       implicit none
       integer, parameter :: nen_e = EIGHT_NODES
       integer, parameter :: nee = NDIM*nen_e
       real, dimension(nee,nee), intent(out) :: me
       real, dimension(nen_e,NDIM), intent(in) :: XYZ
       real, intent(in) :: MatDensity
     end subroutine el05_mass
  end interface

  interface
     subroutine el12_mass(me, XYZ, MatDensity)
       ! Compute element mass matrix for an 8 node hexahedron
       implicit none
       integer, parameter :: nen_e = TWENTYSEVEN_NODES
       integer, parameter :: nee = NDIM*nen_e
       real, dimension(nee,nee), intent(out) :: me
       real, dimension(nen_e,NDIM), intent(in) :: XYZ
       real, intent(in) :: MatDensity
     end subroutine el12_mass
  end interface


  interface

     subroutine el05_stiff_Biot(ke,qs,XYZ,YoungsModulus,PoissonRatio,Qie)
       implicit none
       integer, parameter :: nen_e = EIGHT_NODES
       integer, parameter :: nee = NDIM*nen_e
       real, dimension(nee,nee), intent(out) :: ke
       real, dimension(nee), intent(out)     :: Qs
       real, dimension(nen_e,NDIM), intent(in)   :: XYZ, Qie
       real, intent(in) :: YoungsModulus,PoissonRatio
     end subroutine el05_stiff_Biot

  end interface

  interface

     subroutine el12_stiff_Biot(ke,qs,XYZ,YoungsModulus,PoissonRatio,Qie)
       implicit none
       integer, parameter :: nen_e = TWENTYSEVEN_NODES
       integer, parameter :: nee = NDIM*nen_e
       real, dimension(nee,nee), intent(out) :: ke
       real, dimension(nee), intent(out)     :: Qs
       real, dimension(nen_e,NDIM), intent(in)   :: XYZ, Qie
       real, intent(in) :: YoungsModulus,PoissonRatio
     end subroutine el12_stiff_Biot

  end interface

  interface

     subroutine el05_stiff_kirch(ke,qs,XYZ,YoungsModulus,PoissonRatio,Qie)
       implicit none
       integer, parameter :: nen_e = EIGHT_NODES
       integer, parameter :: nee = NDIM*nen_e
       real, dimension(nee,nee), intent(out) :: ke
       real, dimension(nee), intent(out)     :: Qs
       real, dimension(nen_e,NDIM), intent(in)   :: XYZ, Qie
       real, intent(in) :: YoungsModulus,PoissonRatio
     end subroutine el05_stiff_kirch

  end interface

  interface

     subroutine el12_stiff_kirch(ke,qs,XYZ,YoungsModulus,PoissonRatio,Qie)
       implicit none
       integer, parameter :: nen_e = TWENTYSEVEN_NODES
       integer, parameter :: nee = NDIM*nen_e
       real, dimension(nee,nee), intent(out) :: ke
       real, dimension(nee), intent(out)     :: Qs
       real, dimension(nen_e,NDIM), intent(in)   :: XYZ, Qie
       real, intent(in) :: YoungsModulus,PoissonRatio
     end subroutine el12_stiff_kirch

  end interface

  interface

     subroutine el05_IntForce_kirch(qs,XYZ,YoungsModulus,PoissonRatio,Qie)
       implicit none
       integer, parameter :: nen_e = EIGHT_NODES
       integer, parameter :: nee = NDIM*nen_e
       real, dimension(nee), intent(out)         :: Qs
       real, dimension(nen_e,NDIM), intent(in)   :: XYZ, Qie
       real, intent(in) :: YoungsModulus,PoissonRatio
     end subroutine el05_IntForce_kirch

  end interface
  
  
  interface

     subroutine el12_IntForce_kirch(qs,XYZ,YoungsModulus,PoissonRatio,Qie)
       implicit none
       integer, parameter :: nen_e = TWENTYSEVEN_NODES
       integer, parameter :: nee = NDIM*nen_e
       real, dimension(nee), intent(out)         :: Qs
       real, dimension(nen_e,NDIM), intent(in)   :: XYZ, Qie
       real, intent(in) :: YoungsModulus,PoissonRatio
     end subroutine el12_IntForce_kirch

  end interface

  interface

     subroutine el05_IntForce_biot(qs,XYZ,YoungsModulus,PoissonRatio,Qie)
       implicit none
       integer, parameter :: nen_e = EIGHT_NODES
       integer, parameter :: nee = NDIM*nen_e
       real, dimension(nee), intent(out)         :: Qs
       real, dimension(nen_e,NDIM), intent(in)   :: XYZ, Qie
       real, intent(in) :: YoungsModulus,PoissonRatio
     end subroutine el05_IntForce_biot

  end interface
  
  
  interface

     subroutine el12_IntForce_biot(qs,XYZ,YoungsModulus,PoissonRatio,Qie)
       implicit none
       integer, parameter :: nen_e = TWENTYSEVEN_NODES
       integer, parameter :: nee = NDIM*nen_e
       real, dimension(nee), intent(out)         :: Qs
       real, dimension(nen_e,NDIM), intent(in)   :: XYZ, Qie
       real, intent(in) :: YoungsModulus,PoissonRatio
     end subroutine el12_IntForce_biot

  end interface
  
  interface 
     subroutine getTri_Norm(Positions,Elements,normals,centers,areas, volume,&
          NNodes,NEle)
       
#include "Flash.h"
       
       implicit none
       interface
          function cross_product(a,b)   
            real ::a(3),b(3),cross_product(3)        
          end function cross_product
       end interface
       integer,intent(IN) ::  Nele, Nnodes
       integer,intent(IN) :: Elements(3,Nele)
       real,intent(IN)    :: Positions(Nnodes*3)
       real,intent(OUT)   :: normals(Nele*NDIM),centers(Nele*NDIM)
       real,intent(OUT)   :: areas(Nele),volume
     end subroutine getTri_Norm
  end interface
  
  interface
     subroutine sm_3DFlexible_getElement_EvalMax(lambda_max, body, e)
       use SolidMechanics_data, only: sm_structure
       implicit none
       real,               intent(out) :: lambda_max
       type(sm_structure), intent(in)  :: body
       integer,            intent(in)  :: e
     end subroutine sm_3DFlexible_getElement_EvalMax
  end interface


  ! Rigid Body Functions
  interface
     subroutine el15_mass(maxdofs,ibd,M)
       implicit none
       integer, intent(in) :: maxdofs,ibd
       real, intent(out) :: M(maxdofs,maxdofs)
     end subroutine el15_mass
  end interface


  ! Misc:
  interface
     subroutine el10_lengths(XYZ,Qie,lengths)
       implicit none
       integer, parameter :: nen_e = NINE_NODES
       real, dimension(nen_e,NDIM), intent(in)  :: XYZ, Qie
       real, dimension(2),          intent(out) :: lengths
     end subroutine el10_lengths
  end interface

  interface
     subroutine sm_el01_countParticles(body, e, Dmin, nXi, nEta, ptelem, flag)
       use SolidMechanics_data, only: sm_structure
       implicit none
       ! IO Variables
       type(sm_structure)   :: body     ! entire body structure
       integer, intent(in)  :: e        ! element number
       real, intent(in)     :: Dmin     ! Eulerian Distance
       integer, intent(out) :: nXi, nEta, ptelem ! number of points in Xi, Eta, and total
       integer, intent(in)  :: flag
     end subroutine sm_el01_countParticles
  end interface

  interface
     subroutine sm_el02_countParticles(body, e, Dmin, nXi, nEta, ptelem, flag)
       use SolidMechanics_data, only: sm_structure
       implicit none
       type(sm_structure)   :: body     ! entire body structure
       integer, intent(in)  :: e        ! element number
       real, intent(in)     :: Dmin     ! Eulerian Distance
       integer, intent(out) :: nXi, nEta, ptelem ! number of points in Xi, Eta, and total
       integer, intent(in)  :: flag
     end subroutine sm_el02_countParticles
  end interface

  interface
     subroutine sm_el10_countParticles(body, e, Dmin, nXi, nEta, ptelem, flag)
       use SolidMechanics_data, only: sm_structure
       implicit none
       type(sm_structure)   :: body     ! entire body structure
       integer, intent(in)  :: e        ! element number
       real, intent(in)     :: Dmin     ! Eulerian Distance
       integer, intent(out) :: nXi, nEta, ptelem ! number of points in Xi, Eta, and total
       integer, intent(in)  :: flag
     end subroutine sm_el10_countParticles
  end interface

  interface
     subroutine sm_el01_mapParticles(body, e, ptelem,  &
                                     xpos,ypos,        &
                                     xvel,yvel,        &
                                     xacc,yacc,        &
                                     xnrm,ynrm,        &
                                     areai, loc_num )
       use SolidMechanics_data, only: sm_structure
       implicit none

       ! IO Variables
       type(sm_structure)   :: body     ! entire body structure
       integer, intent(in)  :: e        ! element number
       integer, intent(in)  :: ptelem
       real, dimension(ptelem) :: xpos, ypos, xvel, yvel, xacc, yacc, xnrm, ynrm, areai
       integer, dimension(ptelem) :: loc_num
     end subroutine sm_el01_mapParticles
  end interface

  interface 
     subroutine sm_el02_mapParticles(body, e, ptelem, xpos,ypos,zpos, &
                                     xvel,yvel,zvel,xacc,yacc,zacc,   &
                                     xnrm,ynrm,znrm,areai,loc_num )
       use SolidMechanics_data, only: sm_structure
       implicit none
       type(sm_structure)   :: body     ! entire body structure
       integer, intent(in)  :: e        ! element number
       integer, intent(in)  :: ptelem
       real, dimension(ptelem) :: xpos, ypos, zpos, xvel, yvel, zvel, xacc, yacc, zacc, xnrm, ynrm, znrm, areai
       integer, dimension(ptelem) :: loc_num       
     end subroutine sm_el02_mapParticles
  end interface

  interface 
     subroutine sm_el10_mapParticles(body, e, ptelem, xpos,ypos,zpos, &
                                     xvel,yvel,zvel,xacc,yacc,zacc,   &
                                     xnrm,ynrm,znrm,areai,loc_num )
       use SolidMechanics_data, only: sm_structure
       implicit none
       type(sm_structure)   :: body     ! entire body structure
       integer, intent(in)  :: e        ! element number
       integer, intent(in)  :: ptelem
       real, dimension(ptelem) :: xpos, ypos, zpos, xvel, yvel, zvel, xacc, yacc, zacc, xnrm, ynrm, znrm, areai
       integer, dimension(ptelem) :: loc_num       
     end subroutine sm_el10_mapParticles
  end interface

  interface
     subroutine sm_el10_FluidForce( body, e, particle, Hsp_pres, Hsp_visc )
       use SolidMechanics_data, only: sm_structure
       implicit none

       ! constants
       integer, parameter :: nen_e = NINE_NODES
       integer, parameter :: nee   = NDIM*nen_e

       ! IO Variables
       type(sm_structure),  intent(in) :: body
       integer, intent(in) :: e  ! element number on the wet surface
       real, dimension(NPART_PROPS), intent(in)  :: particle
       real, dimension(nee), intent(out) :: Hsp_pres, Hsp_visc
     end subroutine sm_el10_FluidForce
  end interface

  interface
     subroutine sm_el10_ga_stiff(body, e, particle, beta, dt, kf )
       use SolidMechanics_data, only: sm_structure
       implicit none
       integer, parameter :: nen_e = NINE_NODES
       integer, parameter :: nee   = NDIM*nen_e
       type(sm_structure),           intent(in)  :: body
       integer,                      intent(in)  :: e  ! element number on the wet surface
       real, dimension(NPART_PROPS), intent(in)  :: particle
       real,                         intent(in)  :: beta, dt
       real, dimension(nee,nee),     intent(out) :: kf
     end subroutine sm_el10_ga_stiff
  end interface

  interface 
     subroutine sm_el02_FluidForce(body, e, particle, Hsp_pres, Hsp_visc )
       use SolidMechanics_data, only: sm_structure
       
       implicit none
       
       ! constants
       integer, parameter :: nen_e = THREE_NODES
       integer, parameter :: nee   = NDIM*nen_e
       
       ! IO Variables
       type(sm_structure),  intent(in) :: body
       integer, intent(in) :: e  ! element number on the wet surface
       real, dimension(NPART_PROPS), intent(in)  :: particle
       real, dimension(nee), intent(out) :: Hsp_pres, Hsp_visc

     end subroutine sm_el02_FluidForce
  end interface
  
  interface
     subroutine sm_el12_COM(COM, mass, XYZ, Qie,MatDensity)
       implicit none
       integer, parameter :: nen_e = 27
       integer, parameter :: nee = NDIM*nen_e
       real, dimension(NDIM), intent(out) :: COM
       real                 , intent(out) :: mass
       real, dimension(nen_e,NDIM), intent(in) :: XYZ, Qie
       real, intent(in) :: MatDensity
     end subroutine sm_el12_COM
  end interface
 
  interface
     subroutine el03_IntForce_rbc(p1,p2,p3,  &
          area,Aoj,At, Aot, &
          Vt, Vot, normal,center, tforces)

#include "Flash.h"
#include "SolidMechanics.h"
       implicit none

       Interface
          function cross_product(a,b)
            real ::a(3),b(3),cross_product(3)
          end function cross_product
       end interface

       real,dimension(NDIM), intent(IN) :: p1,p2,p3
       real,dimension(NDIM), intent(IN) :: normal,center
       real,dimension(MAXNODERBC,NDIM),intent(OUT)::tforces
       real,intent(IN) ::Aoj,Vt,area,At
       real,intent(IN) ::Aot,Vot
     end subroutine el03_IntForce_rbc
  end interface

  interface
     subroutine el02_IntForce_rbc(p1,p2,     &
     v1,v2,     &
     ks1,ks2,Lmax_M,modelType, &
     sforces)

#include "Flash.h"
#include "SolidMechanics.h"

       implicit none
       real,intent(IN)  :: p1(NDIM),p2(NDIM)
       real,dimension(MAXNODERBC,NDIM),intent(OUT) :: sforces
       real,intent(IN)  :: ks1,ks2,Lmax_M,v1(NDIM),v2(NDIM)
       integer, intent(IN) :: modelType
     end subroutine el02_IntForce_rbc
  end interface

  interface
     subroutine el04_IntForce_rbc (p1,p2,p3,p4,  &
     n1,n2,c1,c2,a1,a2,th,        &
     tho,cos_tho,sin_tho,    &
     bforces)
#include "Flash.h"

       implicit none

       interface
          function cross_product(a,b)
            real ::a(3),b(3),cross_product(3)
          end function cross_product
       end interface

       real,dimension(NDIM),intent(IN)::p1,p2,p3,p4
       real,dimension(NDIM),intent(IN)::n1,n2,c1,c2
       real, intent(IN)  :: a1, a2
       real, intent(IN)  :: tho,cos_tho,sin_tho
       real, intent(OUT) :: bforces(4,NDIM),th
     end subroutine el04_IntForce_rbc
  end interface

end module sm_element_interface

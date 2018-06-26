
#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

module SolidMechanics_data

implicit none

integer, save :: sm_MeshMe, sm_NumProcs, sm_meshComm

integer, parameter, dimension(19)    :: sm_nen = (/ TWO_NODES,     &       ! eltype=1
                                                    THREE_NODES,   &       ! eltype=2, TRIANGLE
                                                    FOUR_NODES,    &       ! eltype=3
                                                    FOUR_NODES,    &       ! eltype=4
                                                    EIGHT_NODES,   &       ! eltype=5
                                                    SIX_NODES,     &       ! eltype=6, PRISM
                                                    FIVE_NODES,    &       ! eltype=7
                                                    THREE_NODES,   &       ! eltype=8
                                                    SIX_NODES,     &       ! eltype=9
                                                    NINE_NODES,    &       ! eltype=10
                                                    TEN_NODES,     &       ! eltype=11
                                                    TWENTYSEVEN_NODES, &   ! eltype=12
                                                    EIGHTEEN_NODES,&       ! eltype=13 
                                                    FOURTEEN_NODES,&       ! eltype=14
                                                    ONE_NODE,      &       ! eltype=15
                                                    EIGHT_NODES,   &       ! eltype=16
                                                    TWENTY_NODES,  &       ! eltype=17
                                                    FIFTEEN_NODES, &       ! eltype=18
                                                    THIRTEEN_NODES /)      ! eltype=19


! Gravity components:
real, save :: sm_gravX,sm_gravY,sm_gravZ 

! Number of Bodies:
integer, save :: sm_NumBodies

! Group hangle for global comm in SolidMechanics:
integer, save :: sm_GlobalGroup

! Restraint Array for the structures, this is by degree of freedom.
type sm_restr
   integer :: restnode        ! Node restrained in Global numbering.
   integer :: restdof         ! Dof Being restrained in local node numbering.
   integer :: restype,nparam 
   real, allocatable, dimension(:) :: param !(NMAXPARAMRES)
end Type sm_restr

! Restraint array for restrained surfaces.
type sm_restr_surf
      integer :: kinematics_idx
      integer :: nfix
      integer :: restdof   ! DOF direction being restrained (0=All, 1=IAXIS, 2=JAXIS, 3=KAXIS)
      integer, allocatable, dimension(:) :: node_list ! Global node numbers
end type sm_restr_surf

! Transformation matrices data structure
Type trmat
   integer :: fixed_angvar(MDIM)
   real :: TPB(MDIM,MDIM),TNB(MDIM,MDIM),TNBo(MDIM,MDIM),NBB(MDIM,MDIM),NBB_dot(MDIM,MDIM)
   real :: PwB(MDIM,1),NwB_DF(MDIM,1),NwB_N(MDIM,1)
   real :: PaB(MDIM,1),NaB_DF(MDIM,1),NaB_N(MDIM,1)
end Type trmat


! Structural Data per body:
type sm_structure
   
   ! Body Master, we assume for now that the integration of structural equations of each 
   ! body is done serially by the processor BodyMaster. This means global matrix/vector
   ! assemblies, and integration step are done by the master.
   integer :: BodyMaster

   ! Body Type
   integer :: BodyType

   ! Meta body, case of ibd containing only part of of bodys wet surface.
   integer :: MetaBody

   ! Communicator of of MetaBody this body belongs to:
   integer :: mbcomm

   ! Rank of processor on MetaBody Group:
   integer :: mbcomm_Group
   integer :: mbcomm_Me 
   integer :: mbcomm_NumProcs

   ! Integrator Type (def's in SolidMechanics.h)
   integer :: IntegMethod
   integer :: IC_flag

   ! Number of nodes
   integer :: nnp
   ! Number of elmenets
   integer :: nel,max_eltype
   ! Number of links and adjacent triangles
   integer ::nlinks,Nele
   ! Number of general dofs, dofs, and restrained dofs
   integer :: ndofs, neq, ng, max_dofs_per_node
   
   ! XYZ
   real, allocatable, dimension(:) :: x,y,z

   ! ElEM IEN
   integer, allocatable, dimension(:,:) :: IEN
   integer, allocatable, dimension(:,:) :: Tri,TRILinks
   integer, allocatable, dimension(:,:) :: Links
   integer, allocatable, dimension(:,:) :: bend_pts
   ! ELEM Properties
   ! ID
   integer, allocatable, dimension(:,:) :: ID, LM
   ! Constraints

   ! We only do Isotropic.    
   real, allocatable, dimension(:) :: MatDensity
   real, allocatable, dimension(:) :: YoungsModulus
   real, allocatable, dimension(:) :: PoissonsRatio
   real, allocatable, dimension(:) :: Normals,Centers

   ! Material Type:   
   integer, allocatable, dimension(:) :: MatType
   ! Elemenet Type (GMSH numbered)
   integer, allocatable, dimension(:) :: eltype,new_e
   
   !
   ! Structural degrees of freedom
   !
   ! {q} = { qf(1:neq), v(1:ng) }
   real, allocatable, dimension(:)  :: qi
   
   
   ! RBC Material data, ks, Kp, 
   real :: p, p_M, tho, cos_tho, sin_tho, Dist_avg
   real :: Aot,Vot, r_mult, C1, Vol, At, Vt
   real ,dimension(:),allocatable :: kp, ks, Lo
   real ,dimension(:),allocatable :: Dist,kp_M, ks_M, Lo_M, Lmax_M
   real ,dimension(:),allocatable :: Aoj,Areas,comm
   
   !
   ! Structural Matrices in Compressed Sparse Format
   !
   ! Constrained Boundary Condition Information for csc K
   ! [K] = [ K_qq, T_qv 
   !         T_vq, K_vv  ]
   ! {q} = { qf(1:neq), v(1:ng) }
   ! 
   ! Free-Free Part of Matrix
   integer :: qq_nnz
   real, allocatable, dimension(:) :: K, M, Damp
   integer, allocatable, dimension(:) :: qq_IA, qq_JA

   ! Free-Constrained Part of Matrix
   integer :: qv_nnz, qv_nfix
   real, allocatable, dimension(:) :: Kqv, Mqv, Dampqv
   integer, allocatable, dimension(:) :: qv_IA, qv_JA, fix_list_A

   ! Fancy compressed builder mapping
   integer, allocatable, dimension(:,:,:) :: LM_cs
   integer, allocatable, dimension(:,:,:) :: ws_LM_cs
   
   ! Internal forces
   real, allocatable, dimension(:) :: Qs
   
   ! RBC internal forces 
   real, allocatable, dimension(:,:):: Qs_rbc
   ! Proportional Damping Info
   !
   integer :: damping_flag
   real :: damping_kappa_m, damping_kappa_k0

   ! Gravity parameters
   integer :: gravity_flag
   real, dimension(MDIM) :: gravity_vec 
   
   !
   ! Wetted surface
   !
   integer :: ws_nel,ws_max_eltype
   integer, allocatable, dimension(:,:) :: ws_IEN
   integer, allocatable, dimension(:)   :: ws_eltype
   integer, allocatable, dimension(:)   :: ws_nXi, ws_nEta, ws_ptelem
   integer, allocatable, dimension(:)   :: ws_nxL,ws_nyL,ws_nzL       ! Normals outside at wet surface points.

   !
   ! Variables for Dynamic Time marching
   !
   ! Internal and external forces and n and n-1
   real, allocatable, dimension(:) :: Qsn, Hs, Hsn
   ! Special Containers to keep the fluid forces on each dof
   real, allocatable, dimension(:) :: Hs_pres, Hs_visc, Hsi_pres, Hsi_visc

   ! work array for the right hand side
   real, allocatable, dimension(:) :: dyn_rhs

   ! position, velocity, and accel of q at time n 
   real, allocatable, dimension(:) :: qn,qdn,qddn

   ! Variables for Predictor-Corrector integration, only allocated if 
   ! PC integrator is selected:
   ! qn, qdn, qddn used for predicted solution
   ! qi, qdi, qddi used for corrected solution
   real, allocatable, dimension(:) :: qdi,qddi ! qi is already defined above

   ! Variables for steps n-1,n-2,n-3,n-nmult
   real, allocatable, dimension(:,:) :: qms, qdms, qddms
   real, allocatable, dimension(:,:) :: ey, edy
   
   ! Contact force from step n-1 
   real, allocatable, dimension(:,:) :: Cddms

   ! Restraints derived Type
   type(sm_restr),  allocatable, dimension(:) :: restraints ! Allocate with number of Restrained structural Nodes in Body
   integer :: nrcoords,maxrestparams

   ! Restraints on surfaces
   type(sm_restr_surf), allocatable, dimension(:) :: restraints_surf
   integer :: nrsurf

   ! SuperLU Stuff
   ! Variable used by SuperLU for mass matrix LU decomposition
   integer*8 :: lu_factors

   ! Constant Mass Matrix Flag
   integer :: flag_constMass
 
   ! Rigid Body Properties:
   real, allocatable, dimension(:) :: xB,yB,zB
   integer :: ix,ex,ia,ea,iw,ew
   integer :: trmatrix
   type(trmat) :: RBMAT
   integer :: borigin_node                         ! The origin of the Body Frame node is the Center of Mass
   real, allocatable, dimension(:) :: stiff
   real :: mass                                    ! Body Mass.
   real :: Volume                                  ! Body Volume.
   real, allocatable, dimension(:,:) :: I_body     ! Inertia in the reference body axes.
   real, allocatable, dimension(:,:) :: I_newton   ! Inertia Moment in Newtonian axes.
   real, allocatable, dimension(:,:) :: M_rigid    ! Unrestrained Mass matrix
   real, allocatable, dimension(:,:) :: Mqv_rigid  ! Restrained Mass Matrix

   ! Analytical body description:
   integer :: flag_forceinside
   integer :: annbody_type,annbody_nparam
   real, allocatable, dimension(:) :: annbody_param 

   ! Containers for fluid integrals on body volume:
   real :: int_rhog(MDIM), int_xrhog(MDIM)
   real :: int_rhou(MDIM,-1:CONSTANT_ONE)
   real :: int_xrhou(MDIM,-1:CONSTANT_ONE)

   ! Center of Mass Location
   real, dimension(MDIM) :: COM

   ! Is Body crossing Boundaries? Periodic BCs
   logical, dimension(LOW:HIGH,MDIM) :: OnBoundary  
  
   ! Use minimum or maximum delta xi to define distance
   ! among markers in mapping
   integer :: map_use_dmin = SM_TRUE

end type sm_structure

! One to one correspondence with gr_sbBodyInfo for aero-grids.
type(sm_structure), save, dimension(:), pointer :: sm_BodyInfo

! Number of MetaBodies:
integer, save :: sm_metaBodies

! Metabody data array:
type sm_metabod

   ! Array to contain how many Processors (or body entries in sm_bodyInfo) per metabody:
   ! The sum of entries in the sm_metaBody_Procs array must be equal to sm_NumBodies,
   ! in case of using metaBodies (sm_metaBodies > 0).
   integer :: NumProcs

   ! Number of Bodies on this metaBody:
   integer :: NumBods 

   ! Array that contains Processor Group handle for each metaBody:
   integer :: Group 

   ! Communicator for metabody, value also given in sm_bodyInfo%mbcomm for each
   ! metabody component ibd:
   integer :: comm

   ! Processors that handle MetaBody, global sm_meshMe numbering:
   integer, allocatable, dimension(:) :: Members

end type sm_metabod

! One to one correspondence with sm_metaBodInfo and sm_BodyInfo%MetaBody.
type(sm_metabod), save, dimension(:), pointer :: sm_metaBodInfo



end module SolidMechanics_data

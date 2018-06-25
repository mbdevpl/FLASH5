!     


!

subroutine sm_ioRead_rigid(ibd,file)

#include "SolidMechanics.h"
#include "constants.h"
#include "Flash.h"
#include "sm_pk.h"
    
  use SolidMechanics_Data, only : sm_MeshMe, sm_structure, sm_BodyInfo, sm_nen, &
                                  sm_gravX, sm_gravY, sm_gravZ
  use sm_iointerface, only: sm_io_checkHdfErr

  use Driver_interface, only : Driver_abortFlash

  USE HDF5
    
  implicit none
    
  integer, intent(IN) :: ibd
  integer(HID_T), intent(IN) :: file

  INTEGER :: h5err
  INTEGER(HID_T) :: dset
  INTEGER(HSIZE_T), DIMENSION(3) :: dims ! size read/write buffer
  integer :: irest
  integer, allocatable, dimension(:) :: fix_list_A
  integer :: nee
  logical :: link_exists

  type(sm_structure), pointer :: body
    
  integer :: i,A,idof,inod,nfix,ircoord,Tot_coord
  integer, allocatable, dimension(:) :: auxvect_restnod
  real, allocatable, dimension(:,:)  :: auxvect_params

  real :: xi,yi,zi

  integer :: metaProcs, metaMype, ws_totnel
  integer, allocatable, dimension(:,:) :: ws_totIEN
  integer, allocatable, dimension(:)   :: ws_toteltype
  integer :: numel,iel,ielloc,istart,iend

  real, parameter :: eps = 1.e-12 

  body => sm_BodyInfo(ibd)

  !
  ! Read in scalars
  !
  dims = (/1,1,1/)
  ! nnp:
  CALL h5dopen_f (file, "mesh/nnp", dset, h5err) ! dset handle to datase required
  call sm_io_checkHdfErr(h5err, 'mesh/nnp') ! error check for nnp
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%nnp, dims, h5err) !use handle, dump integer in nnp
  CALL h5dclose_f(dset , h5err)
  ! nel:
  CALL h5dopen_f (file, "body/nel", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/nel')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%nel, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! max_eltype: Element Type with max number of nodes.
  CALL h5dopen_f (file, "body/max_eltype", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/max_eltype')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%max_eltype, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! dof_per_node:
  CALL h5dopen_f (file, "body/DOFS_per_node", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/DOFS_per_node')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%max_dofs_per_node, dims, h5err) ! All nodes have these dofs
  CALL h5dclose_f(dset , h5err)
!  ndmax = body%max_dofs_per_node
    
  !
  ! Read in nodal points, Material Reference.
  !
  dims(1) = Body%nnp
  allocate(Body%xB(Body%nnp),Body%x(Body%nnp))
  allocate(Body%yB(Body%nnp),Body%y(Body%nnp))
#if NDIM == MDIM
  allocate(Body%zB(Body%nnp),Body%z(Body%nnp))
#endif
  
  ! x:
  CALL h5dopen_f (file, "mesh/x", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'mesh/x')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%xB, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! y:
  CALL h5dopen_f (file, "mesh/y", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'mesh/y')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%yB, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! z:
#if NDIM == MDIM
  CALL h5dopen_f (file, "mesh/z", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'mesh/z')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%zB, dims, h5err)
  CALL h5dclose_f(dset , h5err)
#endif  

  !
  ! Read in Connectivity Arrays
  !
  ! IEN: a.k.a. ELEM
  dims = (/sm_nen(Body%max_eltype),Body%nel,1/)
  allocate(Body%IEN(sm_nen(Body%max_eltype),Body%nel))
  CALL h5dopen_f (file, "body/IEN", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/IEN')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%IEN, dims, h5err)
  CALL h5dclose_f(dset , h5err)

  ! eltype:
  dims = (/Body%nel,1,1/)
  allocate(Body%eltype(Body%nel))
  CALL h5dopen_f (file, "body/eltype", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/eltype')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%eltype, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  
  !
  ! Material Information
  !
#if NDIM == MDIM
  allocate(Body%Stiff(NDIM)) 
  Body%Stiff = 0.0
  dims = (/NDIM,1,1/)
  ! Rigid Body Stiffness in x,y,z:
  CALL h5dopen_f (file, "body/Stiffness", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/Stiffness')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%Stiff, dims, h5err)
  CALL h5dclose_f(dset , h5err) 

#else
  allocate(Body%Stiff(NDIM+1)) ! [kx ky ktheta]
  Body%Stiff = 0.0
  dims = (/NDIM+1,1,1/)
  ! Rigid Body Stiffness in x,y,z:
  CALL h5dopen_f (file, "body/Stiffness", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/Stiffness')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%Stiff, dims, h5err)
  CALL h5dclose_f(dset , h5err)
#endif

  ! Rigid Body Mass:
  dims = (/1,1,1/)
  CALL h5dopen_f (file, "body/Mass", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/Mass')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%Mass, dims, h5err)
  CALL h5dclose_f(dset , h5err) 
  body%flag_constMass = SM_FALSE ! Always build mass matrix.

  ! Rigid Body Volume:
  call h5lexists_f( file, "body/Volume", link_exists, h5err )
  if( link_exists ) then
    CALL h5dopen_f (file, "body/Volume", dset, h5err)
    call sm_io_checkHdfErr(h5err, 'body/Volume')
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%Volume, dims, h5err)
    CALL h5dclose_f(dset , h5err)
  else
    Body%Volume = 0.0
  endif

  ! Inertia in Body Reference:
  allocate(Body%I_body(NDIM,NDIM),Body%I_newton(NDIM,NDIM)) 
  Body%I_body = 0.0
  dims = (/NDIM,NDIM,1/)
  ! Rigid Body Stiffness in x,y,z:
  CALL h5dopen_f (file, "body/I_body", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/I_body')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%I_body, dims, h5err)
  CALL h5dclose_f(dset , h5err) 

  ! Gravity flag:
  dims = (/1,1,1/)
  CALL h5dopen_f (file, "body/gravity_flag", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/gravity_flag')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%gravity_flag, dims, h5err)
  CALL h5dclose_f(dset , h5err)  

  ! Abort if body has no gravity force but flash.par has nonzero grav components:
  if ( (Body%gravity_flag .eq. SM_FALSE) .and. & 
       (abs(sm_gravX*sm_gravY*sm_gravZ) .gt. eps) ) then
     write(*,*) 'Body=',ibd,', Body%gravity_flag=',Body%gravity_flag
     write(*,*) 'Gravity=',sm_gravX,sm_gravY,sm_gravZ   
     call Driver_abortFlash("Body gravity_flag set to false, when flash.par gravity components non-zero.")
  endif

  ! Gravity vector:
  dims = (/NDIM,1,1/)
  call h5lexists_f( file, "body/gravity", link_exists, h5err )
  if( link_exists ) then
    ! Read body gravity:
    CALL h5dopen_f (file, "body/gravity", dset, h5err)
    call sm_io_checkHdfErr(h5err, 'body/gravity')
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%gravity_vec, dims, h5err)
    CALL h5dclose_f(dset , h5err)   

    ! Check body gravity and flash.par gravity vectors are the same:
    if (abs(Body%gravity_vec(IAXIS)-sm_gravX) .gt. eps) &
      call Driver_abortFlash("Body X direction gravity component not equal to flash.par gravX.")

    if (abs(Body%gravity_vec(JAXIS)-sm_gravY) .gt. eps) &
      call Driver_abortFlash("Body Y direction gravity component not equal to flash.par gravY.")
#if NDIM == MDIM
    if (abs(Body%gravity_vec(KAXIS)-sm_gravZ) .gt. eps) &
      call Driver_abortFlash("Body Z direction gravity component not equal to flash.par gravZ.")
#endif

  elseif (Body%gravity_flag .eq. SM_TRUE) then
    ! load flash.par gravity values to Body%gravity_vec
    ! sm_gravX,sm_gravY,sm_gravZ already loaded in Solidmechanics_init.F90
    Body%gravity_vec(IAXIS) = sm_gravX
    Body%gravity_vec(JAXIS) = sm_gravY
    Body%gravity_vec(KAXIS) = sm_gravZ
  endif  

  ! Rotation Matrix type:
  dims = (/1,1,1/)
  CALL h5dopen_f (file, "body/trmatrix", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/trmatrix')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%trmatrix, dims, h5err)
  CALL h5dclose_f(dset , h5err)   

  ! Assign values of beginning and en of displacement, orientation and anv vel variables:
#if NDIM == MDIM
  Body % ix = 1; Body % ex = 3;
#else
  Body % ix = 1; Body % ex = 2;
#endif

  select case( Body%trmatrix )

  case( RB_IDENTITY )

     if(Body%max_dofs_per_node .ne. 6) then
        write(*,*) 'Body=',ibd,', max dofs per node=',Body%max_dofs_per_node 
        call Driver_abortFlash("sm_ioRead_rigid: for RB_IDENTITY max dofs per node .ne. 6")
     endif

     Body % ia = -1; Body % ea = -1;
     Body % iw =  4; Body % ew =  6;

  case( RB_EULER321 )

     if(Body%max_dofs_per_node .ne. 9) then
        write(*,*) 'Body=',ibd,', max dofs per node=',Body%max_dofs_per_node
        call Driver_abortFlash("sm_ioRead_rigid: for RB_EULER321 max dofs per node .ne. 9")
     endif

     Body % ia = 4; Body % ea = 6;
     Body % iw = 7; Body % ew = 9;

  case( RB_QUATERNN )

     if(Body%max_dofs_per_node .ne. 10) then
        write(*,*) 'Body=',ibd,', max dofs per node=',Body%max_dofs_per_node 
        call Driver_abortFlash("sm_ioRead_rigid: for RB_EULER321 max dofs per node .ne. 10")
     endif

     Body % ia = 4; Body % ea = 7;
     Body % iw = 8; Body % ew =10;

  case( RB_TWODIM )

     if(Body%max_dofs_per_node .ne. 4) then
        write(*,*) 'Body=',ibd,', max dofs per node=',Body%max_dofs_per_node
        call Driver_abortFlash("sm_ioRead_rigid: for RB_TWODIM max dofs per node .ne. 4")
     endif

     Body % ia = 3; Body % ea = 3;
     Body % iw = 4; Body % ew = 4;


  case default
     
     call Driver_abortFlash("sm_ioRead_rigid: error, body transformation matrix not known.")

  end select


  ! Body system origin point, always first point:
  Body % Borigin_node = BODYFRAME_NODE


  ! Read if body has analytical description an internal IB forcing is required:
  dims = (/ 1, 1, 1 /)

  ! flag for forcing inside IBody:
  Body % flag_forceinside = SM_FALSE
  Body % annbody_type     = -99999
  Body % annbody_nparam   = 0
  call h5lexists_f( file, "body/flag_forceinside", link_exists, h5err )  

  if( link_exists ) then

     ! Value of flag:
     CALL h5dopen_f (file, "body/flag_forceinside", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'body/flag_forceinside')
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%flag_forceinside, dims, h5err)
     CALL h5dclose_f(dset , h5err)

     if (Body%flag_forceinside .eq. SM_TRUE) then

        ! Get annbody_type and annbody_nparam:
        CALL h5dopen_f (file, "body/annbody_type", dset, h5err)
        call sm_io_checkHdfErr(h5err, 'body/annbody_type')
        CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%annbody_type, dims, h5err)
        CALL h5dclose_f(dset , h5err)

        CALL h5dopen_f (file, "body/annbody_nparam", dset, h5err)
        call sm_io_checkHdfErr(h5err, 'body/annbody_nparam')
        CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%annbody_nparam, dims, h5err)
        CALL h5dclose_f(dset , h5err)

        ! Allocate vector of shape parameters for analytical body:
        allocate ( body % annbody_param  ( Body%annbody_nparam ) )

        ! Read analytical body parameters:
        dims = (/Body%annbody_nparam,1,1/)        
        CALL h5dopen_f (file, "body/annbody_param", dset, h5err)
        call sm_io_checkHdfErr(h5err, 'body/annbody_param')
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%annbody_param, dims, h5err)
        CALL h5dclose_f(dset , h5err)   

     else

        call h5lexists_f( file, "body/annbody_type", link_exists, h5err )

        if( link_exists ) then

        ! Get annbody_type and annbody_nparam:
        CALL h5dopen_f (file, "body/annbody_type", dset, h5err)
        call sm_io_checkHdfErr(h5err, 'body/annbody_type')
        CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%annbody_type, dims, h5err)
        CALL h5dclose_f(dset , h5err)

        CALL h5dopen_f (file, "body/annbody_nparam", dset, h5err)
        call sm_io_checkHdfErr(h5err, 'body/annbody_nparam')
        CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%annbody_nparam, dims, h5err)
        CALL h5dclose_f(dset , h5err)

        ! Allocate vector of shape parameters for analytical body:
        allocate ( body % annbody_param  ( Body%annbody_nparam ) )

        ! Read analytical body parameters:
        dims = (/Body%annbody_nparam,1,1/)
        CALL h5dopen_f (file, "body/annbody_param", dset, h5err)
        call sm_io_checkHdfErr(h5err, 'body/annbody_param')
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%annbody_param, dims, h5err)
        CALL h5dclose_f(dset , h5err)
        
        endif

     endif

  endif

  ! Check if max Eulerian dx should be used for markers mapping
  dims = (/ 1, 1, 1/)
  call h5lexists_f( file, "body/map_use_dmin", link_exists, h5err )

  if( link_exists ) then
     ! Value of flag:
     CALL h5dopen_f (file, "body/map_use_dmin", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'body/map_use_dmin')
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%map_use_dmin, dims, h5err)
     CALL h5dclose_f(dset , h5err)
  endif

  ! Restrained Nodes:
  dims = (/ 1, 1, 1 /)
  ! ng: Restrained dofs
  call h5lexists_f( file, "/RestSurface/nrestsurf", link_exists, h5err )

  if( link_exists ) then


     ! How Many Restrained coordinates:
     CALL h5dopen_f (file, "RestSurface/nrestsurf", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestSurface/nrestsurf')
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%nrsurf, dims, h5err) ! Restrained Dofs
     CALL h5dclose_f(dset , h5err)

     if (Body%nrsurf .eq. 1) then
        irest = Body%nrsurf
     else
        call Driver_abortFlash("sm_ioRead_rigid: one restrained surface maximum for now.")
     endif

     ! fix_list_A: list of global node #'s on constrained surface
     dims = (/1, 1, 1 /)
     CALL h5dopen_f (file, "RestSurface/nfix", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestSurface/nfix')
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%qv_nfix, dims, h5err)
     CALL h5dclose_f(dset , h5err)
     dims = (/Body%qv_nfix, 1, 1 /)
     allocate(fix_list_A(Body%qv_nfix))
     CALL h5dopen_f (file, "RestSurface/fix_list", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestSurface/fix_list')
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, fix_list_A, dims, h5err)
     CALL h5dclose_f(dset , h5err)
  
     ! populate the restraint_surf obj
     allocate( body%restraints_surf(irest) ) !!! Assuming that there is only 1 surface per body being 
                                             !!! restrained.  can always change later to be more general
     dims = (/1, 1, 1 /)
     CALL h5dopen_f (file, "RestSurface/kinematics_idx", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestSurface/kinematics_idx') ! Type of Restraint
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, body%restraints_surf(irest)%kinematics_idx, dims, h5err)
     CALL h5dclose_f(dset , h5err)
     allocate( body%restraints_surf(irest)%node_list( body%qv_nfix ) )
     body%restraints_surf(irest)%node_list = fix_list_A
     body%restraints_surf(irest)%nfix = body%qv_nfix

     deallocate(fix_list_A)

  else

     Body % nrsurf = 0

  endif


  ! Restrained Nodes:
  dims = (/ 1, 1, 1 /)
  ! ng: Restrained dofs
  call h5lexists_f( file, "/RestNodes/nrestcoords", link_exists, h5err )


  if( link_exists ) then

     ! How Many Restrained coordinates:
     CALL h5dopen_f (file, "RestNodes/nrestcoords", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestNodes/nrestcoords')
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%nrcoords, dims, h5err) ! Restrained Dofs
     CALL h5dclose_f(dset , h5err)

     if ( Body%nrcoords .gt. 0 ) then 

     allocate ( body%restraints  ( Body%nrcoords ) )
     allocate ( auxvect_restnod  ( Body%nrcoords ) )

     ! Maximum restraint Parameters
     CALL h5dopen_f (file, "RestNodes/maxrestparams", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestNodes/maxrestparams')
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%maxrestparams, dims, h5err) ! Restrained Dofs
     CALL h5dclose_f(dset , h5err)

     allocate ( auxvect_params(Body%maxrestparams,Body%nrcoords) )

     ! Node number in global coordinates for restraint irnode:
     dims = (/ Body%nrcoords, 1, 1 /)
     call h5dopen_f (file, "RestNodes/nodes", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestNodes/nodes')          
     call h5dread_f(dset, H5T_NATIVE_INTEGER, auxvect_restnod, dims, h5err)
     call h5dclose_f(dset , h5err)

     do ircoord = 1,Body%nrcoords
        body%restraints(ircoord)%restnode = auxvect_restnod(ircoord)
     enddo

     ! Integer with restrained dof for the node:
     dims = (/ Body%nrcoords,1,1 /)
     call h5dopen_f (file, "RestNodes/restdof", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestNodes/restdof')          
     call h5dread_f(dset, H5T_NATIVE_INTEGER, auxvect_restnod, dims, h5err)
     call h5dclose_f(dset , h5err)     

     do ircoord = 1,Body%nrcoords
        body%restraints(ircoord)%restdof = auxvect_restnod(ircoord)
     enddo

     ! Integer with type of restraint for dof of the node:
     dims = (/ Body%nrcoords,1,1 /)
     call h5dopen_f (file, "RestNodes/restype", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestNodes/restype')          
     call h5dread_f(dset, H5T_NATIVE_INTEGER,  auxvect_restnod, dims, h5err)
     call h5dclose_f(dset , h5err)     

     do ircoord = 1,Body%nrcoords
        body%restraints(ircoord)%restype = auxvect_restnod(ircoord)
     enddo
     

     ! Integer with number of parameters per dof restraint of the node:
     dims = (/ Body%nrcoords,1,1 /)
     call h5dopen_f (file, "RestNodes/nparam", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestNodes/nparam')          
     call h5dread_f(dset, H5T_NATIVE_INTEGER, auxvect_restnod, dims, h5err)
     call h5dclose_f(dset , h5err)  

     do ircoord = 1,Body%nrcoords
        body%restraints(ircoord)%nparam = auxvect_restnod(ircoord)
     enddo


     ! Read the restraint parameters for the node:
     dims = (/ Body%maxrestparams,Body%nrcoords,1 /)
     call h5dopen_f (file, "RestNodes/param", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestNodes/param')          
     call h5dread_f(dset, H5T_NATIVE_DOUBLE, auxvect_params, dims, h5err)
     call h5dclose_f(dset , h5err) 

     do ircoord = 1,Body%nrcoords
        allocate (  body%restraints(ircoord)%param(Body%maxrestparams) )
        body%restraints(ircoord)%param(1:Body%maxrestparams)=auxvect_params(1:Body%maxrestparams,ircoord)
     enddo
  
     deallocate(auxvect_restnod)
     deallocate(auxvect_params)   

     endif
 
  else

     Body%nrcoords = 0

  end if

  ! Which orientation vars are constrained?
!!$  body % RBMAT % fixed_angvar(:) = 0
!!$  do ircoord=1,Body%nrcoords   
!!$
!!$     if ( body%restraints(ircoord)%restnode .eq. Body % Borigin_node) then
!!$
!!$        ! Set fixed angular var flags: 
!!$        ! Case of using Euler angles for Orientation:
!!$        select case(body%restraints(ircoord)%restdof)
!!$        case(MDIM+1) ! Phi
!!$           body % RBMAT % fixed_angvar(ANGLE_1) = 1
!!$        case(MDIM+2) ! Theta
!!$           body % RBMAT % fixed_angvar(ANGLE_2) = 1
!!$        case(MDIM+3) ! Psi
!!$           body % RBMAT % fixed_angvar(ANGLE_3) = 1
!!$        end select
!!$
!!$        ! Change corresponding angular velocity var constraint to
!!$        ! SM_PK_CANGVEL
!!$        ! Case of using Euler angles for Orientation:
!!$        if (body%restraints(ircoord)%restdof .eq. body%iw  ) &
!!$           body%restraints(ircoord)%restype = SM_PK_CANGVEL
!!$        if (body%restraints(ircoord)%restdof .eq. body%iw+1) &
!!$           body%restraints(ircoord)%restype = SM_PK_CANGVEL
!!$        if (body%restraints(ircoord)%restdof .eq. body%iw+2) &
!!$           body%restraints(ircoord)%restype = SM_PK_CANGVEL
!!$
!!$     endif
!!$
!!$  enddo


  !
  ! Load Wet Surface info
  !
  ! ws_max_eltype
  dims = (/1,1,1/)
  CALL h5dopen_f (file, "WetSurface/max_eltype", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'WetSurface/max_eltype')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%ws_max_eltype, dims, h5err)
  CALL h5dclose_f(dset , h5err)      
  

  ! Split wet surface among metabody components:
  if (Body%MetaBody .gt. CONSTANT_ZERO) then

    metaProcs = Body%mbcomm_NumProcs
    metaMype  = Body%mbcomm_Me

    ! Read Total number of surface elements:
    ! ws_totnel
    dims = (/1,1,1/)
    CALL h5dopen_f (file, "WetSurface/nel", dset, h5err)
    call sm_io_checkHdfErr(h5err, 'WetSurface/nel')
    CALL h5dread_f(dset, H5T_NATIVE_INTEGER, ws_totnel, dims, h5err)
    CALL h5dclose_f(dset , h5err)          

    ! ws_totIEN
    dims = (/sm_nen(Body%ws_max_eltype),ws_totnel,1/)
    allocate(ws_totIEN(sm_nen(Body%ws_max_eltype),ws_totnel))
    CALL h5dopen_f (file, "WetSurface/IEN", dset, h5err)
    call sm_io_checkHdfErr(h5err, 'WetSurface/IEN')
    CALL h5dread_f(dset, H5T_NATIVE_INTEGER, ws_totIEN, dims, h5err)
    CALL h5dclose_f(dset , h5err)       

    ! ws_toteltype
    dims = (/ws_totnel,1,1/)
    allocate(ws_toteltype(ws_totnel))
    CALL h5dopen_f (file, "WetSurface/eltype", dset, h5err)
    call sm_io_checkHdfErr(h5err, 'WetSurface/eltype')
    CALL h5dread_f(dset, H5T_NATIVE_INTEGER, ws_toteltype, dims, h5err)
    CALL h5dclose_f(dset , h5err)
     
    ! Now figure out start and end indexes to be read into Body%ws_IEN:
    numel = ws_totnel / metaProcs

    if (numel .le. CONSTANT_ZERO) &
      call Driver_abortFlash("sm_ioRead_rigid : ws_totnel less than metaProcs.")
    
    ! Last processor takes care of the remainder elements:
    istart = metaMype*numel + 1
    if ((metaMype+1) .eq. metaProcs) then ! Last Processor
      iend = ws_totnel
    else
      iend = (metaMype+1)*numel
    endif

    ! Final numel:
    numel = iend - istart + 1

    ! Allocate and Assign:
    Body%ws_nel = numel
    allocate(Body%ws_IEN(sm_nen(Body%ws_max_eltype),Body%ws_nel))
    allocate(Body%ws_eltype(Body%ws_nel))
    ielloc = 0
    do iel = istart , iend 
      ielloc = ielloc + 1
      Body%ws_IEN(:,ielloc)  = ws_totIEN(:,iel)
      Body%ws_eltype(ielloc) = ws_toteltype(iel)
    enddo

    !write(*,*) Body%MetaBody,Body%BodyMaster,metaMype,metaProcs,Body%ws_nel,ws_totnel

    deallocate(ws_totIEN,ws_toteltype)

  else

    ! ws_nel
    dims = (/1,1,1/)
    CALL h5dopen_f (file, "WetSurface/nel", dset, h5err)
    call sm_io_checkHdfErr(h5err, 'WetSurface/nel')
    CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%ws_nel, dims, h5err)
    CALL h5dclose_f(dset , h5err)    
    ! ws_IEN:
    dims = (/sm_nen(Body%ws_max_eltype),Body%ws_nel,1/)
    allocate(Body%ws_IEN(sm_nen(Body%ws_max_eltype),Body%ws_nel))
    CALL h5dopen_f (file, "WetSurface/IEN", dset, h5err)
    call sm_io_checkHdfErr(h5err, 'WetSurface/IEN')
    CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%ws_IEN, dims, h5err)
    CALL h5dclose_f(dset , h5err)        
    ! ws_eltype:
    dims = (/Body%ws_nel,1,1/)
    allocate(Body%ws_eltype(Body%ws_nel))
    CALL h5dopen_f (file, "WetSurface/eltype", dset, h5err)
    call sm_io_checkHdfErr(h5err, 'WetSurface/eltype')
    CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%ws_eltype, dims, h5err)
    CALL h5dclose_f(dset , h5err)
  
  endif



  ! allocate other ws_* arrays
  allocate( body%ws_nXi (   body%ws_nel ), &
            body%ws_nEta(   body%ws_nel ), &
            body%ws_ptelem( body%ws_nel )   )
  

  !
  ! Read in info about Integration Method of choice
  !
  call h5lexists_f( file, "/dynamic/IntegMethod", link_exists, h5err )
  if( link_exists ) then
      dims = (/1,1,1/)
      call h5dopen_f(file, "/dynamic/IntegMethod", dset, h5err )
      call h5dread_f(dset, H5T_NATIVE_INTEGER, body%IntegMethod, dims, h5err)
      call h5dclose_f(dset, h5err)
  else
      body%IntegMethod = SOLIDINTEG_PREDCORR
  end if

  ! 
  ! See if IC's are to be loaded
  !
  call h5lexists_f( file, "/dynamic/IC_flag", link_exists, h5err )
  if( link_exists ) then
      dims = (/1,1,1/)
      call h5dopen_f(file, "/dynamic/IC_flag", dset, h5err )
      call h5dread_f(dset, H5T_NATIVE_INTEGER, body%IC_flag, dims, h5err)
      call h5dclose_f(dset, h5err)
  else
      body%IC_flag = SM_IC_ZERO
  end if

   
  !
  ! DOF storage
  !


  ! Total Number of Gen Coordinates:
  Tot_coord = body%max_dofs_per_node * body%nnp

  ! Allocate ID:
  allocate( body%ID(body%max_dofs_per_node,body%nnp) )

  ! Set ID01
  body%ID = 1

  ! Set to zero restraints
  if(Body%nrcoords .gt. 0) then ! We have Dof Restraints

    do ircoord = 1,Body%nrcoords
       idof = body%restraints(ircoord)%restdof
       inod = body%restraints(ircoord)%restnode 
       body%ID(idof,inod) = 0;
    enddo

  endif

  if(Body%nrsurf .gt. 0) then ! We have Surface Restraints

     do irest=1,Body%nrsurf     
        nfix = body % qv_nfix
        do i=1,nfix
          inod = body%restraints_surf(irest)%node_list(i)
          body%ID(:,inod) = 0;
        enddo
     enddo

  endif

  ! Number Dofs-Restraints:
  Body%neq = 0
  Body%ng  = 0
  do A=1,Body%nnp
     do i=1,Body%max_dofs_per_node

        if(body%ID(i,A) .eq. 1) then ! Dof
           Body%neq = Body%neq + 1 
           body%ID(i,A) = Body%neq
        elseif(body%ID(i,A) .eq. 0) then
           body%ID(i,A) = Tot_coord - Body%ng
           Body%ng = Body%ng + 1
        endif

     enddo
  enddo

  !write(*,*) 'NEQ   =',Body%neq
  !write(*,*) 'NRESTR=',Body%ng
  !write(*,*) 'TOTAL =',Tot_coord

  if (Tot_coord .ne. (Body%ng+Body%neq)) &
      call Driver_abortFlash("sm_ioRead_rigid : Eqs + restraints not adding to NCOORDS")

!!$  do A=1,Body%nnp
!!$     write(*,*) A,body%ID(:,A)
!!$  enddo

  Body%ndofs = Body%ng + Body%neq
  allocate( Body%Qs(Body%neq) ) ! This can be a zero size allocated array - OO - MV
  Body%Qs = 0.                  

    
  allocate( Body%qi(Body%ndofs), Body%qn(Body%ndofs) )
  Body%qi = 0.
  Body%qn = 0.
  allocate( Body%qdn(Body%ndofs), Body%qddn(Body%ndofs) )
  Body%qdn  = 0.
  Body%qddn = 0.

  ! System Matrices
  allocate(Body%M_rigid(Body%neq,Body%neq),Body%Mqv_rigid(Body%neq,Body%ng))

!!$    allocate( Body%qdn1(Body%ndofs), Body%qddn1(Body%ndofs) )
!!$    Body%qdn1  = 0.
!!$    Body%qddn1 = 0.
!!$    
!!$  ! Initialize the large matrices
!!$  allocate( Body%K(Body%qq_nnz), Body%M(Body%qq_nnz), Body%Kqv(Body%qv_nnz), Body%Mqv(Body%qv_nnz) )
!!$  Body%K = 0.
!!$  Body%M = 0.
!!$  Body%Kqv = 0.
!!$  Body%Mqv = 0.
!!$  
!!$  ! Init Damping stuff
!!$  if( Body%damping_flag == 1) then
!!$     allocate(Body%Damp(Body%qq_nnz), Body%Dampqv(Body%qv_nnz))
!!$     Body%Damp   = 0.
!!$     Body%Dampqv = 0.
!!$  end if
!!$  

  ! Initialize the Dyn storage:
  allocate( Body%Hs(Body%neq),Body%Qsn(Body%neq),Body%Hsn(Body%neq),Body%dyn_rhs(Body%neq) )
  Body%Qsn     = 0.
  Body%Hs      = 0.
  Body%Hsn     = 0.
  Body%dyn_rhs = 0.    

  ! Initialize viscous and pressure Hs
  allocate( Body%Hs_visc(Body%ndofs),Body%Hs_pres(Body%ndofs) )
  Body%Hs_visc = 0.
  Body%Hs_pres = 0.

  ! This body belongs to me:
  Body%BodyMaster = sm_MeshMe

  !write(*,*) 'Nodes=',Body%nnp
  !write(*,*) 'Gen Coords=',Body%ndofs
  !write(*,*) 'Eqn=',Body%neq,', rest=',Body%ng
  

 
!!$  write(*,'(A)') 'complete.'
  
end subroutine sm_ioRead_rigid






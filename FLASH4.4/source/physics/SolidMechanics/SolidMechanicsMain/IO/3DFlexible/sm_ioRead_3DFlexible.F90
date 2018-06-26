!!****if* source/physics/SolidMechanics/SolidMechanicsMain/IO/3DFlexible/sm_ioRead_3DFlexible
!!
!!
!! NAME
!!
!!
!!
!!
!! SYNOPSIS
!!  Load in a 3DFlexible body
!!
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_ioRead_3DFlexible(ibd,file)

  use SolidMechanics_Data, only : sm_MeshMe,sm_structure,sm_BodyInfo, sm_nen
  use sm_iointerface, only: sm_io_checkHdfErr, sm_ioRead_3DFlexible_loadIC

  USE HDF5

  implicit none

  ! IO Variables
  integer, intent(IN) :: ibd
  integer(HID_T), intent(in) :: file

  ! Internal variables for reading
  INTEGER :: h5err, group_id
  INTEGER(HID_T) :: dset
  INTEGER(HSIZE_T), DIMENSION(3) :: dims ! size read/write buffer
  character(len=100) ::   group_name

  ! Internal working variables
  type(sm_structure), pointer :: body
  integer :: irest
  integer :: nee, ndmax, nfix, ws_nee
  logical :: link_exists

  ! Meta-body variables
  integer :: metaProcs, metaMype, ws_totnel
  integer, allocatable, dimension(:,:) :: ws_totIEN
  integer, allocatable, dimension(:)   :: ws_toteltype
  integer :: numel, iel, ielloc, istart, iend
  character(len=100) :: metabod_str
  integer :: metabod_unit = 998

  ! Local pointer
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
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%max_dofs_per_node, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ndmax = body%max_dofs_per_node

  !
  ! Read in nodal points, Material Reference.
  !
  dims(1) = Body%nnp
  allocate(Body%x(Body%nnp))
  allocate(Body%y(Body%nnp))
#if NDIM == MDIM
  allocate(Body%z(Body%nnp))
#endif

  ! x:
  CALL h5dopen_f (file, "mesh/x", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'mesh/x')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%x, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! y:
  CALL h5dopen_f (file, "mesh/y", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'mesh/y')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%y, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! z:
#if NDIM == MDIM
  CALL h5dopen_f (file, "mesh/z", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'mesh/z')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%z, dims, h5err)
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

  ! ID
  dims = (/ndmax,Body%nnp,1/)
  allocate(Body%ID(ndmax,Body%nnp))
  CALL h5dopen_f (file, "body/ID", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/ID')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%ID, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! LM
  nee  = ndmax*sm_nen(Body%max_eltype) ! MAXDOFS per element
  dims = (/nee,Body%nel,1/)
  allocate(Body%LM(nee,Body%nel))
  CALL h5dopen_f (file, "body/LM", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/LM')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%LM, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! eltype: Type of element as in Gmsh2
  dims = (/Body%nel,1,1/)
  allocate(Body%eltype(Body%nel))
  CALL h5dopen_f (file, "body/eltype", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/eltype')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%eltype, dims, h5err)
  CALL h5dclose_f(dset , h5err)

  !
  ! Material Information
  !
  dims(1) = Body%nel
  ! Material Type: One per element.
  allocate(Body%MatType(Body%nel))
  CALL h5dopen_f (file, "body/MatType", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/MatType')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%MatType, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! Young's Modulus of Elasticity:
  allocate(Body%YoungsModulus(Body%nel))
  CALL h5dopen_f (file, "body/YoungsModulus", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/YoungsModulus')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%YoungsModulus, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! Poisson's Ratio:
  allocate(Body%PoissonsRatio(Body%nel))
  CALL h5dopen_f (file, "body/PoissonRatio", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/PoissonRatio')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%PoissonsRatio, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! Material Density:
  allocate(Body%MatDensity(Body%nel))
  CALL h5dopen_f (file, "body/MatDensity", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/MatDensity')
  CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%MatDensity, dims, h5err)
  CALL h5dclose_f(dset , h5err)

  !
  ! Proportional Damping
  !
  call h5lexists_f( file, "body/ProportionalDamping/flag", link_exists, h5err )
  if( link_exists ) then
     ! Damping Flag:
     dims = (/1,1,1/)
     CALL h5dopen_f (file, "body/ProportionalDamping/flag", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'body/ProportionalDamping/flag')
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%damping_flag, dims, h5err)
     CALL h5dclose_f(dset , h5err)
     ! Damping Factors:
     if( body%damping_flag == SM_TRUE ) then
        CALL h5dopen_f (file, "body/ProportionalDamping/kappa_M", dset, h5err)
        call sm_io_checkHdfErr(h5err, 'body/ProportionalDamping/kappa_M')
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%damping_kappa_M, dims, h5err)
        CALL h5dclose_f(dset , h5err)
        CALL h5dopen_f (file, "body/ProportionalDamping/kappa_K0", dset, h5err)
        call sm_io_checkHdfErr(h5err, 'body/ProportionalDamping/kappa_K0')
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%damping_kappa_K0, dims, h5err)
        CALL h5dclose_f(dset , h5err)
     end if
  else
     write(*,'(A)',advance='no') ' *** No damping information found'
     Body%damping_flag = SM_FALSE
  end if

  !
  ! Matrix Information
  !
  ! LM_csc(q1,q2,e) -> K(b) in csc format:
  dims = (/ nee, nee, Body%nel/)
  allocate( Body%LM_cs(nee, nee, Body%nel) )
  CALL h5dopen_f (file, "body/csc/LM_csc", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/csc/LM_csc')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%LM_cs, dims, h5err)
  CALL h5dclose_f(dset , h5err)

  !
  ! Matrix partitioning masks for free-free, rest.-free parts
  !
  dims = (/ 1, 1, 1 /)
  ! ng: Restrained dofs
  CALL h5dopen_f (file, "RestSurface/ng", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'RestSurface/ng')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%ng, dims, h5err) ! Restrained Dofs
  CALL h5dclose_f(dset , h5err)
  ! nq: Dofs
  CALL h5dopen_f (file, "RestSurface/neq", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'RestSurface/neq')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%neq, dims, h5err) ! Dofs
  CALL h5dclose_f(dset , h5err)
  ! nnz_qq:
  CALL h5dopen_f (file, "body/csc/qq_nnz", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/csc/qq_nnz')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%qq_nnz, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! nnz_qv:
  CALL h5dopen_f (file, "RestSurface/csc/qv_nnz", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'RestSurface/csc/qv_nnz')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%qv_nnz, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! Masks for the free-free part
  ! qq_JA:
  dims = (/Body%neq+1, 1, 1 /)
  allocate(Body%qq_JA(Body%neq+1))
  CALL h5dopen_f (file, "body/csc/qq_JA", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/csc/qq_JA')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%qq_JA, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! qq_IA
  dims = (/Body%qq_nnz, 1, 1 /)
  allocate(Body%qq_IA(Body%qq_nnz))
  CALL h5dopen_f (file, "body/csc/qq_IA", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/csc/qq_IA')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%qq_IA, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! Masks for the free-constrained part
  ! qv_JA:
  dims = (/Body%ng+1, 1, 1 /)
  allocate(Body%qv_JA(Body%ng+1))
  CALL h5dopen_f (file, "RestSurface/csc/qv_JA", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'RestSurface/csc/qv_JA')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%qv_JA, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  ! qv_IA
  dims = (/Body%qv_nnz, 1, 1 /)
  allocate(Body%qv_IA(Body%qv_nnz))
  CALL h5dopen_f (file, "RestSurface/csc/qv_IA", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'RestSurface/csc/qv_IA')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%qv_IA, dims, h5err)
  CALL h5dclose_f(dset , h5err)

  ! Restrained Nodes:
  Body%nrcoords = 0

  ! Restrained Surfaces:
  dims = (/ 1, 1, 1 /)
  call h5lexists_f( file, "/RestSurface/nrestsurf", link_exists, h5err )
  if( link_exists ) then

     CALL h5dopen_f (file, "RestSurface/nrestsurf", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'RestSurface/nrestsurf')
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%nrsurf, dims, h5err) ! Restrained Dofs
     CALL h5dclose_f(dset , h5err)

     ! populate the restraint_surf obj
     allocate( body%restraints_surf(body%nrsurf) )
     Body%qv_nfix = 0
     do irest = 1,body%nrsurf

        ! open the group
        group_name(:) = ''
        write(group_name,"(A,I1,A)") '/RestSurface/',irest,'/'
        call h5gopen_f( file, trim(group_name), group_id, h5err )

        ! fix_list_A: list of global node #'s on constrained surface
        dims = (/1, 1, 1 /)
        CALL h5dopen_f (group_id, "nfix", dset, h5err)
        call sm_io_checkHdfErr(h5err, 'RestSurface/nfix')
        CALL h5dread_f(dset, H5T_NATIVE_INTEGER, nfix, dims, h5err)
        CALL h5dclose_f(dset , h5err)

        Body%restraints_surf(irest)%nfix = nfix
        Body%qv_nfix = Body%qv_nfix + nfix

        allocate( body%restraints_surf(irest)%node_list( nfix ) )

        ! Read in fix_list
        dims = (/nfix,1,1/)
        CALL h5dopen_f (group_id, "fix_list", dset, h5err)
        call sm_io_checkHdfErr(h5err, 'RestSurface/fix_list')
        CALL h5dread_f(dset, H5T_NATIVE_INTEGER, body%restraints_surf(irest)%node_list, dims, h5err)
        CALL h5dclose_f(dset , h5err)

        ! Read in the data
        dims = (/1,1,1/)
        CALL h5dopen_f (group_id, "kinematics_idx", dset, h5err)
        call sm_io_checkHdfErr(h5err, 'RestSurface/kinematics_idx')
        CALL h5dread_f(dset, H5T_NATIVE_INTEGER, body%restraints_surf(irest)%kinematics_idx, dims, h5err)
        CALL h5dclose_f(dset , h5err)


        ! Flag to specify if only certain dof of the nodes on the surface are to be constrained
        call h5lexists_f( group_id, "restdof", link_exists, h5err )
        if( link_exists ) then
           
           ! if the flag is present in the bodyfile then load the dof (=1,2,3)
           dims = (/ 1, 1, 1 /)
           call h5dopen_f (group_id, "restdof", dset, h5err)
           call sm_io_checkHdfErr(h5err, 'RestSurface/restdof')
           call h5dread_f(dset, H5T_NATIVE_INTEGER, body%restraints_surf(irest)%restdof, dims, h5err)
           call h5dclose_f(dset , h5err)

        else
           ! the flag is not present so default to all the dof on the node
           body%restraints_surf(irest)%restdof = ALL_DOF
        end if

        ! close the group
        call h5gclose_f(group_id, h5err)

     end do

  else
     body%nrsurf = 0
     write(*,*) 'No surface restraints.'
  end if

  !
  ! Load Wet Surface info
  !
  ! ws_max_eltype
  dims = (/1,1,1/)
  CALL h5dopen_f (file, "WetSurface/max_eltype", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'WetSurface/max_eltype')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%ws_max_eltype, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  if (Body%MetaBody .gt. CONSTANT_ZERO) then
     !*** Split wet surface among metabody components:

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

     ! Write the start and stop numbers of the ws_IEN as seen for each metabody
     write(metabod_str,"(A,I5.5,A)") 'IOData/sm_particles_wsIEN.',ibd,'.dat'
     open(unit=metabod_unit, file=trim(metabod_str), action="write", status='replace')
     write(metabod_unit, "(A, I)") 'ibd = ', ibd
     write(metabod_unit, "(A)") 'Global element start and stop numbers'
     write(metabod_unit, "(A,I)") "istart = ", istart
     write(metabod_unit, "(A,I)") "iend = ", iend
     close(metabod_unit)

     deallocate(ws_totIEN,ws_toteltype)

  else
     !*** Use just a single meta-body (old way)
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
     ! ws_LM_cs(q1,q2,e) -> K(b) in csc format:
     ws_nee = sm_nen(Body%ws_max_eltype)*NDIM
     dims = (/ ws_nee, ws_nee, Body%ws_nel/)
     allocate( Body%ws_LM_cs(ws_nee, ws_nee, Body%ws_nel) )
     CALL h5dopen_f (file, "WetSurface/LM_load_csc", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'WetSurface/ws_LM_csc')
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%ws_LM_cs, dims, h5err)
     CALL h5dclose_f(dset , h5err)

  end if

  ! allocate other ws_* arrays
  allocate( body%ws_nXi (   body%ws_nel ), &
       body%ws_nEta(   body%ws_nel ), &
       body%ws_ptelem( body%ws_nel )   )

  !
  ! Read in info about constant mass matrix
  !
  call h5lexists_f( file, "/dynamic/flag_ConstMass", link_exists, h5err )
  if( link_exists ) then
     dims = (/1,1,1/)
     call h5dopen_f(file, "/dynamic/flag_ConstMass", dset, h5err )
     call h5dread_f(dset, H5T_NATIVE_INTEGER, body%flag_ConstMass, dims, h5err)
     call h5dclose_f(dset, h5err)
  else
     body%flag_ConstMass = SM_FALSE
  end if

  !
  ! Gravity Flag:
  !
  call h5lexists_f( file, "/body/gravity_flag", link_exists, h5err )
  if( link_exists ) then
     dims = (/1,1,1/)
     CALL h5dopen_f (file, "body/gravity_flag", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'body/gravity_flag')
     CALL h5dread_f(dset, H5T_NATIVE_INTEGER, body%gravity_flag, dims, h5err)
     CALL h5dclose_f(dset , h5err)
     ! Gravity Vector:
     dims = (/NDIM,1,1/)
     CALL h5dopen_f (file, "body/gravity", dset, h5err)
     call sm_io_checkHdfErr(h5err, 'body/gravity')
     CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%gravity_vec, dims, h5err)
     CALL h5dclose_f(dset , h5err)
  else
     body%gravity_flag = SM_FALSE
     body%gravity_vec(:) = 0.
  end if

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
  !write(*,'(A)',advance='no') 'Allocating Main Arrays...'
  Body%ndofs = Body%ng + Body%neq

  allocate( Body%Qs(Body%neq) )
  Body%Qs = 0.

  allocate( Body%qi(Body%ndofs), Body%qn(Body%ndofs) )
  Body%qi = 0.
  Body%qn = 0.
  allocate( Body%qdn(Body%ndofs), Body%qddn(Body%ndofs) )
  Body%qdn  = 0.
  Body%qddn = 0.
!!$    allocate( Body%qdn1(Body%ndofs), Body%qddn1(Body%ndofs) )
!!$    Body%qdn1  = 0.
!!$    Body%qddn1 = 0.

  ! Initialize the large matrices
  allocate( Body%K(Body%qq_nnz), Body%M(Body%qq_nnz), Body%Kqv(Body%qv_nnz), Body%Mqv(Body%qv_nnz) )
  Body%K = 0.
  Body%M = 0.
  Body%Kqv = 0.
  Body%Mqv = 0.

  ! Init Damping stuff
  if( Body%damping_flag == SM_TRUE ) then
     allocate(Body%Damp(Body%qq_nnz), Body%Dampqv(Body%qv_nnz))
     Body%Damp   = 0.
     Body%Dampqv = 0.
  end if

  ! Initialize the Dyn storage:
  allocate( Body%dyn_rhs(Body%neq) )
  Body%dyn_rhs = 0.
  allocate( Body%Hs(Body%neq) )
  Body%Hs      = 0.
  !allocate( Body%Hsn(Body%neq) )
  !Body%Hsn     = 0.
  !allocate( Body%Qsn(Body%neq) )
  !Body%Qsn     = 0.

  ! Fluid Storage of every dof (not just the unconstrained ones)
  allocate( Body%Hs_pres(Body%ndofs), Body%Hs_visc(Body%ndofs) )
  Body%Hs_pres = 0.
  Body%Hs_visc = 0.
  allocate( Body%Hsi_pres(Body%ndofs), Body%Hsi_visc(Body%ndofs) )
  Body%Hsi_pres = 0.
  Body%Hsi_visc = 0.

  ! Load in the IC's from external file: sm_body.IC.$ibd$.h5
  if( body%IC_flag == SM_IC_EXTFILE ) then
     call sm_ioRead_3DFlexible_loadIC(ibd)
  end if

  ! This body belongs to me:
  Body%BodyMaster = sm_MeshMe

end subroutine sm_ioRead_3DFlexible



!!$    !
!!$    ! Read in Generalized Alpha Stuff:
!!$    !
!!$    dims = (/ 1, 1, 1 /)
!!$    ! dyn_rhoinf:
!!$    CALL h5dopen_f (file, "dynamic/rho_inf", dset, h5err)
!!$    call check_hdferr(h5err, 'dynamic/rho_inf')
!!$    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, dyn_rhoinf, dims, h5err)
!!$    CALL h5dclose_f(dset , h5err)
!!$    ! delta_t
!!$    CALL h5dopen_f (file, "dynamic/timestep", dset, h5err)
!!$    call check_hdferr(h5err, 'dynamic/timestep')
!!$    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, delta_t, dims, h5err)
!!$    CALL h5dclose_f(dset , h5err)
!!$    ! dyn_epsilon
!!$    CALL h5dopen_f (file, "dynamic/epsilon", dset, h5err)
!!$    call check_hdferr(h5err, 'dynamic/epsilon')
!!$    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, dyn_epsilon, dims, h5err)
!!$    CALL h5dclose_f(dset , h5err)
!!$    ! dyn_ntime
!!$    CALL h5dopen_f (file, "dynamic/ntime", dset, h5err)
!!$    call check_hdferr(h5err, 'dynamic/ntime')
!!$    CALL h5dread_f(dset, H5T_NATIVE_INTEGER, dyn_ntime, dims, h5err)
!!$    CALL h5dclose_f(dset , h5err)
!!$    ! dyn_maxiter
!!$    CALL h5dopen_f (file, "dynamic/max_sub_iter", dset, h5err)
!!$    call check_hdferr(h5err, 'dynamic/max_sub_iter')
!!$    CALL h5dread_f(dset, H5T_NATIVE_INTEGER, dyn_maxiter, dims, h5err)
!!$    CALL h5dclose_f(dset , h5err)
!!$    ! dyn_IC_flag

!!$    ! Initialize the Dyn storage:
!!$    allocate( qd_alpha_f(ndofs) , qdd_alpha_m(ndofs), Qsn(neq), Hs(neq), Hsn(neq), dyn_rhs(neq) )
!!$    qd_alpha_f = 0.
!!$    qdd_alpha_m = 0.
!!$    Qsn = 0.
!!$    Hs  = 0.
!!$    Hsn = 0.
!!$    dyn_rhs = 0.

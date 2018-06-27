!     
! File:   init.F95
! Author: tim
! Edit by Hussein Ezzeldin hmezz@gwu.edu Jan 2013 

subroutine sm_ioRead_rbc(ibd,file)

#include "SolidMechanics.h"
#include "Flash.h"
#include "constants.h"

  use SolidMechanics_Data, only : sm_MeshMe,sm_structure,sm_BodyInfo,sm_nen, sm_Numbodies
  use sm_Verlet_data, only: sm_Verlet_type, sm_Verlet_info
  use sm_iointerface, only: sm_io_checkHdfErr,sm_io_rbcinit
  use Simulation_data, Only : Npmax, Npmin, Ypmax, Nplus, &
                              Nminus,per
  use SolidMechanics_rbc_data,Only: Fext,Fext_M, N_M,stretching_exp
  USE HDF5
    
  implicit none
  
  ! Argument list
  integer, intent(IN) :: ibd
  integer(HID_T),intent(in) :: file

  !local variables
  character(LEN=25) :: filename ,file_ibd; 
  integer :: h5err,i
  integer(HID_T) :: dset
  integer(HSIZE_T), DIMENSION(3) :: dims ! size read/write buffer
  integer :: nee,irest, sp_modelType
  integer, allocatable, dimension(:) :: fix_list_A
  integer :: Nele, nLinks, Extpts
  type(sm_structure), pointer :: body
  integer,save :: nnp
  type(sm_Verlet_type), pointer :: integ
  integer, allocatable :: indx(:), indy(:)
  integer, allocatable :: inds(:)
  call sm_io_rbcinit()
  
  body => sm_BodyInfo(ibd)
  ! Check if sm_Verlet_info has been "allocated"
  if( .not. associated( sm_Verlet_info ) ) then
      allocate( sm_Verlet_info( sm_NumBodies ) )
  end if
  
  ! set integ pointer:
  integ => sm_Verlet_info(ibd)
  
  !
  ! Read in scalars
  !
  dims = (/1,1,1/)
  ! nnp:
  CALL h5dopen_f (file, "mesh/nnp", dset, h5err) ! dset handle to datase required
  call sm_io_checkHdfErr(h5err, 'mesh/nnp') ! error check for nnp
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%nnp, dims, h5err) !use handle, dump integer in nnp
  CALL h5dclose_f(dset , h5err)
  nnp=Body%nnp;
  ! nel:
  CALL h5dopen_f (file, "body/nel", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/nel')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%nel, dims, h5err)
  CALL h5dclose_f(dset , h5err)

  ! integ%dt
   CALL h5dopen_f (file, "body/dt", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/dt')
  CALL h5dread_f(dset,H5T_NATIVE_DOUBLE, integ%dt, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  write(*,*) 'body: ',ibd,'integ%dt: ',integ%dt
  ! max_eltype:
  CALL h5dopen_f (file, "body/nLinks", dset, h5err)
  CALL sm_io_checkHdfErr(h5err, 'body/nLinks')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%nLinks, dims, h5err)
  CALL h5dclose_f(dset , h5err)
    
  ! Nele
  CALL h5dopen_f (file, "body/nele", dset, h5err)
  CALL sm_io_checkHdfErr(h5err, 'body/nele')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%Nele, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  !
  ! Setting the wet surface variables 
  !

  
  allocate(inds(body%nnp ))
  !write(*,*) 'nnp',nnp
  inds(:)=(/(i,i=1,body%nnp*3,3)/)
  body%neq = NDIM* Body%nnp
  !
  ! Read in nodal points
  !
  dims(1) = Body%nnp
  allocate(Body%x(Body%nnp))
  allocate(Body%y(Body%nnp))
  allocate(indx(body%nnp),indy(body%nnp))
  ! Allocate Positions, velocities and forces at time n-1 (qms) , n (qn) , and intermediate (qi)
  allocate (Body%qdms(Body%nnp*NDIM,CONSTANT_ONE),& 
            Body%qddms(Body%nnp*NDIM,CONSTANT_ONE),&
            Body%Cddms(Body%nnp*NDIM,CONSTANT_ONE))
  allocate (Body%qn(Body%nnp*NDIM), Body%qdn(Body%nnp*NDIM), Body%qddn(Body%nnp*NDIM))
  allocate (Body%qdi(Body%nnp*NDIM))
  body%max_dofs_per_node= NDIM;
  Body%qdms(:,1)  = 0.;
  if (ibd==1) then
     Body%qdms(inds,1) = 2.;
  elseif (ibd==2) then
     Body%qdms(inds,1) = -2.;
     end if
     Body%qddms(:,1) = 0.0 
     Body%qn(:)      = 0.;
     Body%qdn(:)     = 0.;
     Body%qddn(:)    = 0.;
     Body%qdi(:)     = 0.;
     
#if NDIM == MDIM
  allocate(Body%z(Body%nnp))  
#endif

  allocate (Body%qms(Body%nnp*NDIM,CONSTANT_ONE))
  ! x:
  call h5dopen_f (file, "mesh/x", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'mesh/x')
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%x, dims, h5err)
  call h5dclose_f(dset , h5err)


  !y:
  call h5dopen_f (file, "mesh/y", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'mesh/y')
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%y, dims, h5err)
  call h5dclose_f(dset , h5err)

  do i=1,body%nnp
     body%qms(NDIM*i-2,1)=body%x(i);
     body%qms(NDIM*i-1,1)=body%y(i);
     body%qn(NDIM*i-2)=body%x(i);
     body%qn(NDIM*i-1)=body%y(i);
  end do
  
  write(*,*) 'The bounds of body ',ibd
  write(*,*) 'Xmin= ',minval(body%x),' Xmax= ',maxval(body%x)
  write(*,*) 'Ymin= ',minval(body%y),' Ymax= ',maxval(body%y)

  ! z:
#if NDIM == MDIM
  call h5dopen_f (file, "mesh/z", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'mesh/z')
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, Body%z, dims, h5err)
  call h5dclose_f(dset , h5err)
  write(*,*) 'Zmin= ',minval(body%z),' Zmax= ',maxval(body%z(:))
  do i=1,body%nnp     
     body%qms(NDIM*i,1)=body%z(i);
     body%qn(NDIM*i)=body%z(i);
  end do
  
#endif  


 
 
  !
  ! Read in Connectivity Arrays
  !
  ! IEN: a.k.a. ELEM or Tri
  
  allocate(Body%Tri(3,Body%nele))
  allocate(Body%Links(2,Body%nlinks))
  allocate(Body%TRILinks(2,Body%nlinks))
  allocate(Body%bend_pts(4,Body%nlinks))
  allocate(Body%IEN(4,Body%nel))
  allocate(Body%ID(NDIM,Body%nnp));
  allocate(Body%MatType(Body%nel));
  allocate(Body%elType(Body%nel));
  allocate(Body%new_e(Body%nel));

!!$  allocate(Body%ws_IEN())
!!$  allocate(ws_eltype())
!!$  allocate(ws_nXi(), ws_nEta(), ws_ptelem())
  
  dims = (/2,Body%nLinks,1/)
  call h5dopen_f (file, "body/Links", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/Links')
  call h5dread_f(dset, H5T_NATIVE_INTEGER, Body%Links, dims, h5err)
  call h5dclose_f(dset , h5err) 

  call h5dopen_f (file, "body/TRILinks", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/TRILinks')
  call h5dread_f(dset, H5T_NATIVE_INTEGER, Body%TRILinks, dims, h5err)
  call h5dclose_f(dset , h5err) 

!!$  do i=1,10
!!$     write(*,*) Body%Trilinks(1,i),Body%Trilinks(2,i)
!!$  end do
!  stop
  dims = (/4,Body%nLinks,1/)
  CALL h5dopen_f (file, "body/Bend_pts", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/Bend_pts')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%Bend_pts, dims, h5err)
  CALL h5dclose_f(dset , h5err) 
  
  dims = (/3,Body%nele,1/)
  CALL h5dopen_f (file, "body/Tri", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/Tri')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%Tri, dims, h5err)
  CALL h5dclose_f(dset , h5err)  
  
  dims = (/4,Body%nel,1/)
  CALL h5dopen_f (file, "body/IEN", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/IEN')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%IEN, dims, h5err)
  CALL h5dclose_f(dset , h5err)
!!$  do i=1,10
!!$     write(*,*) body%IEN(1,i), body%IEN(2,i), body%IEN(3,i), body%IEN(4,i)
!!$  end do
  !stop
  dims = (/3,Body%nnp,1/)
  CALL h5dopen_f (file, "body/ID", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/ID')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%ID, dims, h5err)
  CALL h5dclose_f(dset , h5err)
!!$   do i=1,15
!!$     write(*,*) body%ID(1,i), body%ID(2,i), body%ID(3,i)
!!$  end do
  !stop
  !
  dims = (/Body%nel,1,1/)
  CALL h5dopen_f (file, "body/MatType", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/MatType')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%MatType, dims, h5err)
  CALL h5dclose_f(dset , h5err)
  
  
  dims = (/Body%nel,1,1/)
  CALL h5dopen_f (file, "body/eltype", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/eltype')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%elType, dims, h5err)
  CALL h5dclose_f(dset , h5err)

  CALL h5dopen_f (file, "body/new_e", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'body/new_e')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%new_e, dims, h5err)
  CALL h5dclose_f(dset , h5err)


  ! 
  !    Wet surface info
  ! 
   !
  ! Load Wet Surface info
  !
  ! ws_max_eltype
  dims = (/1,1,1/)
  CALL h5dopen_f (file, "WetSurface/max_eltype", dset, h5err)
  call sm_io_checkHdfErr(h5err, 'WetSurface/max_eltype')
  CALL h5dread_f(dset, H5T_NATIVE_INTEGER, Body%ws_max_eltype, dims, h5err)
  CALL h5dclose_f(dset , h5err)
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
  ! allocate other ws_* arrays
  allocate( body%ws_nXi (   body%ws_nel ), &
            body%ws_nEta(   body%ws_nel ), &
            body%ws_ptelem( body%ws_nel )   )



  write(*,'(A)') 'Loading the HDF5 file is complete.'
        
  !*******************************************************************
  !*******************************************************************


  !
  ! DOF storage
  !
  write(*,'(A)',advance='no') 'Allocating Main Arrays...'
  
  allocate( Body%Qs(Body%nnp*NDIM) )   
  allocate( Body%Hs(Body%nnp*NDIM) )  
  allocate( Body%Hs_pres(Body%nnp*NDIM) )  
  allocate( Body%Hs_visc(Body%nnp*NDIM) )  

  Body%Qs = 0.
  body%Hs = 0.
  body%Hs_pres = 0.
  body%Hs_visc = 0.

  allocate( Body%normals(Body%nele*NDIM),Body%centers(Body%nele*NDIM) )     
  Body%Normals = 0.;
  Body%Centers = 0.;

  
  allocate(Body%comm(NDIM))

  sp_modelType= Body%MatType(1+Body%Nele);
  Nele= Body%Nele;
  nLinks=Body%nLinks;

  select case (sp_modelType)

     case (MATERIAL_WLC_POW)
        allocate (Body%kp(CONSTANT_ONE)  ,Body%kp_M(CONSTANT_ONE))
        
     case (MATERIAL_WLC_CQ)

     case (MATERIAL_WLC_FENE)
        allocate (Body%ks_M(CONSTANT_ONE),Body%ks(CONSTANT_ONE)) 

     case (MATERIAL_FENE_POW)
        allocate (Body%kp(CONSTANT_ONE)  ,Body%kp_M(CONSTANT_ONE))
        allocate (Body%ks_M(CONSTANT_ONE),Body%ks(CONSTANT_ONE)) 
        allocate (Body%Lmax_M(CONSTANT_ONE),Body%Lo_M(CONSTANT_ONE))
         allocate (Body%Lo(CONSTANT_ONE))


     case (MATERIAL_SF_FENE_POW) 
        allocate (Body%Dist(nLinks),Body%Lo_M(nLinks))
        allocate (Body%ks_M(CONSTANT_ONE),Body%ks(CONSTANT_ONE)) 
        allocate (Body%kp(nLinks)  ,Body%kp_M(nLinks))
      
     case (MATERIAL_SF_WLC_POW)
        
        allocate (Body%Dist(nLinks),Body%Lo_M(nLinks))
        allocate (Body%kp(nLinks)  ,Body%kp_M(nLinks))
        allocate (Body%Lmax_M(nLinks))
                
     case (MATERIAL_SF_WLC_CQ)
        allocate (Body%Dist(nLinks),Body%Lo_M(nLinks))
        
     case (MATERIAL_SF_WLC_FENE)
        allocate (Body%Dist(nLinks),Body%Lo_M(nLinks))
        allocate (Body%ks_M(nLinks),Body%ks(nLinks)) 

     case default
        call Driver_abortFlash(" No defined spring model was selected. &
        &  Check SolidMechanics.h and the MATERIAL Type section" )
     end select


     sp_modelType= Body%MatType(1);
     select case (sp_modelType) 
        
     case(MATERIAL_SF_AREA_VOL, MATERIAL_AREA_VOL)
        allocate(body%Aoj(Nele))
      
        !case(MATERIAL_AREA, MATERIAL_VOL)
     case default
        call Driver_abortFlash(" No defined area and volume  model was selected. &
             &  Check SolidMechanics.h and the MATERIAL Type section" )
     end select
     allocate(body%areas(Nele))
          
     write(*,*)'Now calculating the rbc parameters'
     
     call sm_rbc_Calc_parameters(ibd )
     write(*,*)'Done calculating the rbc parameters'
     ! initialize all the body masters to MASTERPE
     ! This body belongs to me:
     Body%BodyMaster = sm_MeshMe
  
     
     body%IntegMethod=3;   
     write(*,*) 'Stretching_exp',stretching_exp
     if (stretching_exp.eq.1) then
        CALL sortx(body%x,indx,body%nnp)
        CALL sorty(body%y,indy,body%nnp)
        !Getting the vertices to calculate the axial and transverse
              Npmax=indx(body%nnp)
              WRITE(*,*)' '
              WRITE(*,*)' '
              write(*,*) 'Maximum x coordinate is ',Npmax

              Npmin=indx(1)
              write(*,*) 'Minimum x coordinate is ',Npmin

              Ypmax=indy(body%nnp)
              write(*,*) 'Maximum y coordinate is ',Ypmax
              write(*,*) 'Percent of nodes=',per
              ExtPts=ceiling(per*real(body%nnp))
              write(*,*)'Number of points where the external force will be applied',ExtPts

              ALLOCATE(Nplus(ExtPts),Nminus(ExtPts))

              Nminus=indx(1:ExtPts)
              Nplus=indx(body%nnp-ExtPts+1:body%nnp)

              write(*,*) 'Nplus  is : ',Nplus
              write(*,*) 'Nminus is: ',Nminus
              Fext_M=Fext*1.e-12/(N_M*real(ExtPts));
        
        
     end if     
     write(*,*) 'The parameters for Solidbody ',ibd,' are calculated...'
     !stop 
   end subroutine sm_ioRead_rbc
   
   


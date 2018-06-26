! Subroutine outtotecplot
!
! Subroutine to write out to Paraview in hdf5 format.
!
! ---------------------------------------------------------------------------


  subroutine outtoParaview(filename,ib,mype,time,dt,istep,count, &
                          timer,blockList,blockCount,firstfileflag)
#include "constants.h"
#include "Flash.h"

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
                 Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
                 Grid_getBlkBoundBox,Grid_getBlkCenterCoords

  use ins_interface, only : ins_velgradtensor
  use HDF5
  use SolidMechanics_data, ONLY : sm_structure, sm_bodyInfo, sm_NumBodies, sm_meshMe,sm_meshComm
   
#ifdef FLASH_GRID_PARAMESH
  use physicaldata, only : interp_mask_unk, interp_mask_unk_res
#endif
  implicit none
#include "Flash_mpi.h"
  integer, intent(in)   :: mype,istep,count,firstfileflag,ib
  integer, intent(in)   :: blockCount
  integer, intent(in)   :: blockList(MAXBLOCKS)
  real, intent(in)      :: time,dt,timer
  character(len=100),intent(in)  :: filename

  ! Local variables    
  integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc
  character(6) :: index_lb,index_mype

  real:: xedge(NXB+1),xcell(NXB+1)
  real:: yedge(NYB+1),ycell(NYB+1)
  real:: zedge(NZB+1),zcell(NZB+1)
  real:: intsx(NXB+1), intsy(NYB+1), intsz(NZB+1)
  real(4) :: sxedge(NXB+1),syedge(NYB+1),szedge(NZB+1)
  
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  
  real :: facevarxx(NXB+2*NGUARD+1,NYB+2*NGUARD,NZB+2*NGUARD),&
       facevaryy(NXB+2*NGUARD,NYB+2*NGUARD+1,NZB+2*NGUARD),&
       facevarzz(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD+1)

  real, dimension(NXB+1,NYB+1,NZB+1) :: tpu,tpv,tpw,tpp,&
            tpdudxcorn, tpdudycorn, tpdudzcorn,& 
            tpdvdxcorn, tpdvdycorn, tpdvdzcorn,&
            tpdwdxcorn, tpdwdycorn, tpdwdzcorn,&
            vortx,vorty,vortz,omg,             &
            Sxy,Syz,Sxz,Oxy,Oyz,Oxz,Qcr,divpp,TVtpp
 
  real, dimension(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD) :: tpdudxc,&
        tpdudyc,tpdudzc,tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,&
        tpdwdyc,tpdwdzc

  integer ::  blockID,bodycounter

  real ::  del(3),dx,dy,dz
  real, dimension(MDIM)  :: coord,bsize
  real ::  boundBox(2,MDIM)  
  
  integer :: h5err, h5_outfile_id,group_id
  integer(HID_T) :: dset_id, space_id
  integer(HSIZE_T), DIMENSION(3) :: dims ! size read/write buffer
  character(len=100) :: filenameh5,bodyname,gridName,filenamexmf
  type(sm_structure),pointer :: body
  
  integer :: Nele, Nnodes, BodyType,ibd
  real,allocatable :: uVel(:),vVel(:),wVel(:);
!-----------------------------------------------------------------------

  body => sm_BodyInfo(ib)
  Nele=body%nele;
  Nnodes=body%nnp;
  allocate(uVel(Nnodes),vVel(Nnodes),wVel(Nnodes))

  ! write solution data to data.XXXX.XXXX
  write(*,*) 
  write(filenameh5,'(A14,I4.4,".",I4.4,".h5")') TRIM(filename),count,mype 
  write(*,*) 'The hdf5 file is ',filenameh5
  ! -- data.XXXX.XX --
  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1
 
  ! open(unit=22,file=filename,status='replace')  

  intsx    = (/ (real(i), i=0,NXB) /)
  intsy    = (/ (real(i), i=0,NYB) /)
  intsz    = (/ (real(i), i=0,NZB) /)

  ! ---------------------------
  ! ----   Structure 
  ! ---------------------------
  if (sm_NumBodies.gt.0) then
     bodycounter=0;
     do ibd = 1,sm_NumBodies
        
        if (sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then
           bodycounter=bodycounter+1;
           if (bodycounter.eq.1) then 
              
              ! open  the XMF and HDF5 files
              write(filenamexmf,'(A14,I4.4,".",I4.4,".xmf")') TRIM(filename),count,mype
              open(unit=30,file=filenamexmf,STATUS='UNKNOWN')
              write(30,'(A)') '<?xml version="1.0" ?>'
              write(30,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
              write(30,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2">'
              write(30,'(A,f6.4,A)') '<Time TimeType="Single" Value="',time,'"/>'
              write(30,'(A)') '<Domain>'
              
              ! Open the HDF5 File
              call H5Fcreate_f(trim(filenameh5),H5F_ACC_TRUNC_F, h5_outfile_id, h5err)
              call sm_io_checkHDFerr(h5err, 'failure to open  file')
              group_id = h5_outfile_id 
           end if
           
           write(30,'(A,I3,A)') '<Grid Name="RBC[',bodycounter,']" GridType="Uniform">'
                
           ! Writing the connectivity
           ! XMF
           write(bodyname,'(A,I3.3)')'CONN_',bodycounter
           !write(*,'(A,I3.3)')'CONN_',bodycounter
           write(30,'(A,I5,A)') '<Topology TopologyType="Triangle" Dimensions="',Nele,'" BaseOffset="1">'
           write(30,'(A,I7,A)') '<DataItem Dimensions="',Nele*NDIM,'" NumberType="Int" Precision="4" Format="HDF">'
           write(30,'(A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem></Topology>'
           ! HDF5
           dims = (/Nele*NDIM,1,1/)
           call H5Screate_simple_f(1, dims, space_id, h5err)
           call H5Dcreate_f(group_id,TRIM(adjustl(bodyname)), H5T_NATIVE_INTEGER, space_id, dset_id, h5err)
           call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, body%Tri, dims, h5err)
           call H5Dclose_f(dset_id, h5err)
           call H5Sclose_f(space_id, h5err)
           ! Mesh
           !XMF
           write(bodyname,'(A,I3.3)')'RBC_',bodycounter
           write(30,'(A)') '<Geometry GeometryType="XYZ">'
           write(30,'(A,I7,A)') '<DataItem Dimensions="',Nnodes*NDIM,'" NumberType="Float" Precision="8" Format="HDF">'
           write(30,'(A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem></Geometry>'
           ! Writing the mesh
           ! HDF5 
           dims = (/Nnodes*NDIM,1,1/)
           call H5Screate_simple_f(1, dims, space_id, h5err)
           call H5Dcreate_f(group_id, TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
           call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, body%qms(:,1), dims, h5err)
           call H5Dclose_f(dset_id, h5err)
           call H5Sclose_f(space_id, h5err)
           ! Writing the velocities
           
           dims(:) = (/Nnodes,1,1/)
           do i=1,Nnodes
              uVel(i) = body%qdms(3*i-2,1);
              vVel(i) = body%qdms(3*i-1,1);
              wVel(i) = body%qdms(3*i  ,1);
           end do

           write(bodyname,'(A,I3.3)')'uVel_',bodycounter
           write(30,'(A)') '<Attribute AttributeType="Scalar" Center="Node"'  
           write(30,'(A,A,A)') 'Name="',TRIM(bodyname),'">'
           write(30,'(A,I6,A)') '<DataItem DataType="Float" Precision="8" Dimensions="',Nnodes,'  1" Format="HDF">'
           write(30,'(A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem></Attribute>'
           call H5Screate_simple_f(1, dims, space_id, h5err)
           call H5Dcreate_f(group_id, TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
           call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,uVel, dims, h5err)
           call H5Dclose_f(dset_id, h5err)
           call H5Sclose_f(space_id, h5err)
           
           write(bodyname,'(A,I3.3)')'vVel_',bodycounter
           write(30,'(A)') '<Attribute AttributeType="Scalar" Center="Node"'  
           write(30,'(A,A,A)') 'Name="',TRIM(bodyname),'">'
           write(30,'(A,I6,A)') '<DataItem DataType="Float" Precision="8" Dimensions="',Nnodes,'  1" Format="HDF">'
           write(30,'(A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem></Attribute>'
           !Nodes=body%qdms(2:3:Nnodes*3-1,1);
           call H5Screate_simple_f(1, dims, space_id, h5err)
           call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
           call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,vVel, dims, h5err)
           call H5Dclose_f(dset_id, h5err)
           call H5Sclose_f(space_id, h5err)
           !Nodes=body%qdms(3:3:Nnodes*3,1);
           write(bodyname,'(A,I3.3)')'wVel_',bodycounter
           write(30,'(A)') '<Attribute AttributeType="Scalar" Center="Node"'  
           write(30,'(A,A,A)') 'Name="',TRIM(bodyname),'">'
           write(30,'(A,I6,A)') '<DataItem DataType="Float" Precision="8" Dimensions="',Nnodes,'  1" Format="HDF">'
           write(30,'(A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem></Attribute>'
           write(30,'(A)') '</Grid>'
           call H5Screate_simple_f(1, dims, space_id, h5err)
           call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
           call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,wVel, dims, h5err)
           call H5Dclose_f(dset_id, h5err)
           call H5Sclose_f(space_id, h5err)
        end if
     end do
     !write(*,*)'bodycounter',bodycounter
     
     if (bodycounter.gt.0) then
        
        write(30,'(A)') '</Domain>'
        write(30,'(A)') '</Xdmf>'
        close(30)
        ! Close the group
        if( group_id /= h5_outfile_id ) then
           call H5Gclose_f (group_id, h5err)    
        end if
        ! close the file
        
        call H5Fclose_f(h5_outfile_id, h5err)    
        call sm_io_checkHDFerr(h5err,'failure to close HDF5 file')
     end if
  end if
  deallocate(uVel,vVel,wVel)

  ! -------------------------------
  ! ----  FLUID    
  !--------------------------------
  do lb = 1,blockcount
     
     blockID =  blockList(lb)      


     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)
     dx = del(IAXIS)
     dy = del(JAXIS)
     dz = del(KAXIS)

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     tpu = 0.
     tpv = 0.
     tpw = 0.
     tpp = 0.
     omg = 0.

     xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
     xcell = xedge(:) + dx/2.0;
    
     yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
     ycell = yedge(:) + dy/2.0;
    
     zedge = coord(KAXIS) - bsize(KAXIS)/2.0 + dz*intsz;
     zcell = zedge(:) + dz/2.0; 

     
     facevarxx = facexData(VELC_FACE_VAR,:,:,:)
     facevaryy = faceyData(VELC_FACE_VAR,:,:,:)
     facevarzz = facezData(VELC_FACE_VAR,:,:,:)          
 
     ! U velocity: u(NXB+1,NYB+1,NZB+1)
     ! --------------------------------
     tpu = 0.25*(facevarxx(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD:nzc-1)+ &
                 facevarxx(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc)+ &
                 facevarxx(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc)+ &
                 facevarxx(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1))
    
     ! V velocity: v(NXB+1,NYB+1,NZB+1)
     ! --------------------------------                           
     tpv = 0.25*(facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD:nzc-1) + &
                 facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD:nzc-1) + &
                 facevaryy(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc) + &
                 facevaryy(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));
    
!         

     ! W velocity: w(NXB+1,NYB+1,NZB+1)
     ! --------------------------------
     tpw = 0.25*(facevarzz(NGUARD:nxc-1,NGUARD:nyc-1,NGUARD+1:nzc)   + &
                 facevarzz(NGUARD+1:nxc,NGUARD:nyc-1,NGUARD+1:nzc)   + &
                 facevarzz(NGUARD+1:nxc,NGUARD+1:nyc,NGUARD+1:nzc)   + &
                 facevarzz(NGUARD:nxc-1,NGUARD+1:nyc,NGUARD+1:nzc));
      
     ! P pressure: p(NXB+1,NYB+1,NZB+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc, &
                             solnData(PRES_VAR,:,:,:),tpp)

     ! TV : TV(NXB+1,NYB+1,NZB+1)
     ! -------------------------------
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             solnData(TVIS_VAR,:,:,:),TVtpp)


     ! Divergence: 
     ! ----------
     solnData(DUST_VAR,:,:,:) = 0.
     solnData(DUST_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1) =      &
             (facevarxx(NGUARD+2:nxc  ,NGUARD+1:nyc-1,NGUARD+1:nzc-1)     - &
              facevarxx(NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/dx + &
             (facevaryy(NGUARD+1:nxc-1,NGUARD+2:nyc  ,NGUARD+1:nzc-1)     - &
              facevaryy(NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/dy + &
             (facevarzz(NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+2:nzc  )     - &
              facevarzz(NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/dz 
    

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             solnData(DUST_VAR,:,:,:),divpp)


     ! Velocity derivatives:
     ! -------- -----------            
     call ins_velgradtensor(NGUARD,facexData,faceyData,facezData, &
                               dx,dy,dz,tpdudxc,tpdudyc,tpdudzc,&
                tpdvdxc,tpdvdyc,tpdvdzc,tpdwdxc,tpdwdyc,tpdwdzc)

     ! Extrapolation of center derivatives to corners, the values
     ! of derivatives in guardcells next to edges are obtained 
     ! from real velocities and linearly extrapolated velocities
     ! to edge points.

     ! U derivatives:
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdudxc,tpdudxcorn)
            
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdudyc,tpdudycorn)

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdudzc,tpdudzcorn)

     ! V derivatives:
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdvdxc,tpdvdxcorn)
            
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdvdyc,tpdvdycorn)

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdvdzc,tpdvdzcorn)


     ! W derivatives:
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdwdxc,tpdwdxcorn)
            
     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdwdyc,tpdwdycorn)

     call centervals2corners(NGUARD,NXB,NYB,NZB,nxc,nyc,nzc,&
                             tpdwdzc,tpdwdzcorn)
         
     ! VORTICITY:
     ! ---------
     ! Corner values of vorticity:
     vortx = tpdwdycorn - tpdvdzcorn
     vorty = tpdudzcorn - tpdwdxcorn
     vortz = tpdvdxcorn - tpdudycorn

     ! Q criterion:
     ! First Sij:
     !Sxx = tpdudxcorn
     Sxy = 0.5*(tpdudycorn + tpdvdxcorn)
     Sxz = 0.5*(tpdudzcorn + tpdwdxcorn)
     !Syy = tpdvdycorn
     Syz = 0.5*(tpdvdzcorn + tpdwdycorn)
     !Szz = tpdwdzcorn

     ! Then Oij:
     Oxy = 0.5*(tpdudycorn - tpdvdxcorn)
     Oxz = 0.5*(tpdudzcorn - tpdwdxcorn)
     Oyz = 0.5*(tpdvdzcorn - tpdwdycorn)

     ! Calculate Q:
     Qcr = 0.5*(2.*(Oxy*Oxy+Oyz*Oyz+Oxz*Oxz) -                     &
               (tpdudxcorn*tpdudxcorn + tpdvdycorn*tpdvdycorn +    &
                tpdwdzcorn*tpdwdzcorn +                            &
                2.*(Sxy*Sxy+Syz*Syz+Sxz*Sxz)))                      

     !---------   WRITE THE HDF5 FILE -----------
     ! ------------------------------------------
     ! create file name
     !write(*,*) 'lb', lb
     if (lb.eq.1) then
        ! HDF5
        write(filenameh5,'("data.",I4.4,".",I4.4,".h5")') count,mype
        call H5Fcreate_f(trim(adjustl(filenameh5)),H5F_ACC_TRUNC_F, h5_outfile_id, h5err)
        group_id = h5_outfile_id
        
        ! XDMF
        write(filenamexmf,'("data.",I4.4,".",I4.4,".xmf")') count,mype
        open(unit=30,file=filenamexmf,STATUS='UNKNOWN')
        write(30,'(A)') '<?xml version="1.0" ?>'
        write(30,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
        write(30,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2">'
        write(30,'(A,I7,A)')'<Information Name="Iteration" Value="',istep,'"/>'
        write(30,'(A,f6.4,A)') '<Time TimeType="Single" Value="',time,'"/>'
        write(30,'(A)') '<Domain>'
        write(30,'(A)') '<Grid GridType="Collection" CollectionType="Spatial">'
     end if
     
    
     write(30,'(A,I3,A)') '<Grid Name="Fluidmesh[',lb-1,']">'
     write(30,'(A,I4,I4,I4,A)') '<Topology TopologyType="3DRectMesh" NumberOfElements="',NZB+1,NYB+1,NXB+1,'"/>'
     write(30,'(A)') '<Geometry GeometryType="VXVYVZ">'
     ! Write the x
     sxedge = xedge
     syedge = yedge
     szedge = zedge
     dims=(/(NXB+1),1,1/)
     write(bodyname,'("X_",I3.3)')lb
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_REAL, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, sxedge, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  call H5Sclose_f(space_id, h5err)
  write(30,'(A,I4,A)') '  <DataItem Dimensions="',(NXB+1),'" NumberType="Float" Precision="4" Format="HDF">'
  write(30,'(A,A,A,A,A)') '  ',TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
  
  write(bodyname,'("Y_",I3.3)')lb
  dims=(/(NYB+1),1,1/)
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_REAL, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, syedge, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  write(30,'(A,I4,A)') '<DataItem Dimensions="',(NYB+1),'" NumberType="Float" Precision="4" Format="HDF">'
  write(30,'(A,A,A,A,A)')'   ', TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
  
  ! write the z
  dims=(/(NZB+1),1,1/)
  write(bodyname,'("Z_",I3.3)')lb
  call H5Screate_simple_f(1, dims, space_id, h5err)
  call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_REAL, space_id, dset_id, h5err)
  call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, szedge, dims, h5err)
  call H5Dclose_f(dset_id, h5err)
  write(30,'(A,I4,A)') '  <DataItem Dimensions="',NZB+1,'" NumberType="Float" Precision="4" Format="HDF">'
  write(30,'(A,A,A,A)')' ',TRIM(filenameh5),':/',TRIM(bodyname)
  write(30,'(A)')'  </DataItem>'
  write(30,'(A)')'</Geometry>'
  
  ! Writing the velocities
  ! --------- u: ----------
  dims (1) = (NXB+1)
  dims (2) = (NYB+1)
  dims (3) = (NZB+1)
  
   write(bodyname,'("U_",I3.3)')lb
   call H5Screate_simple_f(3, dims, space_id, h5err)
   call H5Dcreate_f(group_id, TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
   call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tpu, dims, h5err)
   call H5Dclose_f(dset_id, h5err)
   write(30,'(A,A,A)') '<Attribute Name="','U','" AttributeType="Scalar" Center="Node">'
   write(30,'(A,I4,I4,I4,A)') '<DataItem Dimensions="',NZB+1,NYB+1,(NXB+1),'" NumberType="Float" Precision="8" Format="HDF">'
   write(30,'(A,A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
   write(30,'(A)') '</Attribute>'

     !---  v: -----
     write(bodyname,'("V_",I3.3)')lb
     call H5Screate_simple_f(3, dims, space_id, h5err)
     call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tpv, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
     write(30,'(A,A,A)') '<Attribute Name="','V','" AttributeType="Scalar" Center="Node">'
     write(30,'(A,I4,I4,I4,A)') '<DataItem Dimensions="',NZB+1,NYB+1,(NXB+1),'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
     write(30,'(A)') '</Attribute>'
     ! ------- w: -------
     write(bodyname,'("W_",I3.3)')lb
     call H5Screate_simple_f(3, dims, space_id, h5err)
     call H5Dcreate_f(group_id, TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tpw, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
     write(30,'(A,A,A)') '<Attribute Name="','W','" AttributeType="Scalar" Center="Node">'
     write(30,'(A,I4,I4,I4,A)') '<DataItem Dimensions="',NZB+1,NYB+1,(NXB+1),'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
     write(30,'(A)') '</Attribute>'
     !-----  p:
     write(bodyname,'("P_",I3.3)')lb
     call H5Screate_simple_f(3, dims, space_id, h5err)
     call H5Dcreate_f(group_id, TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tpp, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
     write(30,'(A,A,A)') '<Attribute Name="','P','" AttributeType="Scalar" Center="Node">'
     write(30,'(A,I4,I4,I4,A)') '<DataItem Dimensions="',NZB+1,NYB+1,(NXB+1),'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
     write(30,'(A)') '</Attribute>'
     ! ---- omgX:
     write(bodyname,'("omgX_",I3.3)')lb
     call H5Screate_simple_f(3, dims, space_id, h5err)
     call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vortx, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
     write(30,'(A,A,A)') '<Attribute Name="','OmgX','" AttributeType="Scalar" Center="Node">'
     write(30,'(A,I4,I4,I4,A)') '<DataItem Dimensions="',NZB+1,NYB+1,(NXB+1),'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
     write(30,'(A)') '</Attribute>'
     ! ----- omgY:
     write(bodyname,'("omgY_",I3.3)')lb
     call H5Screate_simple_f(3, dims, space_id, h5err)
     call H5Dcreate_f(group_id, TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vorty, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
     write(30,'(A,A,A)') '<Attribute Name="','OmgY','" AttributeType="Scalar" Center="Node">'
     write(30,'(A,I4,I4,I4,A)') '<DataItem Dimensions="',NZB+1,NYB+1,(NXB+1),'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
     write(30,'(A)') '</Attribute>'
     ! ----  omgZ:
     write(bodyname,'("omgZ_",I3.3)')lb
     call H5Screate_simple_f(3, dims, space_id, h5err)
     call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vortz, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
     write(30,'(A,A,A)') '<Attribute Name="','OmgZ','" AttributeType="Scalar" Center="Node">'
     write(30,'(A,I4,I4,I4,A)') '<DataItem Dimensions="',NZB+1,NYB+1,(NXB+1),'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
     write(30,'(A)') '</Attribute>'
     !   
     ! ----- Q:----- the velocity gradient varient
     write(bodyname,'("Qcr_",I3.3)')lb
     call H5Screate_simple_f(3, dims, space_id, h5err)
     call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Qcr, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
    write(30,'(A,A,A)') '<Attribute Name="','Qcr','" AttributeType="Scalar" Center="Node">'
     write(30,'(A,I4,I4,I4,A)') '<DataItem Dimensions="',NZB+1,NYB+1,(NXB+1),'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
     write(30,'(A)') '</Attribute>'
     ! ---- Divergance : -----
     write(bodyname,'("divpp_",I3.3)')lb
     call H5Screate_simple_f(3, dims, space_id, h5err)
     call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, divpp, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
     write(30,'(A,A,A)') '<Attribute Name="','Divpp','" AttributeType="Scalar" Center="Node">'
     write(30,'(A,I4,I4,I4,A)') '<DataItem Dimensions="',NZB+1,NYB+1,(NXB+1),'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
     write(30,'(A)') '</Attribute>'
     ! ---- Turbulent viscosity
     write(bodyname,'("TVtpp_",I3.3)')lb
     call H5Screate_simple_f(3, dims, space_id, h5err)
     call H5Dcreate_f(group_id,TRIM(bodyname), H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, TVtpp, dims, h5err)
     call H5Dclose_f(dset_id, h5err)
     write(30,'(A,A,A)') '<Attribute Name="','TVtpp','" AttributeType="Scalar" Center="Node">'
     write(30,'(A,I4,I4,I4,A)') '<DataItem Dimensions="',NZB+1,NYB+1,(NXB+1),'" NumberType="Float" Precision="8" Format="HDF">'
     write(30,'(A,A,A,A)') TRIM(filenameh5),':/',TRIM(bodyname),'</DataItem>'
     write(30,'(A)') '</Attribute>'
     ! Close the grid in the xmf file
     write(30,'(A)') '</Grid>'
  end do
  

  ! Close the group
  if (blockcount.ge.1) then
     write(30,'(A)') '</Grid>'
     write(30,'(A)') '</Domain>'
     write(30,'(A)') '</Xdmf>'
     close(30)
     if( group_id /= h5_outfile_id ) then
        call H5Gclose_f (group_id, h5err)    
     end if
     
     !
     ! close the file
     !
     call H5Fclose_f(h5_outfile_id, h5err)    
     call sm_io_checkHDFerr(h5err,'failure to close HDF5 file')
  end if
  if (mype .eq. 0) then
     write(*,*) ''
     write(*,*) '*** Wrote data to HDF5 file to ',filename,' ****'
  endif
  
  return
     
  
End subroutine outtoParaview

   
   
   ! Subroutine centervals2corners:
! Subroutine to obtain corver values of a variable given the center 
! values of it in a 3D structured block, suppossing guardcells already
! filled.
!
! ----------------------------------------------------------------------

  subroutine centervals2corners(ng,nxb,nyb,nzb,nxc,nyc,nzc, &
                                unk1,tpp)

    implicit none

    integer ng,nxb,nyb,nzb,nxc,nyc,nzc
    integer nx1,ny1,nz1,nx2,ny2,nz2
    real*8, intent(in) :: unk1(nxb+2*ng,nyb+2*ng,nzb+2*ng)
    real*8, intent(out) :: tpp(nxb+1,nyb+1,nzb+1)

      
    tpp = .5*.25*(unk1(ng:nxc-1,ng:nyc-1,ng:nzc-1) + &
                  unk1(ng:nxc-1,ng:nyc-1,ng+1:nzc) + &
                  unk1(ng:nxc-1,ng+1:nyc,ng:nzc-1) + &
                  unk1(ng+1:nxc,ng:nyc-1,ng:nzc-1) + &
                  unk1(ng+1:nxc,ng+1:nyc,ng:nzc-1) + &
                  unk1(ng:nxc-1,ng+1:nyc,ng+1:nzc) + &
                  unk1(ng+1:nxc,ng:nyc-1,ng+1:nzc) + &
                  unk1(ng+1:nxc,ng+1:nyc,ng+1:nzc));

    ! Z edges:
    ! Edge: x = 1, y = 1:
    tpp(1,1,:) = .25*(unk1(ng,ng+1,ng:nzc-1) + &
                      unk1(ng,ng+1,ng+1:nzc) + &
                      unk1(ng+1,ng,ng:nzc-1) + &
                      unk1(ng+1,ng,ng+1:nzc));           
            
    ! Edge: x = 1, y = end:
    tpp(1,nyb+1,:) = .25*(unk1(ng,nyc-1,ng:nzc-1) + & 
                          unk1(ng,nyc-1,ng+1:nzc) + &
                          unk1(ng+1,nyc,ng:nzc-1) + &
                          unk1(ng+1,nyc,ng+1:nzc));    

    ! Edge: x = end, y = 1:
    tpp(nxb+1,1,:) = .25*(unk1(nxc-1,ng,ng:nzc-1) + &
                          unk1(nxc-1,ng,ng+1:nzc) + &
                          unk1(nxc,ng+1,ng:nzc-1) + &
                          unk1(nxc,ng+1,ng+1:nzc));

    ! Edge: x = end, y = end:
    tpp(nxb+1,nyb+1,:) = .25*(unk1(nxc-1,nyc,ng:nzc-1) + &
                              unk1(nxc-1,nyc,ng+1:nzc) + &
                              unk1(nxc,nyc-1,ng:nzc-1) + &
                              unk1(nxc,nyc-1,ng+1:nzc));
            
    ! Y edges:
    ! Edge: x = 1, z = 1:
    tpp(1,:,1) = .25*(unk1(ng,ng:nyc-1,ng+1) + &
                      unk1(ng,ng+1:nyc,ng+1) + &
                      unk1(ng+1,ng:nyc-1,ng) + &
                      unk1(ng+1,ng+1:nyc,ng));

    ! Edge: x = 1, z = end:
    tpp(1,:,nzb+1) = .25*(unk1(ng,ng:nyc-1,nzc-1) + &
                          unk1(ng,ng+1:nyc,nzc-1) + &
                          unk1(ng+1,ng:nyc-1,nzc) + &
                          unk1(ng+1,ng+1:nyc,nzc)); 

    ! Edge: x = end, z = 1:
    tpp(nxb+1,:,1) = .25*(unk1(nxc-1,ng:nyc-1,ng) + &
                          unk1(nxc-1,ng+1:nyc,ng) + &
                          unk1(nxc,ng:nyc-1,ng+1) + &
                          unk1(nxc,ng+1:nyc,ng+1));

    ! Edge: x = end, z = end:
    tpp(nxb+1,:,nzb+1) = .25*(unk1(nxc-1,ng:nyc-1,nzc) + &
                              unk1(nxc-1,ng+1:nyc,nzc) + &
                              unk1(nxc,ng:nyc-1,nzc-1) + &
                              unk1(nxc,ng+1:nyc,nzc-1));

    ! X edges:
    ! Edge: y = 1, z = 1:
    tpp(:,1,1) = .25*(unk1(ng:nxc-1,ng,ng+1) + &
                      unk1(ng+1:nxc,ng,ng+1) + &
                      unk1(ng:nxc-1,ng+1,ng) + &
                      unk1(ng+1:nxc,ng+1,ng));            

    ! Edge: y = 1, z = end:
    tpp(:,1,nzb+1) = .25*(unk1(ng:nxc-1,ng,nzc-1) + &
                          unk1(ng+1:nxc,ng,nzc-1) + &
                          unk1(ng:nxc-1,ng+1,nzc) + &
                          unk1(ng+1:nxc,ng+1,nzc));   

    ! Edge: y = end, z = 1:
    tpp(:,nyb+1,1) = .25*(unk1(ng:nxc-1,nyc-1,ng) + & 
                          unk1(ng+1:nxc,nyc-1,ng) + &
                          unk1(ng:nxc-1,nyc,ng+1) + &
                          unk1(ng+1:nxc,nyc,ng+1)); 

    ! Edge: y = end, z = end:
    tpp(:,nyb+1,nzb+1) = .25*(unk1(ng:nxc-1,nyc-1,nzc) + &
                              unk1(ng+1:nxc,nyc-1,nzc) + &
                              unk1(ng:nxc-1,nyc,nzc-1) + &
                              unk1(ng+1:nxc,nyc,nzc-1)); 

    ! Corners
    ! Corner x = 1, y = 1, z = 1:
    tpp(1,1,1) = -0.5*unk1(ng+1,ng+1,ng+1) + &
                  0.5*unk1(ng,ng+1,ng+1)   + &
                  0.5*unk1(ng+1,ng,ng+1)   + &
                  0.5*unk1(ng+1,ng+1,ng)   


    ! Corner x = end, y =1, z = 1:
    tpp(nxb+1,1,1) = -0.5*unk1(nxc-1,ng+1,ng+1) + &
                      0.5*unk1(nxc,ng+1,ng+1)   + &
                      0.5*unk1(nxc-1,ng,ng+1)   + &
                      0.5*unk1(nxc-1,ng+1,ng)   


    ! Corner x = end, y = end, z = 1:
    tpp(nxb+1,nyb+1,1) = -0.5*unk1(nxc-1,nyc-1,ng+1) + &
                          0.5*unk1(nxc,nyc-1,ng+1)   + &
                          0.5*unk1(nxc-1,nyc,ng+1)   + &
                          0.5*unk1(nxc-1,nyc-1,ng)

    ! Corner x = 1, y = end, z = 1:
    tpp(1,nyb+1,1) = -0.5*unk1(ng+1,nyc-1,ng+1) + &
                      0.5*unk1(ng,nyc-1,ng+1)   + &
                      0.5*unk1(ng+1,nyc,ng+1)   + &
                      0.5*unk1(ng+1,nyc-1,ng)   

    ! Corner x = 1, y = 1, z = end:
    tpp(1,1,nzb+1) = -0.5*unk1(ng+1,ng+1,nzc-1) + &
                      0.5*unk1(ng,ng+1,nzc-1)   + &
                      0.5*unk1(ng+1,ng,nzc-1)   + &
                      0.5*unk1(ng+1,ng+1,nzc)               

    ! Corner x = end, y = 1, z = end:
    tpp(nxb+1,1,nzb+1) = -0.5*unk1(nxc-1,ng+1,nzc-1) + &
                          0.5*unk1(nxc,ng+1,nzc-1)   + &
                          0.5*unk1(nxc-1,ng,nzc-1)   + &
                          0.5*unk1(nxc-1,ng+1,nzc)    

    ! Corner x = end, y = end, z = end:
    tpp(nxb+1,nyb+1,nzb+1) = -0.5*unk1(nxc-1,nyc-1,nzc-1) + &
                              0.5*unk1(nxc,nyc-1,nzc-1)   + &
                              0.5*unk1(nxc-1,nyc,nzc-1)   + &
                              0.5*unk1(nxc-1,nyc-1,nzc)

    ! Corner x = 1, y = end, z = end:
    tpp(1,nyb+1,nzb+1) = -0.5*unk1(ng+1,nyc-1,nzc-1) + &
                          0.5*unk1(ng,nyc-1,nzc-1)   + &
                          0.5*unk1(ng+1,nyc,nzc-1)   + &
                          0.5*unk1(ng+1,nyc-1,nzc)          

  End subroutine centervals2corners

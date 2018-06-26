subroutine sm_xdmfwrite(filebase,ibd,mype,time,dt,istep,count, &
                          timer,firstfileflag) 

#include "Flash.h"
#include "SolidMechanics.h"
  use Driver_interface, Only : Driver_abortFlash,Driver_getSimTime
  use sm_iointerface, only: sm_io_checkHdfErr,sm_xdmfwrite
  use SolidMechanics_data, ONLY : sm_structure, sm_bodyInfo, sm_NumBodies
  use HDF5
  
  implicit none
  
  ! Argument list
  character(len=100),intent (IN) :: filebase
  integer, intent(IN) :: ibd
  integer, intent(IN) :: mype,istep,count,firstfileflag
  real, intent(IN)  :: time,timer,dt
  ! Local variables
  type(sm_structure),pointer :: body  
  integer :: Nele, Nnodes, BodyType,i
  integer :: h5err, h5_outfile_id,group_id
  integer(HID_T) :: dset_id, space_id
  integer(HSIZE_T), DIMENSION(3) :: dims ! size read/write buffer
  character(len=100) :: filenameh5,filenamexmf
  

  ! -------- OPEN THE XMF FILE 
  ! write the xmf files
  write(filenamexmf,'(A14,I4.4,".",I4.4,".xmf")') TRIM(filebase),count,mype
  write(filenameh5,'(A14,I4.4,".",I4.4,".h5")') TRIM(filebase),count,mype 
     
     
     open(unit=30,file=filenamexmf,STATUS='UNKNOWN')
     write(30,'(A)') '<?xml version="1.0" ?>'
     write(30,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
     write(30,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2">'
     write(30,'(A)') '<Domain>'
     
     if (ibd.gt.0) then
        body => sm_BodyInfo(ibd)
        BodyType=body%bodytype;
        
        Nele=body%nele;
        Nnodes=body%nnp;
        
        write(30,'(A)') '<Grid Name="RBC Grid" GridType="Uniform">'
        write(30,'(A,I5,A)') '<Topology TopologyType="Triangle" Dimensions="',Nele,'" BaseOffset="1">'
        write(30,'(A,I7,A,A,A)') '<DataItem Dimensions="',Nele*NDIM,'" NumberType="Int" Precision="4" Format="HDF">'
        write(30,'(A,A)') TRIM(filenameh5),':/Connectivity</DataItem></Topology>'
        write(30,'(A)') '<Geometry GeometryType="XYZ">'
        write(30,'(A,I7,A,A,A)') '<DataItem Dimensions="',Nnodes*NDIM,'" NumberType="Float" Precision="8" Format="HDF">',TRIM(filenameh5),':/mesh</DataItem></Geometry>'
        write(30,'(A)') '<Attribute AttributeType="Scalar" Center="Node"'  
        write(30,'(A)') 'Name="UVel">'
        write(30,'(A,I6,A)') '<DataItem DataType="Float" Precision="8" Dimensions="',Nnodes,'  1" Format="HDF">'
        write(30,'(A,A)') TRIM(filenameh5),':/uVel</DataItem></Attribute>'
        write(30,'(A)') '<Attribute AttributeType="Scalar" Center="Node"'  
        write(30,'(A)') 'Name="VVel">'
        write(30,'(A,I6,A)') '<DataItem DataType="Float" Precision="8" Dimensions="',Nnodes,'  1" Format="HDF">'
        write(30,'(A,A)') TRIM(filenameh5),':/vVel</DataItem></Attribute>'
        write(30,'(A)') '<Attribute AttributeType="Scalar" Center="Node"'  
        write(30,'(A)') 'Name="WVel">'
        write(30,'(A,I6,A)') '<DataItem DataType="Float" Precision="8" Dimensions="',Nnodes,'  1" Format="HDF">'
        write(30,'(A,A)') TRIM(filenameh5),':/wVel</DataItem></Attribute>'
     end if
     
     write(30,'(A,f6.4,A)') '<Time TimeType="Single" Value="',time,'"/>'
     write(30,'(A)') '</Grid>'
     write(30,'(A)') '<Grid Name="Fluid mesh" GridType="Uniform">'
     write(30,'(A,I4,I4,I4,A)') '<Topology TopologyType="3DRectMesh" NumberOfElements="',NXB+1,NYB+1,NZB+1,'"/>'
     write(30,'(A)') '<Geometry GeometryType="VXVYVZ">'
     write(30,'(A,I4,A)') '<DataItem Dimensions="',NXB+1,'" NumberType="Float" Precision="4" Format="HDF">'
     write(30,'(A,A)') TRIM(filenameh5),':/x  </DataItem>'
     write(30,'(A,I4,A)') '<DataItem Dimensions="',NYB+1,'" NumberType="Float" Precision="4" Format="HDF">'
     write(30,'(A,A)') TRIM(filenameh5),':/y  </DataItem>'
     write(30,'(A,I4,A)') '<DataItem Dimensions="',NZB+1,'" NumberType="Float" Precision="4" Format="HDF">'
     write(30,'(A,A)') TRIM(filenameh5),':/z  </DataItem>'
     write(30,'(A)') '</Geometry>'
     write(30,'(A)') '<Attribute Name="u" AttributeType="Scalar" Center="Node">'
     write(30,'(A,I4,I4,I4,A)') '<DataItem Dimensions="',NXB+1,NYB+1,NZB+1,'," NumberType="Float" Precision="4" Format="HDF">'
     write(30,'(A)') TRIM(filenameh5),':/u  </DataItem></Attribute>'
     write(30,'(A)') '</Grid>'
     write(30,'(A)') '</Domain>'
     write(30,'(A)') '</Xdmf>'
     close(30)
     
  
end subroutine sm_xdmfwrite


!!$     !write the hdf5 files
!!$     write(filenameh5,'(A14,I9.9,A3)') TRIM(filebase),dr_nstep,'.h5'
!!$     call H5Fcreate_f(trim(filenameh5),H5F_ACC_TRUNC_F, h5_outfile_id, h5err)
!!$     call sm_io_checkHDFerr(h5err, 'failure to open  file')
!!$     
!!$     ! Change this if you want to write to a specific group
!!$     group_id = h5_outfile_id    
!!$     dims = (/Nnodes*NDIM,1,1/)
!!$     ! Writing the mesh
!!$     call H5Screate_simple_f(1, dims, space_id, h5err)
!!$     call H5Dcreate_f(group_id, 'mesh', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
!!$     call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, body%qms(:,1), dims, h5err)
!!$     call H5Dclose_f(dset_id, h5err)
!!$     ! Writing the connectivity
!!$      dims = (/Nele*NDIM,1,1/)
!!$      call H5Screate_simple_f(1, dims, space_id, h5err)
!!$      call H5Dcreate_f(group_id, 'Connectivity', H5T_NATIVE_INTEGER, space_id, dset_id, h5err)
!!$      call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, body%Tri, dims, h5err)
!!$      call H5Dclose_f(dset_id, h5err)
!!$      call H5Sclose_f(space_id, h5err)
!!$      
!!$      ! Writing the velocities
!!$      dims(1:2) = (/Nnodes,1/)
!!$      do i=1,Nnodes
!!$         u(i) = body%qdms(3*i-2,1);
!!$         v(i) = body%qdms(3*i-1,1);
!!$         w(i) = body%qdms(3*i  ,1);
!!$      end do
!!$      
!!$      call H5Screate_simple_f(1, dims, space_id, h5err)
!!$      call H5Dcreate_f(group_id, 'uVel', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
!!$      call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,u, dims, h5err)
!!$      call H5Dclose_f(dset_id, h5err)
!!$      call H5Sclose_f(space_id, h5err)
!!$      !Nodes=body%qdms(2:3:Nnodes*3-1,1);
!!$      call H5Screate_simple_f(1, dims, space_id, h5err)
!!$      call H5Dcreate_f(group_id, 'vVel', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
!!$      call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,v, dims, h5err)
!!$      call H5Dclose_f(dset_id, h5err)
!!$      call H5Sclose_f(space_id, h5err)
!!$      !Nodes=body%qdms(3:3:Nnodes*3,1);
!!$      call H5Screate_simple_f(1, dims, space_id, h5err)
!!$      call H5Dcreate_f(group_id, 'wVel', H5T_NATIVE_DOUBLE, space_id, dset_id, h5err)
!!$      call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,w, dims, h5err)
!!$      call H5Dclose_f(dset_id, h5err)
!!$      call H5Sclose_f(space_id, h5err)
!!$   
!!$     ! Close the group
!!$     if( group_id /= h5_outfile_id ) then
!!$        call H5Gclose_f (group_id, h5err)    
!!$     end if
!!$     
!!$     !
!!$     ! close the file
!!$     !
!!$     call H5Fclose_f(h5_outfile_id, h5err)    
!!$     call sm_io_checkHDFerr(h5err,'failure to close snapshots file')


!!****if* source/physics/SolidMechanics/SolidMechanicsMain/IO/sm_uvTecplot
!!
!! NAME
!! 
!!
!! SYNOPSIS
!!
!!  
!! DESCRIPTION 
!!
!! File written in ASCI format, for speedup of Tecplot load use preplot command. 
!! i.e.: preplot geo.000001.plt geobin.000001.plt
!!
!! ARGUMENTS 
!!
!!***

#include "Flash.h"
#include "SolidMechanics.h"
#include "constants.h"

#if NDIM==MDIM
#define WRITE_BINARY_GEO
#endif

subroutine sm_uvTecplot(file_name,count)

  use SolidMechanics_Data, only : sm_meshMe,sm_meshComm,sm_NumBodies,sm_BodyInfo
  use sm_assemble_interface, only : sm_assemble_COM
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, only : gr_globalDomain
  use gr_sbData, only : gr_sbBodyInfo
  implicit none

  integer, intent(in) :: count
  character(len=*) :: file_name

  real :: xi,yi,zi,nmlx,nmly,nmlz,pres,fvx,fvy,fvz,wx,wy,wz
  real :: ux,uy,uz,upx,upy,upz, uix, uiy, uiz
  character(len=6) ::index_name
  character(len=6) :: index_ibd
  integer :: i,ie,ibd

  integer :: ierr

  integer :: x_sgn,y_sgn,z_sgn,nwx,nwy,nwz
  real    :: Lx, Ly, Lz
  logical :: vflag
  real, allocatable, dimension(:) :: vx,vy,vz,vxm,vym,vzm

  integer :: plot_nel

  real, parameter :: eps = 1.e-12

#ifdef WRITE_BINARY_GEO
  character(28) :: filename
  integer*4 TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd
  integer*4 visdouble,disdouble
  integer*4 Debug,ZoneType
  character*1 NULLCHR
  integer, allocatable, dimension(:,:) :: icon
#endif

#ifndef WRITE_BINARY_GEO
  !!! Write to Tecplot - Put in IO file:
  write(index_name,"(I6.6)") count
  if (sm_meshMe .eq. MASTER_PE) then
     open(unit=110,file='./IOData/'//file_name//'-'//index_name,form='formatted',&
          status='replace')
     close(110)
  endif
#else

  write(*,*) 'Undefined routine sm_uvTecplot for binary file', sm_meshMe
   
#endif /* WRITE_BINARY_GEO */


  call MPI_BARRIER(sm_meshComm,ierr)

  !! Loop over bodies and write:
  do ibd = 1,sm_NumBodies
     if (sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster)then
 
#ifndef WRITE_BINARY_GEO
        ! Write Particles
        open(unit=110,file='./IOData/'//file_name//'-'//index_name,form='formatted',&
             status='old',position='append')
        write(index_ibd,"(I6.6)") ibd
#if NDIM == MDIM
        write(110,'(A)') 'VARIABLES = "X","Y","Z","ux","uy","uz","uix","uiy","uiz","upx","upy","upz"' 
        !write(110,'(A)') 'VARIABLES = "X","Y","Z","ux","uy","uz","uix","uiy","uiz"' 
        write(110,*)'ZONE I=',gr_sbBodyInfo(ibd)%totalpart,', T="Parts'//index_ibd//'"',', F=POINT' 
        do i = 1,gr_sbBodyInfo(ibd)%totalpart
           xi  =gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i) 
           yi  =gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i)
           zi  =gr_sbBodyInfo(ibd)%particles(POSZ_PART_PROP,i)
           ux  =gr_sbBodyInfo(ibd)%particles(VELX_PART_PROP,i)
           uy  =gr_sbBodyInfo(ibd)%particles(VELY_PART_PROP,i)
           uz  =gr_sbBodyInfo(ibd)%particles(VELZ_PART_PROP,i) 
           uix =gr_sbBodyInfo(ibd)%particles(UITP_PART_PROP,i) 
           uiy =gr_sbBodyInfo(ibd)%particles(VITP_PART_PROP,i) 
           uiz =gr_sbBodyInfo(ibd)%particles(WITP_PART_PROP,i) 
           upx =gr_sbBodyInfo(ibd)%particles(UUPD_PART_PROP,i) 
           upy =gr_sbBodyInfo(ibd)%particles(VUPD_PART_PROP,i) 
           upz =gr_sbBodyInfo(ibd)%particles(WUPD_PART_PROP,i) 
           write(110,'(16g19.10)') xi,yi,zi,ux, uy, uz,uix,uiy,uiz,upx,upy,upz        
        enddo
#else
        ! Write Particles 2D:
        write(110,'(A)') 'VARIABLES = "X", "Y", "ux", "uy", "uix", "uiy", "upx", "upy", "deg"'
        !write(110,'(A)') '# VARIABLES = "X", "Y", "ux", "uy", "uix", "uiy","deg"'
        write(110,*)'ZONE I=',gr_sbBodyInfo(ibd)%totalpart,', T="Parts'//index_ibd//'"',', F=POINT'
        do i = 1,gr_sbBodyInfo(ibd)%totalpart
           xi  =gr_sbBodyInfo(ibd)%particles(POSX_PART_PROP,i)
           yi  =gr_sbBodyInfo(ibd)%particles(POSY_PART_PROP,i)
           ux =gr_sbBodyInfo(ibd)%particles(VELX_PART_PROP,i)
           uy =gr_sbBodyInfo(ibd)%particles(VELY_PART_PROP,i)
           uix =gr_sbBodyInfo(ibd)%particles(UITP_PART_PROP,i) 
           uiy =gr_sbBodyInfo(ibd)%particles(VITP_PART_PROP,i) 
           upx =gr_sbBodyInfo(ibd)%particles(UUPD_PART_PROP,i) 
           upy =gr_sbBodyInfo(ibd)%particles(VUPD_PART_PROP,i) 
           write(110,'(12g19.10)') xi,yi,ux,uy,uix, uiy, upx, upy, 360.0*real(i)/gr_sbBodyInfo(ibd)%totalpart !, upx, upy
        enddo
#endif  /* NDIM */
        close(110)      

#else
        ! Dont write Particle Information in Binary mode.
#endif /* WRITE_BINARY_GEO */     
     endif
  enddo
!    endif
!  enddo

  return

end subroutine sm_uvTecplot

!
!
!
!
!
!
   subroutine sm_CalcDiameters (ibd,firstfileflag)
#include "Flash.h" 
#include "SolidMechanics.h"
#include "constants.h"

 use Simulation_data, Only : Npmax,Npmin, Ypmax 
 use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs, dr_nbegin,       &
                         dr_nend, dr_dt, dr_tmax, dr_simTime, dr_nstep
 use SolidMechanics_data, ONLY :sm_BodyInfo ,sm_meshMe,  &
                        sm_meshComm,sm_structure
 implicit none
     integer, intent(IN) :: ibd,firstfileflag
     real,dimension(NDIM) :: p1,p2,comm
     real :: yi, zi , Daxial, Dtrans
     type(sm_structure), pointer :: body
     integer :: i
     !integer, allocatable :: Npmax(:),Npmin(:)
     integer :: Nnodes


     body => sm_bodyinfo(ibd)
     NNodes=body%nnp
     
     write(*,*) Npmax, Npmin                
     p1=body%qms(Npmax*NDIM-2:Npmax*NDIM,1);
     p2=body%qms(Npmin*NDIM-2:Npmin*NDIM,1);
     
     Daxial=SQRT(sum((p1-p2)*(p1-p2)));
     
     !              Calculating the center of masses
     COMM(:)=0.
     COMM(1)=SUM(body%qms(1:NDIM:NDIM*NNodes-1,1)) /real(nnodes) 
     COMM(2)=SUM(body%qms(2:NDIM:NDIM*Nnodes-1,1)) /real(NNodes)
     COMM(3)=SUM(body%qms(3:NDIM:NDIM*NNodes  ,1)) /real(NNodes)
                
     
     !       Calculate the transverse diameter
     
     yi=body%qms(Ypmax*NDIM-1,1)
     zi=body%qms(Ypmax*NDIM,1)
     Dtrans=2.*sqrt((yi-comm(2)) * (yi-comm(2)) + &
          (zi-comm(3)) * (zi-comm(3)) )                                 
     
     if (firstfileflag.eq.0) then 
        open(unit=14,file="Diamters.dat",status="replace")
     else 
        open(unit=14,file="Diamters.dat",status="old", &
             position="append")
     end if
     write(14,'(I6,1X,F12.8,1X,2F16.12)')dr_nstep,dr_simtime,Daxial,Dtrans
     close(14)
     
   end subroutine sm_CalcDiameters
   

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_Verlet_advance(ibd,restart)
      use SolidMechanics_data, only: sm_structure, sm_BodyInfo,  sm_Numbodies
      use sm_Verlet_data, only: sm_Verlet_type, sm_Verlet_info
      use sm_integinterface, only: sm_verlet_partone, sm_verlet_parttwo
      use sm_assemble_interface, only: sm_assemble_IntForce_rbc
      use Driver_data, Only : dr_nstep
      
      implicit none
      integer, intent(in) :: ibd
      logical, intent(in) :: restart
      real :: time,dt

      ! Internal Variables
      type(sm_structure), pointer :: body
      type(sm_Verlet_type), pointer :: integ    
      integer :: i
      body => sm_BodyInfo(ibd)
      integ=> sm_Verlet_info(ibd)

      

      call Driver_getSimTime(time)
      
      ! Get the current DT
      call Driver_getDT(dt)
      write(*,*) 'dt=',dt,'integ%dt=',integ%dt
      if (dt<integ%dt) integ%dt=dt;
      
      ! Fluid Forces + External Body Forces Fext
      ! This populates body%Hs
      ! Before the integration of the structure the fluid forces are assembled
     
      call sm_assemble_ExtForce(ibd, time)
      
      ! Advance the positions and the velocities of the structure
      call sm_verlet_partone(ibd,integ%dt)
      !write(*,*) 'No. of Bodies: ',sm_NumBodies


      
      
      if (maxval(body%qn(:))>100) then
         write(*,*) 'Check the integration part one '
         stop
      end if
      
     ! Get the forces at the t + lambda*dt
      call sm_assemble_IntForce_rbc(ibd,SM_IOPT_NMIDSTEP)
      
     
      ! advance the position and intermediate velcocities       
      call sm_verlet_parttwo(ibd,integ%dt)
   
  

end subroutine sm_Verlet_advance


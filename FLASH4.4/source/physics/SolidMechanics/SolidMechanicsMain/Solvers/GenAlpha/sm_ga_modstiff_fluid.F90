!     
! File:   Modify the stiffness matrix to inlude linearized terms from the pressure
! Author: tim
!

#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_ga_modstiff_fluid(ibd)
  use SolidMechanics_data,  only: sm_nen,  sm_BodyInfo, sm_structure
  use sm_GenAlpha_data,     only: sm_GenAlpha_type, sm_GenAlpha_info
  use sm_element_interface, only: sm_el10_ga_stiff
  use gr_sbData,            only: gr_sbBodyInfo
  use Driver_interface,     only: Driver_abortFlash
  implicit none

  !IO Variables
  integer, intent(in)    :: ibd ! body number

  ! Define internal variables
  type(sm_structure), pointer :: body
  type(sm_GenAlpha_type), pointer :: integ
  real, dimension(NPART_PROPS) :: particle
  integer :: nee, e, nen_e, p, a1, a2, i1, i2, p1, p2, idx
  real, allocatable, dimension(:,:) :: kfe

  body => sm_BodyInfo(ibd)
  integ => sm_GenAlpha_info(ibd)

  write(*,*) 'Modify GenAlpha stiffness matrix for pressure'

  ! Loop over particles
  do p = 1,gr_sbBodyInfo(ibd)%totalPart
     
     ! get patch info
     particle = gr_sbBodyInfo(ibd)% particles(1:NPART_PROPS,p)
     e        = int( particle(ELEM_PART_PROP) )
     nen_e    = sm_nen( body%ws_eltype(e) )
     nee      = body%max_dofs_per_node*nen_e  ! number of dofs per element, this is based on disp. dofs
     allocate( kfe(nee, nee) )

     select case( body%ws_eltype(e) )

     case( NINE_NODE_QUADRILATERAL )

        call sm_el10_ga_stiff(body, e, particle, integ%beta, integ%dt, kfe)

     case default

        call Driver_abortFlash("Error: unknown element type.")

     end select

     ! Reduce kfe into global stiff matrix
     do a1 = 1,nen_e
        do i1 = 1,body%max_dofs_per_node

           ! local equation number 
           ! (groups dofs together as [u_1 u_2 ..., v_1 v_2 ..., w_1 w_2 ... ]
           !  this is how things are grouped everywhere in 3Dflexible)
           p1 = a1 + nen_e*(i1-1)

           if( body%ID(i1,body%ws_IEN(a1,e)) <= body%neq ) then

              do a2 = 1,nen_e
                 do i2 = 1,body%max_dofs_per_node

                    if( body%ID(i2,body%ws_IEN(a2,e)) <= body%neq ) then

                       p2 = a2 + nen_e*(i2-1)
                       idx = body%ws_LM_cs(p1,p2,e)            
                       ! minus since since this term is coming from an external load
                       body%K(idx) = body%K(idx) - kfe(p1,p2)  
                     
                    end if

                 end do
              end do
           
           end if
           
        end do
     end do

     deallocate(kfe)

  enddo

  return

end subroutine sm_ga_modstiff_fluid

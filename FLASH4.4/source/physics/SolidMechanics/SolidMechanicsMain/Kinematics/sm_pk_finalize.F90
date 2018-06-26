!     
! File:   sm_pk_finalize.F90
! Author: tim
!
! 
subroutine sm_pk_finalize
        use sm_pk_data
        implicit none
        integer :: i
	
	! kill all objects
        if( sm_pk_NumKinematics > 0 ) then
           do i = 1,sm_pk_NumKinematics
              if( allocated( sm_pk_info(i)% params ) ) then
                 deallocate( sm_pk_info(i)%params )
              end if
           end do
        end if

end subroutine sm_pk_finalize

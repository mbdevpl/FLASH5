subroutine rand_gen (zij)
  
  !***********************************
  !     RANDOM NUMBER GENERATOR
  !     HUSSEIN EZZELDIN UMCP FEB 2011
  !************************************
  IMPLICIT none
  
  INTEGER  ::   i_seed, dt_seed(8)
  INTEGER,ALLOCATABLE  :: a_seed(:)
  
  real, intent(INOUT) :: zij
  
  i_seed = 8;
  
  call  RANDOM_SEED(size=i_seed)
  allocate(a_seed(1:i_seed))
  
  CALL RANDOM_SEED(get=a_seed)
  
  CALL DATE_AND_TIME(values=dt_seed)
  a_seed(i_seed)=dt_seed(1);
  a_seed(1)=dt_seed(4)*dt_seed(7)*dt_seed(6)
  
  CALL RANDOM_SEED(put=a_seed)
  DEALLOCATE(a_seed)
  
  CALL random_number(zij)
  
  
END SUBROUTINE RAND_GEN

FUNCTION cross_product(a,b)

  real :: a(3),b(3),cross_product(3)

  cross_product(1)=(a(2)*b(3))   &
                  -(a(3)*b(2))

  cross_product(2)=-((a(1)*b(3)) &
                    -(a(3)*b(1)))

  cross_product(3)=(a(1)*b(2))   &
                  -(a(2)*b(1))

END FUNCTION cross_product


subroutine getTri_Norm(Positions,Elements,normals,centers,areas, volume,& 
     NNodes,NEle)

#include "Flash.h"

  implicit none
  interface
     function cross_product(a,b)   
       real ::a(3),b(3),cross_product(3)        
     end function cross_product
  end interface

  integer,intent(IN) ::  Nele, Nnodes
  integer,intent(IN) :: Elements(3,Nele)
  real,intent(IN)    :: Positions(Nnodes*NDIM)
  real,intent(OUT)   :: normals(Nele*NDIM),centers(Nele*NDIM)
  real,intent(OUT)   :: areas(Nele),volume
  
  ! local variables
  integer :: i
  real  ,dimension(NDIM)::Pos21,Pos31,comm,vol_cen
  integer :: n1,n2,n3
  
  Volume=0.;
  comm(1)=sum(Positions(1:3:NDIM*Nnodes))/real(nnodes) ;
  comm(2)=sum(Positions(2:3:NDIM*Nnodes))/real(nnodes) ;
  comm(3)=sum(Positions(3:3:NDIM*Nnodes))/real(nnodes) ;
 
  do i=1,Nele
     n1=Elements(1,i);
     n2=Elements(2,i);
     n3=Elements(3,i);
 
     ! Calculate all the normals of all the elements
     Pos31 = Positions(NDIM*n3-CONSTANT_TWO:NDIM*n3)- Positions(NDIM*n1-CONSTANT_TWO:NDIM*n1)
     Pos21 = Positions(NDIM*n2-CONSTANT_TWO:NDIM*n2)- Positions(NDIM*n1-CONSTANT_TWO:NDIM*n1)
     
     normals(3*i-2:3*i) = cross_product(Pos21,Pos31)
!!$     normals(3*i-2)= (Pos21(2)*Pos31(3)-Pos21(3)*Pos31(2))
!!$     normals(3*i-1)=-(Pos21(1)*Pos31(3)-Pos21(3)*Pos31(1))
!!$     normals(3*i)  = (Pos21(1)*Pos31(2)-Pos21(2)*Pos31(1))
     
     centers(3*i-2)   =(Positions(NDIM*n1-2)+Positions(NDIM*n2-2)+Positions(NDIM*n3-2))/3.
     centers(3*i-1)   =(Positions(NDIM*n1-1)+Positions(NDIM*n2-1)+Positions(NDIM*n3-1))/3.
     centers(3*i)     =(Positions(NDIM*n1)  +Positions(NDIM*n2)  +Positions(NDIM*n3))  /3.

     ! Calculate the areas and volume of the RBC
     ! area is 1/2 the magnitude of the normal vector
     Areas(i)      = 0.5*sqrt(sum(normals(3*i-2:3*i)*normals(3*i-2:3*i)))
     vol_cen= centers(3*i-2:3*i)-comm;
     !write(*,*) 'volume center= ',vol_cen
     
     Volume= Volume + sum( normals(3*i-2:3*i)* centers(3*i-2:3*i) )
  end do
  
  Volume=Volume/6.;
  write(*,*)'volume=',Volume

  return;
end subroutine getTri_Norm


subroutine  MACHINEEPSILON (eps)
  
  implicit none
  double precision,intent(OUT) :: eps
  double precision  :: MACHEPS=1.D0
  
  MACHEPS = MACHEPS / 2.D0
  !MACHEPS = 1.D0
  do while (1.D0 + MACHEPS / 2.D0 .NE. 1.D0)
     MACHEPS = MACHEPS / 2.D0
     !  IF ( 1.D0 + MACHEPS / 2.D0 .EQ. 1.D0 ) GOTO 110
  end do
  !GO TO 100
  !110  CONTINUE
  eps= MACHEPS;
  print*, MACHEPS
  return
end subroutine MACHINEEPSILON

Subroutine int2char(i,strng)
  implicit none
  integer,intent(IN) :: i
  character(6),intent(OUT) ::  strng
  
  integer :: k, val, valaux
  real ::  val2
  
  valaux=0
  strng = '000000'
  
  
  do k = 6,1,-1
     
     val2 = (i-valaux) / (10**(k-1))
     val = floor(val2)
     
     valaux = valaux + val*(10**(k-1))
     
     !        write(*,*) 7-k,val,valaux
     
     if (val .GE. 1) then
        
        select case (val)
           
        case (1)
           
           strng(7-k:7-k) = "1"
           
        case (2)
           
           strng(7-k:7-k) = '2'
           
           
        case (3)
           
           strng(7-k:7-k) = '3'
           
        case (4)
           
           strng(7-k:7-k) = '4'
           
        case (5)
           
           strng(7-k:7-k) = '5'
           
        case (6)
           
           strng(7-k:7-k) = '6'
           
        case (7)
           
           strng(7-k:7-k) = '7'
           
        case (8)
           
           strng(7-k:7-k) = '8'
           
        case (9)
           
           strng(7-k:7-k) = '9'
           
        end select
        
     endif
     
  enddo
  
end subroutine int2char



  subroutine sortx(Positions,ind,NNodes)
! 
!
!/////////////////////////////////////////////////////////////


#include "constants.h"


      IMPLICIT none
      integer, intent(IN) ::  NNodes
      real, intent(in)    :: Positions(NNodes) 
      real :: temp
      real    :: arr(NNodes)
      integer :: i,j,n,tempind
      integer, intent(out) :: ind(NNodes)
        

      arr(:)= Positions !storing the x positions in array arr
        

      ! Creating the indcies array
      ind(1:NNodes)=(/(i,i=1,NNodes)/)

      !Arranging the array arr in ascending order

      do i=1,NNodes

         do j=1,NNodes-i
              !write(*,*) i, j

            if (arr(j).gt.arr(j+1)) then
               temp=arr(j)
               arr(j)=arr(j+1)
               arr(j+1)=temp
               tempind=ind(j)
               ind(j)=ind(j+1)
               ind(j+1)=tempind
             endif

         enddo

      enddo

      end subroutine sortx


!
!
!
      subroutine sorty(Positions,ind,Nnodes)
!
!
#include "constants.h"

      implicit none
      integer, INTENT(in):: NNodes
      integer, INTENT(out):: ind(Nnodes)
      real, INTENT(in)  :: Positions(NNodes)
      real  :: temp
      real  :: arr(NNodes)
      integer :: i,j,n,tempind
        


      arr(:)=ABS(Positions) !storing the x positions in array arr


      ! Creating the indcies array
      ind(1:NNodes)=(/(i,i=1,NNodes)/)

      !Arranging the array arr in ascending order
      do i=1,NNodes

         do j=1,NNodes-i

            if (arr(j).gt.arr(j+1)) then
               temp=arr(j)
               arr(j)=arr(j+1)
               arr(j+1)=temp
               tempind=ind(j)
               ind(j)=ind(j+1)
               ind(j+1)=tempind

            endif

         enddo

      enddo
      end subroutine sorty

!!****if* source/Grid/GridMain/paramesh/paramesh4/AmrexTransition/gr_sanitizeDataAfterInterp
!!
!! NAME
!!
!!  gr_sanitizeDataAfterInterp
!!
!!
!! SYNOPSIS
!!
!!  call gr_sanitizeDataAfterInterp(integer(in)          :: ntype,
!!                                  character(len=*)(in) :: info,
!!                                  integer(in)          :: layers(MDIM))
!!
!!
!!
!! DESCRIPTION
!!
!!  Given a block node type specification, loop over all of 
!!  the matching blocks and check whether matching solution data in certain
!!  variableslie in a reasonable range of values.
!!
!!  Energies (ENER_VAR and EINT_VAR) are expected to be .ge. gr_smalle,
!!  and the density (DENS_VAR) is expected to be .ge. gr_smallrho,
!!  where gr_smalle and gr_smallrho are lower bounds coming from the
!!  runtime parameters smalle and smlrho, respectively.
!!
!!  For data that do satisfy these expectations, warning messages are
!!  generated, but the offending data is not modified.
!!
!! ARGUMENTS
!!
!!   ntype - the node type of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS)
!!
!!   layers - number of guardcell layers to be included in the check 
!!
!! NOTES
!!
!!  The checks are based on gr_conserveToPrimitive, which is called to
!!  convert solution data back from conserved form when using the old
!!  (convertToConsvdForMeshCalls) way of ensuring that the mesh
!!  handles data interpolation in conserved form.
!!
!!  This is meant to be called where gr_conserveToPrimitive used to be
!!  called when using the new (convertToConsvdInMeshInterp) way of
!!  ensuring that the mesh handles data interpolation in conserved
!!  form.
!!
!! SEE ALSO
!!
!!  gr_conserveToPrimitive
!!
!!***

!!REORDER(4): solnData

#define DEBUG_CONSCONV

subroutine gr_sanitizeDataAfterInterp(ntype, info, layers)

  use Grid_interface, ONLY : Grid_getBlkPtr
  use Grid_interface, ONLY : Grid_copyF4DataToMultiFabs
  use Grid_data, ONLY : gr_smallrho,gr_smalle, gr_meshMe
  use gr_specificData, ONLY : gr_sanitizeDataMode, gr_sanitizeVerbosity
  use Logfile_interface, ONLY : Logfile_stamp
  use physicaldata, ONLY: gcell_on_cc
  use paramesh_dimensions, ONLY: il_bnd,iu_bnd,jl_bnd,ju_bnd,kl_bnd,ku_bnd, kl_bndi, ndim
  use gr_amrextData, ONLY : gr_amrextUnkMFs
  use block_iterator !, ONLY : block_iterator_t
  use block_metadata, ONLY : block_metadata_t

  implicit none
#include "constants.h"
#undef REAL_FORMAT
#define REAL_FORMAT "(1PG23.16)"
#include "Flash.h"

  integer, intent(IN) :: ntype
  character(len=*), intent(IN) :: info
  integer,dimension(MDIM), intent(IN):: layers

  real, dimension(:,:,:,:), pointer :: solnData
  integer,dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: n, blockID
  integer :: myPe, i,j
  integer :: iskip, jskip, kskip
  integer :: il,iu,jl,ju,kl,ku
  integer :: count, level
  integer :: kwrite,locs(3),kReorder(1:ku_bnd-kl_bnd+1),nReorder
  character(len=32), dimension(4,2) :: block_buff
  character(len=32)                 :: number_to_str
  type(block_iterator_t) :: itor
  type(block_metadata_t) :: blockDesc

111 format (a,a,a1,(1x,a18,'=',a),(1x,a2,'=',a5),(1x,a5,'=',a),(1x,a4,'=',a))
112 format (i3,1x,24(1x,1G8.2))
113 format (' :,',i2,',',i2,1x,24(1x,1G8.2))

  if (gr_sanitizeDataMode == 0) return ! IMMEDIATE RETURN

  call Grid_copyF4DataToMultiFabs(CENTER, gr_amrextUnkMFs, nodetype=ACTIVE_BLKS)

  iskip = NGUARD - layers(IAXIS)
  jskip = (NGUARD - layers(JAXIS)) * K2D
  kskip = (NGUARD - layers(KAXIS)) * K3D

  nReorder = 0
  count = 0
  
  ! DEVNOTE: Is it *always* correct to use CENTER here?
  itor = block_iterator_t(ntype)
  do while (itor%is_valid())
     count = count + 1
     call itor%blkMetaData(blockDesc)
     blockID = blockDesc%id
     if (blockID < 0) blockID = -count
     level = blockDesc%level
     call Grid_getBlkPtr(blockDesc, solnData, localFlag=.TRUE.)
     blkLimits   = blockDesc%localLimits
     blkLimitsGC = blockDesc%localLimitsGC

     il = blkLimitsGC(LOW,IAXIS) + iskip
     iu = blkLimitsGC(HIGH,IAXIS) - iskip
     jl = blkLimitsGC(LOW,JAXIS) + jskip
     ju = blkLimitsGC(HIGH,JAXIS) - jskip
     kl = blkLimitsGC(LOW,KAXIS) + kskip
     ku = blkLimitsGC(HIGH,KAXIS) - kskip


#ifdef DENS_VAR
     if (gcell_on_cc(DENS_VAR) .OR. gr_sanitizeDataMode == 2) then
        ! small limits -- in case the interpolants are not monotonic
        if (any(solnData(DENS_VAR,il:iu,jl:ju,kl:ku) .LT. gr_smallrho)) then
           kwrite = kl_bndi
           if (gr_sanitizeVerbosity .GE. 5 .AND. ndim==3) then
              call set_kReorder
!              print*,'kReorder(1:nReorder)',kReorder(1:nReorder)
              locs = minloc(solnData(DENS_VAR,il:iu,jl:ju,kReorder(1:nReorder)))
!              print*,'LOCS:',locs
              kwrite = kReorder(locs(3))
           end if
           write (block_buff(1,1), '(a18)') 'min. unk(DENS_VAR)'
           !        write (number_to_str, '('//REAL_FORMAT//',a1)') minval(solnData(DENS_VAR,il:iu,jl:ju,kl:ku)), ','
           write (number_to_str, '(G30.22)') minval(solnData(DENS_VAR,il:iu,jl:ju,kl:ku))
           write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(2,1), '(a)') 'PE'
           write (number_to_str, '(i7)') gr_meshMe
           write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(3,1), '(a)') 'block'
           write (number_to_str, '(i7)') blockID
           write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(4,1), '(a)') 'level'
           write (number_to_str, '(i7)') level
           write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

           if (gr_sanitizeVerbosity .GE. 1) call Logfile_stamp( block_buff, 4, 2, 'WARNING '//info)
           if (gr_sanitizeVerbosity .GE. 4) print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,4)

           if (gr_sanitizeVerbosity .GE. 5) then
              do j=ju_bnd,jl_bnd,-1
                 if (kwrite==kl_bndi) then
                 ! For 3D, this prints a slice at the lowest k index that is interior - KW
                    print 112, j, (solnData(DENS_VAR,i,j,kwrite), i=il_bnd,iu_bnd)
                 else
                    print 113, j,kwrite, (solnData(DENS_VAR,i,j,kwrite), i=il_bnd,iu_bnd)
                 end if
              end do
           end if

           if (gr_sanitizeDataMode == 3) then
              solnData(DENS_VAR,il:iu,jl:ju,kl:ku) = max(gr_smallrho,solnData(DENS_VAR,il:iu,jl:ju,kl:ku))
           end if
           if (gr_sanitizeDataMode == 4) call Driver_abortFlash("DENS var below acceptable minimum")
        end if
     end if
#endif

#ifdef ENER_VAR               
     if (gcell_on_cc(ENER_VAR) .OR. gr_sanitizeDataMode == 2) then
        ! energy
        if (any(solnData(ENER_VAR,il:iu,jl:ju,kl:ku) .LT. gr_smalle*0.999999999)) then
           kwrite = kl_bndi
           if (gr_sanitizeVerbosity .GE. 5 .AND. ndim==3) then
              call set_kReorder
              locs = minloc(solnData(ENER_VAR,il:iu,jl:ju,kReorder(1:nReorder)))
              kwrite = kReorder(locs(3))
           end if
           write (block_buff(1,1), '(a)') 'min. unk(ENER_VAR)'
           write (number_to_str, '('//REAL_FORMAT//')') minval(solnData(ENER_VAR,il:iu,jl:ju,kl:ku))
           write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(2,1), '(a)') 'PE'
           write (number_to_str, '(i7)') gr_meshMe
           write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(3,1), '(a)') 'block'
           write (number_to_str, '(i7)') blockID
           write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(4,1), '(a)') 'level'
           write (number_to_str, '(i7)') level
           write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

           if (gr_sanitizeVerbosity .GE. 1) call Logfile_stamp( block_buff, 4, 2, 'WARNING '//info)
           if (gr_sanitizeVerbosity .GE. 4) print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,4)

           if (gr_sanitizeVerbosity .GE. 5) then
              do j=ju_bnd,jl_bnd,-1
                 if (kwrite==kl_bndi) then
                    print 112, j, (solnData(ENER_VAR,i,j,kwrite), i=il_bnd,iu_bnd) 
                 else
                    print 113, j,kwrite, (solnData(ENER_VAR,i,j,kwrite), i=il_bnd,iu_bnd)
                 end if
              end do
           end if

           if (gr_sanitizeDataMode == 3) then
              solnData(ENER_VAR,il:iu,jl:ju,kl:ku) = max(gr_smalle,solnData(ENER_VAR,il:iu,jl:ju,kl:ku))
           end if
           if (gr_sanitizeDataMode == 4) call Driver_abortFlash("ENER var below acceptable minimum")
        end if
     end if
#endif
#ifdef EINT_VAR
     if (gcell_on_cc(EINT_VAR) .OR. gr_sanitizeDataMode == 2) then
        if (any(solnData(EINT_VAR,il:iu,jl:ju,kl:ku) .LT. gr_smalle*0.999999999)) then
           kwrite = kl_bndi
           if (gr_sanitizeVerbosity .GE. 5 .AND. ndim==3) then
              call set_kReorder
              locs = minloc(solnData(EINT_VAR,il:iu,jl:ju,kReorder(1:nReorder)))
              kwrite = kReorder(locs(3))
           end if
           write (block_buff(1,1), '(a)') 'min. unk(EINT_VAR)'
           write (number_to_str, '('//REAL_FORMAT//')') minval(solnData(EINT_VAR,il:iu,jl:ju,kl:ku))
           write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(2,1), '(a)') 'PE'
           write (number_to_str, '(i7)') gr_meshMe
           write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(3,1), '(a)') 'block'
           write (number_to_str, '(i7)') blockID
           write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(4,1), '(a)') 'level'
           write (number_to_str, '(i7)') level
           write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

           if (gr_sanitizeVerbosity .GE. 1) call Logfile_stamp( block_buff, 4, 2, 'WARNING '//info)
           if (gr_sanitizeVerbosity .GE. 4) print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,4)

           if (gr_sanitizeVerbosity .GE. 5) then
              do j=ju_bnd,jl_bnd,-1
                 if (kwrite==kl_bndi) then
                    print 112, j, (solnData(EINT_VAR,i,j,kwrite), i=il_bnd,iu_bnd) 
                 else
                    print 113, j,kwrite, (solnData(EINT_VAR,i,j,kwrite), i=il_bnd,iu_bnd)
                 end if
              end do
           end if

           if (gr_sanitizeDataMode == 3) then
              solnData(EINT_VAR,il:iu,jl:ju,kl:ku) = max(gr_smalle,solnData(EINT_VAR,il:iu,jl:ju,kl:ku))
           end if
           if (gr_sanitizeDataMode == 4) call Driver_abortFlash("EINT var below acceptable minimum")
        end if
     end if
#endif

     call itor%next()
  end do

  call itor%destroy_iterator()

  return 

contains
  integer function nodetype(i)
    integer,intent(IN) :: i
    nodetype = 999              !fake for output
  end function nodetype
  subroutine set_kReorder
    integer :: i,j,k
    if (nReorder==0) then       !We need to do this only once per gr_sanitize... call.

       !  The following code sets
       !  kReorder = (/5,6,7,8,9,10,11,12,4,13,3,14,2,15,1,16/) ! for kl:ku = 1:16

       i = 0
       do k = kl, ku
          if (k>NGUARD .AND. k .LE. ku_bnd-NGUARD) then
             i = i+1
             kReorder(i) = k
          end if
       end do
       do j = 1, layers(KAXIS)
          k = NGUARD+1-j
          if (k .LE. NGUARD) then
             i = i+1
             kReorder(i) = k
          end if
          k = ku_bnd-NGUARD+j
          if (k > ku_bnd-NGUARD) then
             i = i+1
             kReorder(i) = k
          end if
       end do
       nReorder = i
    end if
  end subroutine set_kReorder

end subroutine gr_sanitizeDataAfterInterp

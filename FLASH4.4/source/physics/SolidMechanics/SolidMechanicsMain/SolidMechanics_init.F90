!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SolidMechanics_init
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!
!! VARIABLES
!!
!!
!! DESCRIPTION
!!
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine SolidMechanics_init(restart)
  use HDF5
  use SolidMechanics_data, only : sm_MeshMe,sm_NumProcs, sm_NumBodies, sm_BodyInfo, sm_meshComm, &
                                  sm_metaBodies, sm_metaBodInfo, sm_GlobalGroup, &
                                  sm_gravX, sm_gravY, sm_gravZ
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use sm_iointerface, only : sm_ioReadSolid,sm_ioInit
  use sm_pk_interface, only: sm_pk_init
  use sm_integinterface, only: sm_PredCorr_init, sm_GenAlpha_init, sm_Verlet_init
  use Driver_interface, only : Driver_getMype,Driver_getNumProcs, &
                               Driver_abortFlash, Driver_getSimTime, &
                               Driver_getComm
  use sm_contact_interface, only : sm_contact_init

  implicit none
#include "constants.h"
#include "SolidMechanics.h"
#include "sm_integrator.h"
#include "Flash_mpi.h"
#include "Flash.h"

  logical, INTENT(IN) :: restart
  integer             :: h5ferr, ibd, load_pk_flag

  integer, allocatable, dimension(:) :: sm_read_PE

  integer           :: imet, imetpr, nummetabods
  character(LEN=20) :: metaBod
  integer           :: metpr_start, metpr_end

  integer :: deltaProcs,i,iproc,ierr

  ! Read Number of Bodies from runtime file:
  call RuntimeParameters_get("sm_NumBodies", sm_NumBodies)
  allocate(sm_BodyInfo(sm_NumBodies))

  ! Which Processor am I:
  call Driver_getMype(GLOBAL_COMM, sm_MeshMe)

  ! Number of Processors:
  call Driver_getNumProcs(GLOBAL_COMM, sm_NumProcs)

  ! Get Comm
  call Driver_getComm( GLOBAL_COMM, sm_meshComm )

  ! Read MetaBody info:
  call RuntimeParameters_get("sm_metaBodies", sm_metaBodies)

  ! Read Gravity components:
  call RuntimeParameters_get("gravX", sm_gravX)
  call RuntimeParameters_get("gravY", sm_gravY)
#if NDIM == MDIM
  call RuntimeParameters_get("gravZ", sm_gravZ)
#else
  sm_gravZ = 0.0
#endif

  ! Initialize Metabody, = 0 means this ibd is not part of a set.
  do ibd = 1,sm_NumBodies
     sm_BodyInfo(ibd)%MetaBody = CONSTANT_ZERO
  enddo

  ! Open HDF5 fortran library
  call h5open_f(h5ferr)

  ! Disable error printing. Errors are stored in stack when
  ! h5lexists_f is used with /field/variable, and field does not exist - MV
  CALL h5eset_auto_f(0, h5ferr)

  !********
  ! Read in the bodies
  !
  if (sm_metabodies .gt. 0) then ! With metabodies, same structural object distributed
     ! on sevelar bodies that contain parts of the wet surface.

     allocate(sm_metaBodInfo(sm_metabodies))

     ! Loop to load how many processors per metaBody
     nummetabods = 0

     do imet = 1,sm_metabodies

        write(metaBod,"(A,I2.2)") 'sm_metaBodyProcs_',imet

        call RuntimeParameters_get(trim(metaBod),sm_metaBodInfo(imet)%NumProcs)

        ! Body and processor count run together for now, that is one bod -> one new proc.
        sm_metaBodInfo(imet)%NumBods = sm_metaBodInfo(imet)%NumProcs

        nummetabods = nummetabods + sm_metaBodInfo(imet)%NumProcs

     enddo


     ! Test total sm_BodyInfo entries is sum(sm_metaBody_Procs(:))
     if (nummetabods .ne. sm_NumBodies) &
          call Driver_abortFlash('SolidMechanics Init: Metabody Procs sum not equal to sm_NumBodies.')

     if (nummetabods .gt. sm_NumProcs) &
          call Driver_abortFlash('SolidMechanics Init: Metabody Procs sum greater than sm_NumBodies.')

     ! Assign MetaBods
     nummetabods = 0
     do imet = 1,sm_metabodies

        ! Only processors member of the metabody allocate mbcomm_Members:
        ! Allocation takes place on the body iproc is the master:
        metpr_start = nummetabods
        !metpr_end   = nummetabods + sm_metaBodInfo(imet)%NumProcs - 1

        ! Processors that will handle this metabod imet:
        allocate(sm_metaBodInfo(imet)%Members(0:sm_metaBodInfo(imet)%NumProcs-1))
        do i=0,sm_metaBodInfo(imet)%NumProcs-1
           sm_metaBodInfo(imet)%Members(i) = metpr_start + i
        enddo

        nummetabods = nummetabods + sm_metaBodInfo(imet)%NumProcs

     enddo

     ! Build the groups and get communicators:
     call MPI_Comm_group(sm_meshComm,sm_GlobalGroup,ierr)
     ibd = 0
     iproc =-1
     do imet = 1,sm_metabodies

        ! Group of processors managing metabody is created in sm_metaBody_Group
        call MPI_Group_incl(sm_GlobalGroup, sm_metaBodInfo(imet)%NumProcs, &
             sm_metaBodInfo(imet)%Members,               &
             sm_metaBodInfo(imet)%Group, ierr)

        sm_metaBodInfo(imet)%comm = MPI_UNDEFINED
        ! Get communicators for each individual metabody
        call MPI_Comm_create(sm_meshComm, sm_metaBodInfo(imet)%Group, &
             sm_metaBodInfo(imet)%comm, ierr);

        ! Now fill out each body component mbcomm, mbcomm_Me, mbcomm_NumProcs
        do imetpr = 1,sm_metaBodInfo(imet)%NumBods

           ! Body and processor count run together for now, that is one bod -> one new proc.
           ibd   = ibd + 1
           iproc = iproc + 1

           ! Metabody this Body belongs to:
           sm_BodyInfo(ibd)%MetaBody = imet

           ! Master for this body:
           sm_BodyInfo(ibd)%BodyMaster = iproc

           ! Metabody Group:
           sm_BodyInfo(ibd)%mbcomm_Group = sm_metaBodInfo(imet)%Group

           ! Metabody communicator
           sm_BodyInfo(ibd)%mbcomm = sm_metaBodInfo(imet)%comm

           ! Number of processors that share Metabody ibd is part of:
           sm_BodyInfo(ibd)%mbcomm_NumProcs = sm_metaBodInfo(imet)%NumProcs

           ! My rank number on the local group: ? Might not be needed here
           if (sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then

              call MPI_Group_rank (sm_BodyInfo(ibd)%mbcomm_Group, sm_BodyInfo(ibd)%mbcomm_Me,ierr)

           else

              sm_BodyInfo(ibd)%mbcomm_Me = MPI_UNDEFINED

           endif


        enddo

     enddo

     do ibd=1,sm_NumBodies

        !Read Each Body
        if (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then
           call sm_ioReadSolid(ibd)
        endif

     enddo

     ! Test values:
     !     do imet = 1,sm_metabodies
     !       write(*,*) ' '
     !       write(*,*) imet,', Me=',sm_meshMe,', Group',sm_metaBodInfo(imet)%Group,', comm=', sm_metaBodInfo(imet)%comm,MPI_COMM_NULL
     !       !', Members=',sm_metaBodInfo(imet)%Members

     !       call mpi_barrier(sm_meshComm, ierr)
     !     enddo

     !     do ibd=1,sm_NumBodies

     !       if(sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then
     !         write(*,*) ' '
     !         write(*,*) 'Body=',ibd,', Master=',sm_BodyInfo(ibd)%BodyMaster,', MetaBody=',sm_BodyInfo(ibd)%MetaBody
     !         write(*,*) 'MetaProcs=',sm_BodyInfo(ibd)%mbcomm_NumProcs,', mbcomm_Me=',sm_BodyInfo(ibd)%mbcomm_Me
     !         write(*,*) 'Global Comm=',sm_meshComm,', mbcomm=',sm_BodyInfo(ibd)%mbcomm
     !         write(*,*) 'Global Group=',sm_GlobalGroup,', mbcomm_Group=',sm_BodyInfo(ibd)%mbcomm_Group
     !       endif

     !       call mpi_barrier(sm_meshComm, ierr)
     !     enddo
     !     call Driver_abortFlash('SolidMechanics Init: up to End of Metabod Init')

  else ! Case no metabodies, all separate bodies

     ! Which Processor will read which body?
     allocate(sm_read_PE(sm_NumBodies)); sm_read_PE(:) = -99999

!!! HERE code to define processor number to read each body in sm_read_PE(1:sm_NumBodies)
     if (sm_NumBodies .ge. sm_NumProcs) then
        deltaProcs = sm_NumBodies/sm_NumProcs
        do ibd=1,sm_NumProcs
           sm_read_PE((ibd-1)*deltaProcs+1:ibd*deltaProcs) = ibd-1
        enddo
        ! Read remainder of Bodies:
        ibd = 1
        do iproc = deltaprocs*sm_NumProcs+1,sm_NumBodies
           sm_read_PE(iproc) = ibd-1
           ibd = ibd + 1
        enddo
     else
        do ibd=1,sm_NumBodies
           sm_read_PE(ibd) = ibd-1
        enddo
     endif

     ! Now read Bodies:
     if (sm_MeshMe .eq. MASTER_PE) write(*,*) 'Reading Bodies ...Num of bodies',sm_Numbodies
     do ibd=1,sm_NumBodies

        !Read Each Body
        sm_BodyInfo(ibd)%BodyMaster = sm_read_PE(ibd)
        if (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then
           call sm_ioReadSolid(ibd)  !PROC BodyMaster = sm_read_PE(ibd) READS ALL BODIES.
        endif

     enddo

     deallocate(sm_read_PE)

  endif ! Case metabodies

  call mpi_barrier(MPI_COMM_WORLD, ierr)
  if (sm_MeshMe .eq. MASTER_PE) write(*,*) 'Done.'

  ! Load Prescribed Kinematics Library - all procs read.
  load_pk_flag = SM_TRUE
  ! Initalize the prescribed kinematics
  if( load_pk_flag == SM_TRUE ) then
     call sm_pk_init()
  end if

  !
  ! Init the solvers
  !
  do ibd=1,sm_NumBodies
     if (sm_MeshMe .eq. sm_BodyInfo(ibd)%BodyMaster) then


        select case(sm_BodyInfo(ibd)%IntegMethod)

           ! Generalized Alpha method
        case(SOLIDINTEG_GENALPHA)
           call sm_GenAlpha_init(ibd, restart)

           ! Predictor-Corrector method
        case(SOLIDINTEG_PREDCORR)
           call sm_PredCorr_init(ibd,restart)

           ! VERLET
        case(SOLIDINTEG_MODVVERLET)
           call sm_Verlet_init(ibd,restart)

        case default
           call Driver_abortFlash('SolidMechanics: Unknown advance integration scheme.')

        end select

        ! 2D test
        !sm_BodyInfo(ibd)%qn(3)  = 1. ! theta0 = 1
        !sm_BodyInfo(ibd)%qdn(1) = 1. ! theta0 = 1
        !sm_BodyInfo(ibd)%qdn(2) = 1. ! theta0 = 1

     end if
  enddo


  ! Call sm_contact_init if no. of numbers is more than 1
  if (sm_NumBodies>1) call sm_contact_init()


  ! Call io Init routine:
  call sm_ioInit(restart)

  ! Enable error printing. Errors are stored in stack when
  ! h5lexists_f is used with /field/variable, and field does not exist - MV
  CALL h5eset_auto_f(1, h5ferr)

  !call Driver_abortFlash('SolidMechanics Init: up to End of Init')

  return

end subroutine SolidMechanics_init

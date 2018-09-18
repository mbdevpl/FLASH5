!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/GPU/bn_xnetInit
!!
!! NAME
!!
!!  bn_xnetInit
!!
!!
!! SYNOPSIS
!!
!!  bn_xnetInit ( character(*), intent(IN) :: data_dir,
!!                character(80), intent(OUT) :: data_desc )
!!
!! DESCRIPTION
!!
!!  Calls XNet routines to read and broadcast XNet data 
!!  Called from bn_initNetwork.
!!
!! ARGUMENTS
!!
!!   data_dir -- nuclear data directory containing REACLIB-formatted data
!!   data_desc -- brief description of  network
!!
!!***

subroutine bn_xnetInit(data_dir,data_desc)
  use Driver_data, ONLY : dr_globalComm
  use nuclear_data, ONLY : ny, read_nuclear_data
  use reaction_data, ONLY : read_reaction_data
  use xnet_abundances, ONLY : ystart, yo, y, yt, ydot
  use xnet_conditions, ONLY : t, tt, to, tdel, tdel_next, tdel_old, t9t, rhot, yet, &
    t9, rho, ye, t9o, rhoo, yeo, t9dot, cv, etae, detaedt9, nt, ntt, nto, ints, intso, &
    nstart, tstart, tstop, tdelstart, t9start, rhostart, yestart, nh, th, t9h, rhoh, yeh, nhmx
  use xnet_controls, ONLY : idiag, iheat, iscrn, iprocess, iweak, kmon, ktot, nnucout, &
    nnucout_string, inucout, output_nuc, nzbatchmx, lzactive, myid, nproc, mythread, nthread, &
    lun_diag, lun_ev, lun_ts, lun_stdout
  use xnet_eos, ONLY : eos_initialize
  use xnet_flux, ONLY : flux_init
  use xnet_gpu, ONLY : gpu_init
  use xnet_integrate_bdf, ONLY : bdf_init
  use xnet_jacobian, ONLY : read_jacobian_data
  use xnet_match, ONLY : read_match_data
  use xnet_parallel, ONLY : parallel_initialize, parallel_myproc, parallel_nprocs, parallel_IOProcessor
  use xnet_preprocess, ONLY : net_preprocess
  use xnet_screening, ONLY : screening_init
  use xnet_util, ONLY : name_ordered

  !$ use omp_lib

  implicit none

  ! Input variables
  character(*), intent(in) :: data_dir

  ! Output variables
  character(80), intent(out) :: data_desc

  ! Local variables
  character(80) :: diag_file
  character(80), parameter :: diag_file_base = 'xnet_diag'

  Call parallel_initialize(dr_globalComm)
  myid = parallel_myproc()
  nproc = parallel_nprocs()

  Call gpu_init

  ! Initialize MPI/OpenMP identifiers
  mythread = 1
  nthread = 1
  !$omp parallel default(shared)
  !$ mythread = omp_get_thread_num()
  !$omp master
  !$ nthread = omp_get_num_threads()
  !$omp end master
  !$omp end parallel

  !$omp parallel default(shared) private(diag_file)

  ! Open diagnositic output file
  if ( idiag >= 0 ) then
    diag_file = trim(diag_file_base)
    call name_ordered(diag_file,myid,nproc)
    call name_ordered(diag_file,mythread,nthread)
    open(newunit=lun_diag, file=diag_file)
  else
    lun_diag = lun_stdout
  endif

  ! Allocate control arrays
  allocate (lzactive(nzbatchmx))
  allocate (iweak(nzbatchmx),lun_ev(nzbatchmx),lun_ts(nzbatchmx))
  allocate (kmon(2,nzbatchmx),ktot(5,nzbatchmx))
  !$omp end parallel

  allocate (output_nuc(nnucout),inucout(nnucout))
  write(nnucout_string,"(i4)") nnucout
  nnucout_string = adjustl(nnucout_string)

  if ( iprocess > 0 .and. parallel_IOProcessor() ) call net_preprocess( lun_stdout, data_dir, data_dir )

  ! Read and distribute nuclear and reaction data
  call read_nuclear_data(data_dir,data_desc)
  call read_reaction_data(data_dir)

  ! Read data on matching forward and reverse reactions
  call read_match_data(data_dir)

  ! Initialize screening
  Call screening_init

  ! Initialize flux tracking
  call flux_init

  ! Read and distribute jacobian matrix data
  call read_jacobian_data(data_dir)

  ! Initialize EoS for screening or self-heating
  if ( iscrn > 0 .or. iheat > 0 ) call eos_initialize

  ! Initialize BDF integrator
  call bdf_init

  !$omp parallel default(shared)

  ! Allocate abundance arrays
  allocate (y(ny,nzbatchmx),yo(ny,nzbatchmx),yt(ny,nzbatchmx),ydot(ny,nzbatchmx),ystart(ny,nzbatchmx))

  ! Allocate conditions arrays
  allocate (t(nzbatchmx),tt(nzbatchmx),to(nzbatchmx), &
    &       tdel(nzbatchmx),tdel_next(nzbatchmx),tdel_old(nzbatchmx), &
    &       t9(nzbatchmx),t9t(nzbatchmx),t9o(nzbatchmx), &
    &       rho(nzbatchmx),rhot(nzbatchmx),rhoo(nzbatchmx), &
    &       ye(nzbatchmx),yet(nzbatchmx),yeo(nzbatchmx), &
    &       nt(nzbatchmx),ntt(nzbatchmx),nto(nzbatchmx), &
    &       ints(nzbatchmx),intso(nzbatchmx), &
    &       t9dot(nzbatchmx),cv(nzbatchmx),etae(nzbatchmx),detaedt9(nzbatchmx))

  ! Allocate thermo history arrays
  allocate (nh(nzbatchmx),nstart(nzbatchmx), &
    &       tstart(nzbatchmx),tstop(nzbatchmx),tdelstart(nzbatchmx), &
    &       t9start(nzbatchmx),rhostart(nzbatchmx),yestart(nzbatchmx), &
    &       th(nhmx,nzbatchmx),t9h(nhmx,nzbatchmx),rhoh(nhmx,nzbatchmx),yeh(nhmx,nzbatchmx))

  !$omp end parallel

  return
end subroutine bn_xnetInit
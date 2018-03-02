!!****if* source/Grid/GridMain/AMR/Amrex/gr_specificData
!!
!! NAME
!!  gr_specificData
!!
!! SYNOPSIS
!!
!!  use gr_specificData
!!
!! DESCRIPTION
!!
!!  Module that contains some private data specific to a Grid implementation.
!!
!!  This includes data that is used by one or more specific implementations,
!!  is probably not of general enough use to be in the Grid_data module,
!!  but is also not part of the internal implementation of an underlying
!!  Grid kernel.
!!
!!  The current Amrex version is derived from the paramesh version, and
!!  should probably change more.
!!***

!!REORDER(5):scratch, scratch_ctr, scratch_facevar[xyz]

Module gr_specificData

  implicit none

#include "constants.h"
#include "Flash.h"

!!$  integer, save :: gr_iloGc = GRID_ILO_GC
!!$  integer, save :: gr_ihiGc = GRID_IHI_GC
!!$  integer, save :: gr_jloGc = GRID_JLO_GC
!!$  integer, save :: gr_jhiGc = GRID_JHI_GC
!!$  integer, save :: gr_kloGc = GRID_KLO_GC
!!$  integer, save :: gr_khiGc = GRID_KHI_GC
!!$
!!$  integer, save :: gr_ilo = GRID_ILO
!!$  integer, save :: gr_ihi = GRID_IHI
!!$  integer, save :: gr_jlo = GRID_JLO
!!$  integer, save :: gr_jhi = GRID_JHI
!!$  integer, save :: gr_klo = GRID_KLO
!!$  integer, save :: gr_khi = GRID_KHI

!! Define the block information

!  Type define
!  Not currently used with Amrex
!!$  type gridBlock
!!$     !!cornerID is integer coordinates of the lower left cornor
!!$     !! (ie the smallest point) of a block
!!$     integer,dimension(MDIM) :: cornerID
!!$     !! atmost 2 neighbors, 2faces along
!!$     !! each dimension, hence.
!!$     real,dimension(3,GRID_IHI_GC) :: firstAxisCoords
!!$     real,dimension(3,GRID_JHI_GC) :: secondAxisCoords
!!$     real,dimension(3,GRID_KHI_GC) :: thirdAxisCoords
!!$     integer :: blockType
!!$  end type gridBlock
!!$
!!$  type(gridBlock),save,dimension(MAXBLOCKS),target :: gr_oneBlock


!  Scratchy stuff not currently used with Amrex either
#ifdef BSS_GRID_ARRAYS
#if NSCRATCH_GRID_VARS > 0
  real,target,dimension(SCRATCH_GRID_VARS_BEGIN:SCRATCH_GRID_VARS_END,&
                        GRID_IHI_GC+1,&
                        GRID_JHI_GC+1,&
                        GRID_KHI_GC+1,&
                        MAXBLOCKS) :: scratch
#else
  real,target,dimension(1,1,1,1,1) :: scratch
#endif

#if NSCRATCH_CENTER_VARS > 0
  real,target,dimension(SCRATCH_CENTER_VARS_BEGIN:SCRATCH_CENTER_VARS_END,&
                        GRID_IHI_GC,&
                        GRID_JHI_GC,&
                        GRID_KHI_GC,&
                        MAXBLOCKS) :: scratch_ctr
#else
  real,target,dimension(1,1,1,1,1):: scratch_ctr
#endif

#if NSCRATCH_FACEX_VARS > 0
  real,target,dimension(SCRATCH_FACEX_VARS_BEGIN:SCRATCH_FACEX_VARS_END,&
                        GRID_IHI_GC+1,&
                        GRID_JHI_GC,  &
                        GRID_KHI_GC,  &
                        MAXBLOCKS) :: scratch_facevarx
#else
  real, target,dimension(1,1,1,1,1):: scratch_facevarx
#endif

#if NSCRATCH_FACEY_VARS > 0
  real,target,dimension(SCRATCH_FACEY_VARS_BEGIN:SCRATCH_FACEY_VARS_END,&
                        GRID_IHI_GC,  &
                        GRID_JHI_GC+1,&
                        GRID_KHI_GC,  &
                        MAXBLOCKS) :: scratch_facevary
#else
  real, target,dimension(1,1,1,1,1):: scratch_facevary
#endif

#if NSCRATCH_FACEZ_VARS > 0
  real,target,dimension(SCRATCH_FACEZ_VARS_BEGIN:SCRATCH_FACEZ_VARS_END,&
                        GRID_IHI_GC,  &
                        GRID_JHI_GC,  &
                        GRID_KHI_GC+1,&
                        MAXBLOCKS) :: scratch_facevarz
#else
  real, target,dimension(1,1,1,1,1):: scratch_facevarz
#endif

#else
  real, save, target, allocatable :: scratch(:,:,:,:,:)
  real, save, target, allocatable :: scratch_ctr(:,:,:,:,:)
  real, save, target, allocatable :: scratch_facevarx(:,:,:,:,:)
  real, save, target, allocatable :: scratch_facevary(:,:,:,:,:)
  real, save, target, allocatable :: scratch_facevarz(:,:,:,:,:)
#endif

  integer, save, allocatable :: gr_nToLeft(:)

  ! variables for making visible some Grid metainformation in a PARAMESH-like form
  ! to a friendly IO implementation
  integer, save              :: gr_ioLocalNumBlocks
  integer, save, allocatable :: gr_ioBlkLrefine(:)
  integer, save, allocatable :: gr_ioBlkNodeType(:)
  real,    save, allocatable :: gr_ioBlkCoords(:,:)
  real,    save, allocatable :: gr_ioBlkBsize(:,:)
  real,    save, allocatable :: gr_ioBlkBoundBox(:,:,:)

  logical, save :: gr_enableMaskedGCFill
  integer, save :: gr_sanitizeDataMode, gr_sanitizeVerbosity
  integer,save :: gr_globalOffset !stores beginning blk offset for a proc

end Module gr_specificData

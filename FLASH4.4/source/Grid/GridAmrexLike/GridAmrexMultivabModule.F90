module famrex_multivab_module

  use,intrinsic :: iso_c_binding
  use famrex_box_module
#include "Flash.h"
#ifdef FLASH_GRID_PARAMESH
  use tree,ONLY: lnblocks,lrefine,lrefine_max
#else
  use Grid_data, ONLY : lnblocks
#endif
  use Grid_interface, ONLY: Grid_getBlkPtr
  use Grid_interface, ONLY: Grid_getBlkIndexLimits
  use Grid_interface, ONLY: Grid_getBlkCornerID

  implicit none

#define IMPURE_ELEMENTAL
#define CONTIGUOUS_POINTER pointer
#define famrex_abort Driver_abortFlash
#include "constants.h"
  private

  public :: famrex_multivab_build, famrex_multivab_destroy
  public :: famrex_mviter_build, famrex_mviter_destroy

  integer,parameter :: ndims=N_DIM
  integer,parameter :: amrex_real=kind(1.0)

  type, public   :: famrex_multivab
     logical               :: owner = .false.
     type   (c_ptr)        :: p     =  c_null_ptr
     integer(c_int)        :: nc    =  0
     integer(c_int)        :: ng    =  0
     integer               :: nodetypes = LEAF
     integer               :: gds = CENTER
     integer               :: dur = VD_DUR_PERM
     integer               :: lev = -1
     ! different variants for cell index numbering:
     !  1 = normal FLASH convention:     per block, leftmost guard cell = 1
     !  0 = zero-based FLASH convention: per block, leftmost guard cell = 0
     ! -1 = global convention for a refinement level: leftmost guard cell = 1
     ! -2 = global convention for a refinement level: leftmost inner cell = 1
     integer               :: cellIdxBase = -2
     INTEGER :: dm
   contains
     generic   :: assignment(=) => famrex_multivab_assign
     procedure :: dataPtr       => famrex_multivab_dataptr
     procedure, private :: &
          famrex_multivab_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: famrex_multivab_destroy
#endif
  end type famrex_multivab

  interface famrex_multivab_build
     module procedure famrex_multivab_build
     module procedure famrex_multivab_build_l
  end interface famrex_multivab_build

  type, public :: famrex_mviter
     type(c_ptr)      :: p       = c_null_ptr
     type(famrex_multivab),pointer :: mvp      => NULL()
     integer ,private :: counter = -1
     integer ,private :: cur     = 0
     integer ,private :: nodetypes   = ALL_BLKS
     integer ,private :: gds     = CENTER
     integer ,private :: lev=-1
   contains
     generic   :: assignment(=)    => famrex_mviter_assign  ! will abort if called
     procedure :: clear            => famrex_mviter_clear
     procedure :: next             => famrex_mviter_next
     procedure :: tilebox          => famrex_mviter_tilebox
     procedure :: localIndex       => famrex_mviter_localIndex
     procedure, private :: famrex_mviter_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: famrex_mviter_destroy
#endif
  end type famrex_mviter

  interface famrex_mviter_build
     module procedure famrex_mviter_build_r
  end interface famrex_mviter_build


contains

  subroutine famrex_multivab_build (mv, nodetypes, gds, dm, nc, lev, ng, nodal)
    type(famrex_multivab), intent(inout) :: mv

    integer, intent(in ) :: nodetypes
    integer, intent(in )   :: gds
    INTEGER,intent(in )   :: dm
    integer, intent(in) :: nc
    integer, intent(in), optional :: lev
    integer, intent(in), optional :: ng
    logical, intent(in), optional :: nodal(*)
    integer :: inodal(3), dir
    mv%owner = .true.
    mv%nc = nc
    mv%ng = 0
    if (present(ng)) mv%ng = ng
    inodal = 0
    if (present(nodal)) then
       do dir = 1, ndims
          if (nodal(dir)) inodal(dir) = 1
       end do
    end if
    mv%nodetypes = nodetypes
    if (present(lev)) mv%lev = lev
    mv%gds = gds
    mv%dm = dm
  end subroutine famrex_multivab_build

  subroutine famrex_multivab_build_l (mv, blkList, gds, dm, nc, ng, nodal)
    type(famrex_multivab), intent(inout) :: mv

    integer, intent(in ) :: blkList(:)
    integer, intent(in )   :: gds
    INTEGER,intent(in )   :: dm
    integer, intent(in) :: nc
    integer, intent(in), optional :: ng
    logical, intent(in), optional :: nodal(*)
    integer :: inodal(3), dir
    integer :: numBlks
    mv%owner = .FALSE.
    numBlks = size(blkList)
    mv%nc = nc
    mv%ng = 0
    if (present(ng)) mv%ng = ng
    inodal = 0
    if (present(nodal)) then
       do dir = 1, ndims
          if (nodal(dir)) inodal(dir) = 1
       end do
    end if
    mv%gds = gds
    mv%dm = dm
  end subroutine famrex_multivab_build_l

  IMPURE_ELEMENTAL subroutine famrex_multivab_destroy (this)
    type(famrex_multivab), intent(inout) :: this
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine famrex_multivab_destroy

  subroutine famrex_multivab_assign (dst, src)
    class(famrex_multivab), intent(inout) :: dst
    type (famrex_multivab), intent(in   ) :: src
    call famrex_multivab_destroy(dst)
    dst%owner = .false.
    dst%p     = src%p
    dst%nc    = src%nc
    dst%ng    = src%ng
    dst%gds   = src%gds
    dst%dur   = src%dur
    dst%dm    = src%dm
    dst%nodetypes = src%nodetypes
    dst%lev   = src%lev
  end subroutine famrex_multivab_assign

  function famrex_multivab_dataPtr (this, mvi) result(dp)
    class(famrex_multivab), intent(in) :: this
    type(famrex_mviter), intent(in) :: mvi
    real(amrex_real), CONTIGUOUS_POINTER, dimension(:,:,:,:) :: dp
    type(c_ptr) :: cp
    real(amrex_real), CONTIGUOUS_POINTER :: fp(:,:,:,:)
    integer, dimension(LOW:HIGH,MDIM) :: blkLim ,blkLimGC
    integer :: cur
    integer :: cornerID(MDIM), strideUnused(MDIM)
    integer(c_int) :: n(4)
    type(famrex_box) :: bx
    cur = mvi%localIndex()
    call Grid_getBlkIndexLimits(cur,blkLim,blkLimGC)
    call Grid_getBlkPtr(cur,fp,this%gds)

    bx%lo = blkLimGC(LOW ,:)                   !DEV: or blkLimGC(LOW ,:) ?
    if (this%cellIdxBase==-1) then
       call Grid_getBlkCornerID(cur,cornerID,strideUnused)
       cornerID = (cornerID - 1) / 2**(lrefine_max-lrefine(cur)) + 1
       bx%lo(:) = bx%lo(:) - 1 + cornerID(:)
    else if (this%cellIdxBase==-2) then
       call Grid_getBlkCornerID(cur,cornerID,strideUnused)
       cornerID = (cornerID - 1) / 2**(lrefine_max-lrefine(cur)) + 1
       bx%lo(:) = bx%lo(:) - 1 + cornerID(:)
       bx%lo(1:ndims) = bx%lo(1:ndims) - NGUARD
    else if (this%cellIdxBase==0) then
       bx%lo(:) = bx%lo(:) - 1
    end if
#ifdef INDEXREORDER
    dp(bx%lo(1):,bx%lo(2):,bx%lo(3):,1:) => fp
#else
    dp(1:,bx%lo(1):,bx%lo(2):,bx%lo(3):) => fp
#endif
  end function famrex_multivab_dataPtr

!------ imultivab routines removed ------!

!------ MFIter routines ------!

  subroutine famrex_mviter_build_r (mvi, mv, tiling)
    type(famrex_mviter) :: mvi
    type(famrex_multivab),target, intent(in ) :: mv
    logical, intent(in), optional :: tiling
    logical :: ltiling
    integer(c_int) :: t
    ltiling = .false.;  if (present(tiling)) ltiling = tiling
    if (ltiling) then
       t = 1
    else
       t = 0
    end if
    mvi%counter = 0
    mvi%mvp => mv
    mvi%nodetypes = mv%nodetypes
    mvi%gds = mv%gds
    mvi%lev = mvi%lev
  end subroutine famrex_mviter_build_r


  subroutine famrex_mviter_destroy (this)
    type(famrex_mviter) :: this
    call this%clear()
  end subroutine famrex_mviter_destroy

  subroutine famrex_mviter_assign (dst, src)
    class(famrex_mviter), intent(inout) :: dst
    type (famrex_mviter), intent(in   ) :: src
    ! No way to disable it at compile time, so ...
    call famrex_abort("famrex_mviter assignment is disabled")
  end subroutine famrex_mviter_assign

  subroutine famrex_mviter_clear (this)
    class(famrex_mviter) :: this
!!    print*,'CALLED famrex_mviter_clear!'
    this%counter = -1
    this%cur     = -1
    this%p = c_null_ptr
  end subroutine famrex_mviter_clear

  logical function famrex_mviter_next (this)
    use Grid_interface,ONLY: Grid_blockMatch
    class(famrex_mviter),intent(inout) :: this
    integer(c_int) :: isvalid
    logical        :: fisvalid
    integer        :: cur, i
    this%counter = this%counter + 1

    cur = this%cur
    do i=cur+1,lnblocks
       if (Grid_blockMatch(i,this%nodetypes,this%lev)) EXIT !right nodetype?
    end do
    this%cur = i
    fisvalid = (i .LE. lnblocks)
    famrex_mviter_next = fisvalid

  end function famrex_mviter_next

  function famrex_mviter_tilebox (this) result (bx)
    class(famrex_mviter), intent(in) :: this
    type(famrex_box) :: bx
    integer, dimension(LOW:HIGH,MDIM) :: blkLim ,blkLimGC
    integer :: inodal(3)
    integer :: cornerID(MDIM), strideUnused(MDIM)
    inodal = 0
    call Grid_getBlkIndexLimits(this%cur,blkLim,blkLimGC)
    bx%lo = blkLim(LOW ,:)                   !DEV: or blkLimGC(LOW ,:) ?
    bx%hi = blkLim(HIGH,:)                   !DEV: or blkLimGC(HIGH,:) ?
    if (this%mvp%cellIdxBase==-1) then
       call Grid_getBlkCornerID(this%cur,cornerID,strideUnused)
       cornerID = (cornerID - 1) / 2**(lrefine_max-lrefine(this%cur)) + 1
       bx%lo(:) = bx%lo(:) - 1 + cornerID(:)
       bx%hi(:) = bx%hi(:) - 1 + cornerID(:)
    else if (this%mvp%cellIdxBase==-2) then
       call Grid_getBlkCornerID(this%cur,cornerID,strideUnused)
       cornerID = (cornerID - 1) / 2**(lrefine_max-lrefine(this%cur)) + 1
       bx%lo(:) = bx%lo(:) - 1 + cornerID(:)
       bx%hi(:) = bx%hi(:) - 1 + cornerID(:)
       bx%lo(1:ndims) = bx%lo(1:ndims) - NGUARD
       bx%hi(1:ndims) = bx%hi(1:ndims) - NGUARD
    else if (this%mvp%cellIdxBase==0) then
       bx%lo(:) = bx%lo(:) - 1
       bx%hi(:) = bx%hi(:) - 1
    end if
    where (inodal .ne. 0) bx%nodal = .true.  ! note default is false
  end function famrex_mviter_tilebox


  function famrex_mviter_localIndex (this) result (cidx)
    class(famrex_mviter), intent(in) :: this
    integer :: cidx
    cidx = this%cur
  end function famrex_mviter_localIndex

end module famrex_multivab_module


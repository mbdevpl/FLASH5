
module famrex_box_module

  use,intrinsic :: iso_c_binding

  implicit none

  private

  integer,parameter :: ndims=N_DIM

  type, public :: famrex_box
     integer, dimension(3) :: lo    = 1
     integer, dimension(3) :: hi    = 1
     logical, dimension(3) :: nodal = .false.
   contains
     generic   :: grow     => famrex_box_grow_s, famrex_box_grow_v
     procedure, private :: famrex_box_grow_s
     procedure, private :: famrex_box_grow_v
  end type famrex_box

  interface famrex_box
     module procedure famrex_box_build
  end interface famrex_box

contains

  function famrex_box_build (lo, hi, nodal) result(bx)
    integer, intent(in) :: lo(ndims), hi(ndims)
    logical, intent(in), optional :: nodal(ndims)
    type(famrex_box) :: bx
    bx%lo(1:ndims) = lo(1:ndims)
    bx%hi(1:ndims) = hi(1:ndims)
    if (present(nodal)) bx%nodal(1:ndims) = nodal(1:ndims)
  end function famrex_box_build

  subroutine famrex_box_grow_s (this, i)
    class(famrex_box), intent(inout) :: this
    integer, intent(in) :: i
    this%lo = this%lo - i
    this%hi = this%hi + i
  end subroutine famrex_box_grow_s

  subroutine famrex_box_grow_v (this, i)
    class(famrex_box), intent(inout) :: this
    integer, intent(in) :: i(ndims)
    this%lo(1:ndims) = this%lo(1:ndims) - i(1:ndims)
    this%hi(1:ndims) = this%hi(1:ndims) + i(1:ndims)
  end subroutine famrex_box_grow_v

end module famrex_box_module


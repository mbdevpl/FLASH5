subroutine ImBound_setData()
  use IncompNS_data, only : ins_invRe
  use ImBound_data,  only : ib_nu 

  implicit none

  ib_nu = ins_invRe

  return
end subroutine ImBound_setData

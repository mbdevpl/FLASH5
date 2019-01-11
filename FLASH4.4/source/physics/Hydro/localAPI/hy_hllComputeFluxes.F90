Subroutine hy_hllComputeFluxes(tileLimits, Uin, plo, flX, flY, flZ, loFl, del, dt)
  implicit none

  integer, intent(IN)  :: tileLimits(LOW:HIGH, 1:MDIM)
  integer, intent(IN)  :: plo(*)
  real,    intent(IN)  :: UIN(plo(1):,plo(2):,plo(3):,plo(4):)  !CAPITALIZATION INTENTIONAL!
  integer, intent(IN)  :: loFl(*)
  real,    intent(OUT) :: FLX(loFl(1):,loFl(2):,loFl(3):,loFl(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(OUT) :: FLY(loFl(1):,loFl(2):,loFl(3):,loFl(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(OUT) :: FLZ(loFl(1):,loFl(2):,loFl(3):,loFl(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(IN)  :: del(1:MDIM)
  real,    intent(IN)  :: dt

  flX(:, :, :, :) = 0.0
  flY(:, :, :, :) = 0.0
  flZ(:, :, :, :) = 0.0
End Subroutine hy_hllComputeFluxes


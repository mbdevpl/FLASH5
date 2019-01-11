#include "constants.h"

Subroutine hy_hllComputeFluxes(tileLimits, Uin, plo, flX, flY, flZ, loFl, del, dt)
  implicit none

  integer, intent(IN)  :: tileLimits(LOW:HIGH, 1:MDIM)
  integer, intent(IN)  :: plo(*)
  real,    intent(IN)  :: Uin(plo(1):,plo(2):,plo(3):,plo(4):)
  integer, intent(IN)  :: loFl(*)
  real,    intent(OUT) :: flX(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
  real,    intent(OUT) :: flY(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
  real,    intent(OUT) :: flZ(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
  real,    intent(IN)  :: del(1:MDIM)
  real,    intent(IN)  :: dt

  flX(:, :, :, :) = 0.0
  flY(:, :, :, :) = 0.0
  flZ(:, :, :, :) = 0.0
End Subroutine hy_hllComputeFluxes


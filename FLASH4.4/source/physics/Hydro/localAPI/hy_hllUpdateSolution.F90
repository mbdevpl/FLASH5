Subroutine hy_hllUpdateSolution(tileLimits, Uin, plo, Uout, flX, flY, flZ, loFl, del, dt)
  implicit none
  
  integer, intent(IN)  :: tileLimits(LOW:HIGH, 1:MDIM)
  integer, intent(IN)  :: plo(*)
  real,    intent(IN)  :: Uin(plo(1):,plo(2):,plo(3):,plo(4):)
  real,    intent(OUT) :: Uout(plo(1):,plo(2):,plo(3):,plo(4):)
  integer, intent(IN)  :: loFl(*)
  real,    intent(IN)  :: flX(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
  real,    intent(IN)  :: flY(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
  real,    intent(IN)  :: flZ(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
  real,    intent(IN)  :: del(1:MDIM)
  real,    intent(IN)  :: dt
    
  Uout(:, :, :, :) = 0.0
End Subroutine hy_hllUpdateSolution

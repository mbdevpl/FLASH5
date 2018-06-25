! STUB: Subroutine outtotecplot
!
! Subroutine to write out to Paraview in hdf5 format.
!
! ---------------------------------------------------------------------------


  subroutine outtoParaview(filename,ib,mype,time,dt,istep,count, &
                          timer,blockList,blockCount,firstfileflag)

  integer, intent(in)   :: mype,istep,count,firstfileflag,ib
  integer, intent(in)   :: blockCount
  integer, intent(in)   :: blockList(MAXBLOCKS)
  real, intent(in)      :: time,dt,timer
  character(len=100),intent(in)  :: filename

end subroutine outtoParaview

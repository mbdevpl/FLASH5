Module sim_interface
  interface
     subroutine sim_initBlockAnalytical(block)
       use block_metadata, ONLY: block_metadata_t
       type(block_metadata_t),intent(IN) :: block
     end subroutine sim_initBlockAnalytical
  end interface
end Module sim_interface

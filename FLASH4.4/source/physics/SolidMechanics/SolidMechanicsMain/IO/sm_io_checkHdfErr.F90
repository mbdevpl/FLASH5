subroutine sm_io_checkHdfErr(hdferr, info_message)
    Use Driver_Interface, only: Driver_abortFlash
    integer, intent(in) :: hdferr
    CHARACTER, intent(in) :: info_message*(*)
    
    if( hdferr < 0 ) then
        write(*,*) ''
        write(*,*) 'Error: ', hdferr
        write(*,*) '          '//info_message
        write(*,*) '          ending program.'
        
        !Old method:
        ! stop
        
        ! better method
        call Driver_abortFlash("Recheck input file for structure.")
    endif
    
end


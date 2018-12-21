module sim_interface

    interface
        subroutine sim_writeDataPoints(initData, tileDesc, points, values)
            use flash_tile, ONLY : flash_tile_t
            implicit none
            real,               intent(IN), pointer :: initData(:, :, :, :)
            type(flash_tile_t), intent(IN)          :: tileDesc
            real,               intent(IN)          :: points(:, :)
            real,               intent(IN)          :: values(:)
        end subroutine sim_writeDataPoints
    end interface

    interface
        subroutine sim_collectLeaves
            implicit none
        end subroutine sim_collectLeaves
    end interface 

    interface
        subroutine sim_printLeaves(title)
            implicit none
            character(*), intent(IN)    :: title
        end subroutine sim_printLeaves
    end interface 

    interface
        subroutine sim_advance(step, points, values, set_msg, leaf_msg)
            implicit none
            integer,      intent(IN)    :: step
            real,         intent(IN)    :: points(:, :)
            real,         intent(IN)    :: values(:)
            character(*), intent(IN)    :: set_msg
            character(*), intent(IN)    :: leaf_msg
        end subroutine sim_advance
    end interface 

end module sim_interface


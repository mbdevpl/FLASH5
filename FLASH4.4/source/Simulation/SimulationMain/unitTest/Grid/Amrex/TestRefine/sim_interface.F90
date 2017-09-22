module sim_interface

    interface
        subroutine sim_writeDataPoints(initData, block, points, values)
            use block_metadata, ONLY : block_metadata_t
            implicit none
            real,                   intent(IN), pointer :: initData(:, :, :, :)
            type(block_metadata_t), intent(IN)          :: block
            real,                   intent(IN)          :: points(:, :)
            real,                   intent(IN)          :: values(:)
        end subroutine sim_writeDataPoints
    end interface

    interface
        subroutine sim_printLeaves(title)
            implicit none
            character(*), intent(IN) :: title
        end subroutine sim_printLeaves
    end interface 

end module sim_interface


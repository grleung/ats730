module hello_world

    ! This module 'hello_world' contains the subroutine 'hello_world_subroutine'

    implicit none

    contains

    subroutine hello_world_subroutine()

        ! Use the 'constants' module in order to get our favorite number
        use constants

        implicit none

        ! Print text and favorite number to the terminal
        print*,'Hello World'
        print*,'My favorite number is: '
        print*,fav_number

    end subroutine hello_world_subroutine

end module hello_world
program MAC

    ! This file contains the MAC program

    ! Get the initialization subroutines
    use grid, only: init_grid
    use base_state, only: init_base_state

    implicit none

    ! First, let's setup
    call init_grid 
    
    ! Call base state initialization
    call init_base_state

end program MAC
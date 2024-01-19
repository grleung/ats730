program MAC

    ! This file contains the MAC program

    ! Get the initialization subroutines
    use grid, only: init_grid
    use base_state, only: init_base_state
    use io, only: write_base_state

    implicit none

    ! First, let's setup
    call init_grid 
    
    ! Call base state initialization
    call init_base_state

    ! Write output to a simple text file
    call write_base_state('hw1_output.txt')
 
end program MAC
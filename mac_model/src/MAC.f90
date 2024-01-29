program MAC

    ! This file contains the MAC program

    ! Get the initialization subroutines
    use grid, only: init_grid
    use base_state, only: init_base_state
    use cape, only: calculate_parcel_cape
    use io, only: write_parcel_traj !write_base_state, 

    implicit none

    ! First, let's setup
    call init_grid 
    
    ! Call base state initialization
    call init_base_state

    ! Calculate parcel CAPE
    call calculate_parcel_cape

    ! Write output to a simple text file
    !call write_base_state('hw1_output.txt')
    !call write_parcel_traj('hw2_output.txt')
 
end program MAC
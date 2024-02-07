program MAC

    ! This file contains the MAC program

    ! Get the initialization subroutines
    use grid, only: init_grid
    use run_constants, only: base_outpath, parcel_outpath
    use mem, only: allocate_mem, deallocate_mem
    use base_state, only: init_base_state
    use cape, only: calculate_parcel_cape
    use io, only: read_namelist, write_parcel_traj, write_base_state,write_current_state
    use initial_perturb, only: init_perturb

    implicit none
    
    call read_namelist
    
    call allocate_mem

    ! First, let's setup
    call init_grid 
    
    ! Call base state initialization
    call init_base_state
    ! Write output to a simple text file
    call write_base_state

    ! Calculate parcel CAPE
    call calculate_parcel_cape
    call write_parcel_traj

    call init_perturb

    call write_current_state

    call deallocate_mem

end program MAC
program MAC

    ! This file contains the MAC program

    ! Get the initialization subroutines
    use grid, only: init_grid
    use run_constants, only: base_outpath, parcel_outpath,dt,endt
    use mem, only: allocate_mem, deallocate_mem
    use base_state, only: init_base_state
    use cape, only: calculate_parcel_cape
    use io, only: read_namelist, write_parcel_traj, write_base_state,write_current_state
    use initial_perturb, only: init_perturb
    use advection, only: advect
    use model_vars, only: it

    implicit none

    integer :: nt !  total number of timesteps
   
    call read_namelist
    nt = int(endt/dt)
    
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

    ! initialize perturbation from base state
    call init_perturb

    !each timestep
    do it=1,nt
        call advect
        call write_current_state
    enddo

    call deallocate_mem

end program MAC
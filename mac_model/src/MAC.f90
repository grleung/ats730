program MAC

    ! This file contains the MAC program

    ! Get the initialization subroutines
    use grid, only: init_grid
    use run_constants, only: base_outpath, parcel_outpath,dt,endt,outfreq
    use mem, only: allocate_mem, deallocate_mem
    use base_state, only: init_base_state
    use cape, only: calculate_parcel_cape
    use io, only: read_namelist, write_parcel_traj, write_base_state,write_current_state
    use initial_perturb, only: init_perturb
    use advection, only: advect
    use model_vars, only: it
    use solve_prog, only: tendencies
    use boundaries, only: enforce_bounds_x, enforce_bounds_y,enforce_bounds_z
    use timestep, only: step_time
    use cfl, only: check_cfl

    implicit none

    integer :: nt !  total number of timesteps
   
    call read_namelist
    nt = int(endt/dt)
    
    call allocate_mem
    print*, "memory allocated" 
    
    ! First, let's setup
    call init_grid 
    print*, "grid initialized" 

    ! Call base state initialization
    call init_base_state
    print*, "base state initialized"
    ! Write output to a simple text file
    call write_base_state
    print*, "base state written"

    ! Calculate parcel CAPE
    call calculate_parcel_cape
    call write_parcel_traj
    print*, "parcel trajectory calculated" 

    ! initialize perturbation from base state
    call init_perturb
    print*, "perturbation initialized" 
    
    call write_current_state
    print*, "initial timestep written"

    !each timestep
    do it=1,nt
        !call advect
        call enforce_bounds_z
        call enforce_bounds_y
        call enforce_bounds_x
        print*, "bounds enforced"
        call tendencies     
        print*, "tendencies calculated"
        call enforce_bounds_z
        call enforce_bounds_y
        call enforce_bounds_x
        print*, "bounds enforced"
        call step_time
        if (modulo(it,outfreq)==0) then
            call write_current_state
            print*, "timestep written"
        endif

        call check_cfl
        print*,'cfl ok'
    
        print*,it
    enddo

    call deallocate_mem

end program MAC
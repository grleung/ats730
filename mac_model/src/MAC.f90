program MAC

    ! This file contains the MAC program

    ! Get the initialization subroutines
    use grid, only: init_grid
    use run_constants, only: base_outpath, parcel_outpath,dt,endt,outfreq
    use mem, only: allocate_mem, deallocate_mem
    use base_state, only: init_base_state
    use cape, only: calculate_parcel_cape
    use io, only: read_namelist, write_parcel_traj, write_base_state,write_current_state,write_current_micro
    use initial_perturb, only: init_perturb
    use advection, only: advect
    use model_vars, only: it
    use solve_prog, only: tendencies
    use boundaries, only: enforce_bounds_x, enforce_bounds_z
    use timestep, only: step_time
    use aerosol, only: init_aerosol,microphysics_driv

    implicit none

    integer :: nt !  total number of timesteps
   
    call read_namelist
    nt = int(endt/dt)
    
    call allocate_mem
    print*,'allocated mem'

    ! First, let's setup
    call init_grid 
    print*,'grid init'
    
    ! Call base state initialization
    call init_base_state
    print*,'base init'
    ! Write output to a simple text file
    call write_base_state
    print*,'base written'

    call init_aerosol
    !print*,'aerosol init'

    ! Calculate parcel CAPE
    !call calculate_parcel_cape
    !call write_parcel_traj

    ! initialize perturbation from base state
    call init_perturb
    print*,'pert init'

    call write_current_state
    print*,'initial timestep written'

    !each timestep
    do it=1,nt
        !call advect
        call enforce_bounds_z
        call enforce_bounds_x
        print*,'bounds enforced'
        call tendencies 
        print*,'done tend'
        call microphysics_driv
        print*,'done micro'
        call enforce_bounds_z
        call enforce_bounds_x
        print*,'bounds enforced'
        call step_time
        print*,'step time'

        if (modulo(it,outfreq)==0) then
            call write_current_state
            print*,'state write'
            call write_current_micro
            print*,'micro write'
        endif
    
        print*,it
    enddo

    call deallocate_mem

end program MAC
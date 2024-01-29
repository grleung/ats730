module grid
    implicit none

    contains

    subroutine init_grid()

        use grid_constants, only: nz, dz0, dzrat
        use model_vars, only: dzn, zun, zwn

        implicit none

        integer :: iz ! counter for z-coordinate

        ! print statement to check what's happening
        ! print*,'nz'
        ! print*,nz

        ! Set up the delta z's for each vertical level from 1 to nz
        ! depending on the initial delta z (dz0) and stretch ratio (dzrat)
        ! right now everything is uniform so this doesn't do anything 
        ! but set up this way in case we want a different vertical grid later
        dzn(1) = dz0
        do iz = 2,nz
            dzn(iz) = dzn(iz-1)*dzrat
        enddo
        
        ! Set up the vertical heights of the w grid
        zwn(1) = -dzn(1) ! this is a fictitious point that is one grid spacing below surface
        zwn(2) = 0. ! model surface 
        do iz = 3,nz
            zwn(iz) = zwn(iz-1)+dzn(iz)
        enddo

        ! Set up the vertical heights of the scalar/u grid
        ! do this by interpolating between the heights of the w grid
        ! note that ztn(1) and ztn(nz) are fictitious points anyway
        do iz = 1,nz
            zun(iz) = (zwn(iz)+zwn(iz+1))/2
        enddo
        zun(nz) = zun(nz-1)+dzn(nz-1) 
    end subroutine init_grid
end module grid

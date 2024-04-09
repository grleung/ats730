module grid
    implicit none

    contains

    subroutine init_grid()

        use run_constants, only: nz, dz0, dzrat, nx, dx, ny, dy
        use model_vars, only: dzn, zsn, zmn, dxn, xsn, xmn, dyn, ysn, ysn, ymn

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate
        integer :: iy ! counter for y-coordinate

        ! set up z-coordinate

        ! Set up the delta z's for each vertical level from 1 to nz
        ! depending on the initial delta z (dz0) and stretch ratio (dzrat)
        ! right now everything is uniform so this doesn't do anything 
        ! but set up this way in case we want a different vertical grid later
        dzn(1) = dz0
        do iz = 2,nz
            dzn(iz) = dzn(iz-1)*dzrat
        enddo
        
        ! Set up the vertical heights of the w grid
        zmn(1) = -dzn(1) ! this is a fictitious point that is one grid spacing below surface
        zmn(2) = 0. ! model surface 
        do iz = 3,nz
            zmn(iz) = zmn(iz-1)+dzn(iz)
        enddo

        ! Set up the vertical heights of the scalar/u grid
        ! do this by interpolating between the heights of the w grid
        ! note that ztn(1) and ztn(nz) are fictitious points anyway
        do iz = 1,nz
            zsn(iz) = (zmn(iz)+zmn(iz+1))*0.5
        enddo
        zsn(nz) = zsn(nz-1)+dzn(nz-1) 

        ! set up x-coordinate

        ! grid spacing is just evenly spaced in horizontal
        do ix = 1,nx
            dxn(ix) = dx
        enddo

        ! set up horizontal location of u grid, which is staggered with the scalar grid
        xmn(1) = -dxn(1) !this is a fictitious point that is one grid spacing to the west of the boundary (may change for PBC)
        xmn(2) = 0.
        do ix = 3,nx
            xmn(ix) = xmn(ix-1)+dxn(ix)
        enddo

        ! set up horizontal location of scalar/w/v grid
        ! by interpolating between u-grid locations
        do ix = 1,nx
            xsn(ix) = (xmn(ix)+xmn(ix+1))*0.5
        enddo
        xsn(nx) = xsn(nx-1)+dxn(nx-1)

        ! set up y-coordinate

        ! grid spacing is just evenly spaced in horizontal
        do iy = 1,ny
            dyn(iy) = dy
        enddo

        ! set up horizontal location of v grid, which is staggered with the scalar grid
        ymn(1) = -dyn(1) !this is a fictitious point that is one grid spacing to the west of the boundary (may change for PBC)
        ymn(2) = 0.
        do iy = 3,ny
            ymn(iy) = ymn(iy-1)+dyn(iy)
        enddo

        ! set up horizontal location of scalar/w/u grid
        ! by interpolating between v-grid locations
        do iy = 1,ny
            ysn(iy) = (ymn(iy)+ymn(iy+1))*0.5
        enddo
        ysn(ny) = ysn(ny-1)+dyn(ny-1)


        
    end subroutine init_grid
end module grid

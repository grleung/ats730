module initial_perturb

    ! This module contains the subroutine to set initial perturbation from the base state

    implicit none

    contains

    subroutine init_perturb()
        use constants, only: g,cp,trigpi
        use run_constants, only: nz, nx, ny,dz0,radx,rady,radz,amp,zcnt,xcnt,ycnt,pbc_x,pbc_y,pbc_z,pert_wind
        use model_vars, only:zsn,xsn,ysn,thb,thvb,rhoub,thp,pip,pp,up,vp,wp,rvp,rcp,rrp
        use boundaries, only: enforce_bounds_x,enforce_bounds_y,enforce_bounds_z

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate
        integer :: iy ! counter for y-coordinate
        real :: rad ! distance from center of perturbation

        do iz = 2, nz-1
            do iy = 2, ny-1
                do ix = 2, nx-1
                    rad = (((zsn(iz)-zcnt)/radz)**2 + ((xsn(ix)-xcnt)/radx)**2 + ((ysn(iy)-ycnt)/rady)**2)**(0.5)
                    
                    ! construct perturbation
                    if (rad<=1) then
                        if (pert_wind) then 
                            up(ix,iy,iz,2) = 0.5 * amp * (COS(rad*trigpi) + 1)
                        else
                            thp(ix,iy,iz,2) = 0.5 * amp * (COS(rad*trigpi) + 1)
                        endif
                    else
                        if (pert_wind) then 
                            up(ix,iy,iz,2) = 0.
                        else
                            thp(ix,iy,iz,2) = 0.
                        endif
                    endif

                    ! set initial values for other variables
                    vp(ix,iy,iz,2) = 0.
                    wp(ix,iy,iz,2) = 0.
                    rvp(ix,iy,iz,2) = 0.
                    rcp(ix,iy,iz,2) = 0.
                    rrp(ix,iy,iz,2) = 0.

                enddo
            enddo
        enddo
        
        !construct pressure perturbation field to be in balance with theta perturbation
        do iy = 2, ny-1
            do ix = 2, nx - 1
                ! assume that perturbation pressure is 0 near top of model domain
                pip(ix,iy,nz-1,2) = 0.

                ! integrate downward: dpi/dz = g/cp * thp/thb^2
                do iz = nz-2,2,-1
                    pip(ix,iy,iz,2) = pip(ix,iy,iz+1,2) - ((g/cp)*(thp(ix,iy,iz,2)/(thb(iz)**2))*dz0)
                enddo 
            enddo 
        enddo

        call enforce_bounds_z
        call enforce_bounds_y
        call enforce_bounds_x

        ! for now, set the past to be the same as present timestep
        ! since we don't have any information about the past
        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    thp(ix,iy,iz,1) = thp(ix,iy,iz,2)
                    pip(ix,iy,iz,1) = pip(ix,iy,iz,2)
                    up(ix,iy,iz,1) = up(ix,iy,iz,2)
                    vp(ix,iy,iz,1) = vp(ix,iy,iz,2)
                    wp(ix,iy,iz,1) = wp(ix,iy,iz,2)
                    rvp(ix,iy,iz,1) = rvp(ix,iy,iz,2)
                    rcp(ix,iy,iz,1) = rcp(ix,iy,iz,2)
                    rrp(ix,iy,iz,1) = rrp(ix,iy,iz,2)
                enddo
            enddo
        enddo
        
        !calculate perturbation pressure in Pa
        do iz = 1,nz
            do iy=1,ny
                do ix = 1,nx
                    pp(ix,iy,iz) = pip(ix,iy,iz,2)*cp*rhoub(iz)*thvb(iz)
                enddo
            enddo
        enddo

    end subroutine init_perturb

end module initial_perturb
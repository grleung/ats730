module initial_perturb

    ! This module contains the subroutine to set initial perturbation from the base state

    implicit none

    contains

    subroutine init_perturb()
        use constants, only: g,cp,trigpi
        use run_constants, only: nz, nx, dz0,radx,radz,amp,zcnt,xcnt, pbc_x, pbc_z,pert_wind
        use model_vars, only:zsn, xsn,thb,rvb,thvb,rhoub,thp,pip,pp,up,wp,rvp,rcp,rrp,thvp
        use boundaries, only: enforce_bounds_x,enforce_bounds_z
        use thermo_functions, only:calc_thv

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate
        real :: rad ! distance from center of perturbation

        do iz = 2, nz-1
            do ix = 2, nx-1
                rad = (((zsn(iz)-zcnt)/radz)**2 + ((xsn(ix)-xcnt)/radx)**2)**(0.5)
                
                ! construct perturbation
                if (rad<=1) then
                    if (pert_wind) then 
                        up(ix,iz,2) = 0.5 * amp * (COS(rad*trigpi) + 1)
                    else
                        thp(ix,iz,2) = 0.5 * amp * (COS(rad*trigpi) + 1)
                    endif
                endif
            enddo
        enddo
        
        !construct pressure perturbation field to be in balance with theta perturbation
        do ix = 2, nx - 1
            ! assume that perturbation pressure is 0 near top of model domain
            pip(ix,nz-1,2) = 0.

            ! integrate downward: dpi/dz = g/cp * thp/thb^2
            do iz = nz-2,2,-1
                pip(ix,iz,2) = pip(ix,iz+1,2) - ((g/cp)*(thp(ix,iz,2)/(thb(iz)**2))*dz0)
            enddo 
        enddo 

        do ix = 2,nx-1
            do iz = 2,nz-1
                rvp(ix,iz,2) = 0.
                rcp(ix,iz,2) = 0.
                rrp(ix,iz,2) = 0.
                up(ix,iz,2) = 0.
                wp(ix,iz,2) = 0.
            enddo
        enddo

        call enforce_bounds_z
        call enforce_bounds_x

        ! for now, set the past to be the same as present timestep
        ! since we don't have any information about the past
        do iz = 1, nz
            do ix = 1, nx
                thp(ix,iz,1) = thp(ix,iz,2)
                pip(ix,iz,1) = pip(ix,iz,2)
                up(ix,iz,1) = up(ix,iz,2)
                wp(ix,iz,1) = wp(ix,iz,2)
                rvp(ix,iz,1) = rvp(ix,iz,2)
                rcp(ix,iz,1) = rcp(ix,iz,2)
                rrp(ix,iz,1) = rrp(ix,iz,2)
            enddo
        enddo
        
        !calculate perturbation pressure in Pa
        do iz = 1,nz
            do ix = 1,nx
                pp(ix,iz) = pip(ix,iz,2)*cp*rhoub(iz)*thvb(iz)
                thvp(ix,iz) = calc_thv(thb(iz)+thp(ix,iz,2),rvb(iz)+rvp(ix,iz,2))
            enddo
        enddo

    end subroutine init_perturb

end module initial_perturb
module initial_perturb

    ! This module contains the subroutine to set initial perturbation from the base state

    implicit none

    contains

    subroutine init_perturb()
        use constants, only: g,cp
        use run_constants, only: nz, nx, radx,radz,amp,zcnt,xcnt
        use thermo_functions, only: calc_thv, calc_rsat, calc_satfrac
        use model_vars, only:zsn, xsn,thb, thvb, thp,pip,pp,rhoub

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate
        integer :: it ! counter for time coordinate
        real :: rad ! distance from center of perturbation

        do iz = 2, nz-1
            do ix = 2, nx-1
                rad = (((zsn(iz)-zcnt)/radz)**2 + ((xsn(ix)-xcnt)/radx)**2)**(0.5)
                
                ! construct temperature perturbation
                if (rad<=1) then
                    thp(ix,iz,1) = 0.5 * amp * (COS(rad*4*ATAN(1.)) + 1)
                endif
            enddo
        enddo
        
        !construct pressure perturbation field to be in balance with theta perturbation
        do ix = 2, nx - 1
            ! assume that perturbation pressure is 0 near top of model domain
            pip(ix,nz-1,1) = 0.

            ! integrate downward: dpi/dz = g/cp * thp/thb^2
            do iz = nz-2,2,-1
                pip(ix,iz,1) = pip(ix,iz+1,1) + (g/cp)*(thp(ix,iz,1)/thb(iz)**2)*(zsn(iz-1)-zsn(iz))
            enddo 
        enddo 

        ! set fictitous points to be equal to first real point
        do iz = 2,nz-2
            pip(1,iz,1) = pip(2,iz,1)
            pip(nx,iz,1) = pip(nx-1,iz,1)
            thp(1,iz,1) = thp(2,iz,1)
            thp(nx,iz,1) = thp(nx-1,iz,1)
        enddo

        do ix = 2,nx-2
            pip(ix,1,1) = pip(ix,2,1)
            pip(ix,nz,1) = pip(ix, nz-1,1)
            thp(ix,1,1) = thp(ix,2,1)
            thp(ix,nz,1) = thp(ix,nz-1,1)
        enddo

        
        ! for now, set the present to be the same as past timestep
        do iz = 1, nz
            do ix = 1, nx
                thp(ix,iz,2) = thp(ix,iz,1)
                pip(ix,iz,2) = pip(ix,iz,1)
            enddo
        enddo
        
        !calculate perturbation pressure in Pa
        do iz = 1,nz
            do ix = 1,nx
                pp(ix,iz) = pip(ix,iz,1)*cp*rhoub(iz)*thvb(iz)
            enddo
        enddo

    end subroutine init_perturb

end module initial_perturb
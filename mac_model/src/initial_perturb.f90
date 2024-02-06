module initial_perturb

    ! This module contains the subroutine to set initial perturbation from the base state

    implicit none

    contains

    subroutine init_perturb()
        use constants, only: g,cp
        use run_constants, only: nz, nx, radx,radz,amp,zcnt,xcnt
        use thermo_functions, only: calc_thv, calc_rsat, calc_satfrac
        use model_vars, only:zsn, xsn,thb, rvb, thvb, pib, piwb, rhoub,rhowb,thp,rvp,pip,up,wp,pp

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate
        real :: rad ! distance from center of perturbation

        do iz = 2, nz-1
            do ix = 2, nx-1
                rad = (((zsn(iz)-zcnt)/radz)**2 + ((xsn(ix)-xcnt)/radx)**2)**(0.5)
                
                ! construct temperature perturbation
                if (rad<=1) then
                    thp(ix,iz) = 0.5 * amp * (COS(rad*4*ATAN(1.)) + 1)
                endif
            enddo
        enddo
        
        !construct pressure perturbation field
        do ix = 2, nx - 1
            pip(ix,nz-1) = 0.

            do iz = nz-2,2,-1
                pip(ix,iz) = pip(ix,iz+1) - (g/cp)*(thp(ix,iz)/thb(iz)**2)*(zsn(iz-1)-zsn(iz))
            enddo 
        enddo 

        do iz = 2,nz-2
            pip(1,iz) = pip(2,iz)
            pip(nz,iz) = pip(nz-1,iz)
            thp(1,iz) = thp(2,iz)
            thp(nz,iz) = thp(nz-1,iz)
        enddo
        
        do iz = 1,nz
            do ix = 1,nx
                pp(ix,iz) = pip(ix,iz)*cp*rhoub(iz)*thvb(iz)
            enddo
        enddo

        print*,pp(:,2)

    end subroutine init_perturb

end module initial_perturb
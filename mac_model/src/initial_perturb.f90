module initial_perturb

    ! This module contains the subroutine to set initial perturbation from the base state

    implicit none

    contains

    subroutine init_perturb()
        use constants, only: g,cp,trigpi
        use run_constants, only: nz, nx, radx,radz,amp,zcnt,xcnt, pbc_x, pbc_z,pert_wind
        use model_vars, only:zsn, xsn,thb, thvb,rhoub,thp,pip,pp,up,wp

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
                pip(ix,iz,2) = pip(ix,iz+1,2) - (g/cp)*(thp(ix,iz,2)/thb(iz)**2)*(zsn(iz)-zsn(iz-1))
            enddo 
        enddo 

        ! set fictitous points to be equal to first real point
        if (pbc_x==.True.) then 
            do iz = 2,nz-1
                pip(1,iz,2) = pip(nx-1,iz,2)
                pip(nx,iz,2) = pip(2,iz,2)
                thp(1,iz,2) = thp(nx-1,iz,2)
                thp(nx,iz,2) = thp(2,iz,2)
                up(1,iz,2) = up(nx-1,iz,2)
                up(nx,iz,2) = up(2,iz,2)
                wp(1,iz,2) = wp(nx-1,iz,2)
                wp(nx,iz,2) = wp(2,iz,2)
            enddo ! end z loop
        else 
            ! if not using PBCs, just set 1st point to be same as 2nd, last pt to be same as 2nd to last, etc.
            do iz = 2,nz-1
                pip(1,iz,2) = pip(2,iz,2)
                pip(nx,iz,2) = pip(nx-1,iz,2)
                thp(1,iz,2) = thp(2,iz,2)
                thp(nx,iz,2) = thp(nx-1,iz,2)
                up(1,iz,2) = up(2,iz,2)
                up(nx,iz,2) = up(nx-1,iz,2)
                wp(1,iz,2) = wp(2,iz,2)
                wp(nx,iz,2) = wp(nx-1,iz,2)
            enddo
        endif !end x pbc settings

        ! set fictitous points to be equal to first real point
        if (pbc_z==.True.) then 
            do ix = 1,nx
                pip(ix,1,2) = pip(ix,nz-1,2)
                pip(ix,nz,2) = pip(ix, 2,2)
                thp(ix,1,2) = thp(ix,nz-1,2)
                thp(ix,nz,2) = thp(ix,2,2)
                up(ix,1,2) = up(ix,nz-1,2)
                up(ix,nz,2) = up(ix,2,2)
                wp(ix,1,2) = wp(ix,nz-1,2)
                wp(ix,nz,2) = wp(ix,2,2)
            enddo ! end x loop
        else 
            ! if not using PBCs, just set 1st point to be same as 2nd, last pt to be same as 2nd to last, etc.
            do ix = 1,nx
                pip(ix,1,2) = pip(ix,2,2)
                pip(ix,nz,2) = pip(ix, nz-1,2)
                thp(ix,1,2) = thp(ix,2,2)
                thp(ix,nz,2) = thp(ix,nz-1,2)
                up(ix,1,2) = up(ix,2,2)
                up(ix,nz,2) = up(ix,nz-1,2)

                wp(ix,2,2) = 0.
                wp(ix,1,2) = wp(ix,2,2)
                wp(ix,nz,2) = 0.
            enddo ! end x loop
        endif ! end z pbc

        ! for now, set the past to be the same as present timestep
        ! since we don't have any information about the past
        do iz = 1, nz
            do ix = 1, nx
                thp(ix,iz,1) = thp(ix,iz,2)
                pip(ix,iz,1) = pip(ix,iz,2)
                up(ix,iz,1) = up(ix,iz,2)
                wp(ix,iz,1) = wp(ix,iz,2)
            enddo
        enddo
        
        !calculate perturbation pressure in Pa
        do iz = 1,nz
            do ix = 1,nx
                pp(ix,iz) = pip(ix,iz,2)*cp*rhoub(iz)*thvb(iz)
            enddo
        enddo

    end subroutine init_perturb

end module initial_perturb
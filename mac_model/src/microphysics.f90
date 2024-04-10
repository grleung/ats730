module microphysics
    !This module contains subroutines for doing the microphysical calculations

    implicit none

    contains
    
    subroutine calc_autoconversion
        use run_constants, only: nz,nx,ny,cldautothresh,autorate
        use model_vars, only:cld2rain_auto,rcp

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        ! first, reset all tendencies to zero
        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    if (rcp(ix,iy,iz) >= cldautothresh) then
                        cld2rain_auto(ix,iy,iz) = autorate * (rcp(ix,iy,iz,2) - cldautothresh)
                    else
                        cld2rain_auto(ix,iy,iz) = 0.
                    endif
                enddo
            enddo
        enddo

    end subroutine calc_autoconversion

    subroutine calc_accretion
        use run_constants, only: nz,nx,ny,accrrate
        use model_vars, only:cld2rain_accr,rcp,rrp,rhoub

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        ! first, reset all tendencies to zero
        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    cld2rain_accr(ix,iy,iz) = accrrate * rhoub(iz) * rcp(ix,iy,iz,2) * rrp(ix,iy,iz,2)**(7/8)
                enddo
            enddo
        enddo

    end subroutine calc_accretion

    subroutine calc_rainevap
        use constants, only: p00, cp, rd
        use run_constants, only: nz,nx,ny,accrrate
        use model_vars, only:rain2vap,rrp,rvp,rhoub,thb,thp,pib,pip
        use micro_functions, only: rain_vent
        use thermo_functions, only: calc_rsat

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        real :: fvent,rsat,th,pi

        ! first, reset all tendencies to zero
        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    fvent = rain_vent(rhoub(iz),rrp(ix,iy,iz,2))
                    
                    th = thb(iz)+thp(ix,iy,iz,2)
                    pi = pib(iz)+pip(ix,iy,iz,2)
                    rsat = calc_rsat(pi*th, p00* (pi**(cp/rd)))

                    rain2vap(ix,iy,iz) = (1/rhoub(iz)) * ((1-(rvp(ix,iy,iz,2)/rsat))*fvent*(rhoub(iz)*rrp(ix,iy,iz,2))**.525) &
                                                        / (2.03e4 * (9.58e6/(rhoub(iz)*rsat)))
                    
                enddo
            enddo
        enddo

    end subroutine calc_rainevap


end module microphysics
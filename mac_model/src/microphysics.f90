module microphysics
    !This module contains subroutines for doing the microphysical calculations

    implicit none

    contains

    subroutine check_negs
        use run_constants, only: minmix,nx,ny,nz
        use model_vars, only: rvb,rvp,rcp,rrp

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        real :: vapadd=0,cldadd=0,rainadd=0

        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    if ((rvb(iz)+rvp(ix,iy,iz,3))<0.) then
                        !keep track of how much vapor we have removed so we can add it back in later
                        vapadd = vapadd + (rvp(ix,iy,iz,3)+rvb(iz))
                        rvp(ix,iy,iz,3) = -rvb(iz)
                    endif 

                    if (rcp(ix,iy,iz,3) < 0.) then
                        cldadd = cldadd + rcp(ix,iy,iz,3)
                        rcp(ix,iy,iz,3) = 0.
                    endif 

                    if (rrp(ix,iy,iz,3) < 0.) then
                        rainadd = rainadd + rrp(ix,iy,iz,3)
                        rrp(ix,iy,iz,3) = 0.
                    endif 
                enddo
            enddo
        enddo

        ! calculate amount of vapor to add to each grid point
        vapadd = -vapadd/(nx*ny*nz)
        cldadd = -cldadd/(nx*ny*nz)
        rainadd = -rainadd/(nx*ny*nz)

        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    rvp(ix,iy,iz,3) = rvp(ix,iy,iz,3)+vapadd
                    rcp(ix,iy,iz,3) = rcp(ix,iy,iz,3)+cldadd
                    rrp(ix,iy,iz,3) = rrp(ix,iy,iz,3)+rainadd
                enddo
            enddo
        enddo


    end subroutine check_negs

    subroutine add_tend_micro
        use model_vars, only: rvp_tend_total, rcp_tend_total, rrp_tend_total, thp_tend_total,pib &
                            , rain2vap,cld2rain_accr,cld2rain_auto,vap2cld, rcp, rrp,rvp,rvb
        use run_constants, only: nx,ny,nz,rd2t,rdt,dt
        use constants, only: cp, lv

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        real :: cldavail, cldsink, cldrat, rainavail,cldex,rainex

        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    ! not yet balanced

                    ! check that the amount being removed through accretion and autoconversion is less than the existing amount of cloud
                    
                    cldavail = (rcp(ix,iy,iz,1)*rd2t) !+ rcp_tend_total(ix,iy,iz)
                    cldsink = cld2rain_accr(ix,iy,iz) + cld2rain_auto(ix,iy,iz) 

                    if (cldsink > cldavail) then
                        cldrat = cld2rain_accr(ix,iy,iz)/cldsink

                        cld2rain_accr(ix,iy,iz) = cld2rain_accr(ix,iy,iz) - (cldrat * (cldsink - cldavail))
                        cld2rain_auto(ix,iy,iz) = cld2rain_auto(ix,iy,iz) - ((1-cldrat) * (cldsink - cldavail))
                    endif

                    rcp_tend_total(ix,iy,iz) = rcp_tend_total(ix,iy,iz) - cld2rain_accr(ix,iy,iz) - cld2rain_auto(ix,iy,iz) 

                    rainavail = rrp(ix,iy,iz,1)*rd2t !+ rrp_tend_total(ix,iy,iz)

                    if (rain2vap(ix,iy,iz) > rainavail) then
                        rain2vap(ix,iy,iz) = rainavail
                    endif

                    rrp_tend_total(ix,iy,iz) = rrp_tend_total(ix,iy,iz) + cld2rain_accr(ix,iy,iz) + cld2rain_auto(ix,iy,iz)  - rain2vap(ix,iy,iz)

                    thp_tend_total(ix,iy,iz) = thp_tend_total(ix,iy,iz) - ((lv/(cp*pib(iz))) * (rain2vap(ix,iy,iz)))
                    rvp_tend_total(ix,iy,iz) = rvp_tend_total(ix,iy,iz) + rain2vap(ix,iy,iz) 
                    
                enddo
            enddo
        enddo
    end subroutine add_tend_micro


    subroutine sat_adjust
        use model_vars, only: vap2cld, pib,thb,rvb,pip,thp,rvp,rcp          &
                            ,rvp_tend_total,rcp_tend_total,thp_tend_total
        use thermo_functions, only: calc_rsat
        use constants, only: p00, cp, rd,lv
        use run_constants, only: nx,ny,nz,rdt

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        real :: rsat,th,pi,phi,rv

        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx

                    th = thb(iz)+thp(ix,iy,iz,3)
                    pi = pib(iz)+pip(ix,iy,iz,3)
                    rv = rvb(iz)+rvp(ix,iy,iz,3)
                    rsat = calc_rsat(pi*th, p00* (pi**(cp/rd)))

                    if (rv > rsat) then
                        phi = rsat * (17.27 * 237. * lv)/(cp*((pi*th) - 36.)**2)

                        vap2cld(ix,iy,iz) = (rv - rsat)/(1+phi)

                        ! if microphysics wants to condense more water than is available
                        ! set condensed vapor to be maximum amount of vapor
                        if (vap2cld(ix,iy,iz) > rvb(iz)+rvp(ix,iy,iz,1)) then
                            vap2cld(ix,iy,iz) = rvb(iz)+rvp(ix,iy,iz,1)
                        else if (-vap2cld(ix,iy,iz) > rcp(ix,iy,iz,1)) then
                            vap2cld(ix,iy,iz) = -rcp(ix,iy,iz,1)
                        endif

                    else 
                        vap2cld(ix,iy,iz) = 0.
                    endif

                    rvp_tend_total(ix,iy,iz) = rvp_tend_total(ix,iy,iz) - rdt*vap2cld(ix,iy,iz)
                    rcp_tend_total(ix,iy,iz) = rcp_tend_total(ix,iy,iz) + rdt*vap2cld(ix,iy,iz)
                    thp_tend_total(ix,iy,iz) = thp_tend_total(ix,iy,iz) + ((lv/(cp*pib(iz))) * vap2cld(ix,iy,iz) * rdt)
                enddo
            enddo
        enddo

    end subroutine sat_adjust

    
    subroutine calc_autoconversion
        use run_constants, only: nz,nx,ny,cldautothresh,autorate
        use model_vars, only:cld2rain_auto,rcp

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        
        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    if (rcp(ix,iy,iz,2) > cldautothresh) then
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

        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    cld2rain_accr(ix,iy,iz) = accrrate * rhoub(iz) * rcp(ix,iy,iz,2) * rrp(ix,iy,iz,2)**(7/8)

                    if (cld2rain_accr(ix,iy,iz)<0.) then
                        cld2rain_accr(ix,iy,iz)=0.
                    endif

                enddo
            enddo
        enddo

    end subroutine calc_accretion

    subroutine calc_rainevap
        use constants, only: p00, cp, rd
        use run_constants, only: nz,nx,ny,accrrate
        use model_vars, only:rain2vap,rrp,rvp,rvb,rhoub,thb,thp,pib,pip
        use micro_functions, only: rain_vent
        use thermo_functions, only: calc_rsat

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        real :: fvent,rsat,th,pi,rv

        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    fvent = rain_vent(rhoub(iz),rrp(ix,iy,iz,2))

                    th = thb(iz)+thp(ix,iy,iz,2)
                    pi = pib(iz)+pip(ix,iy,iz,2)
                    rv = rvp(ix,iy,iz,2)+rvb(iz)
                    rsat = calc_rsat(pi*th, p00*(pi**(cp/rd)))

                    !condition here in case it's supersaturated so we aren't condensing onto raindrops
                    if (rv < rsat) then
                        rain2vap(ix,iy,iz) = (1/rhoub(iz)) * ((1-(rv/rsat))*fvent*(rhoub(iz)*rrp(ix,iy,iz,2))**.525) &
                                                        / (2.03e4 * (9.58e6/(rhoub(iz)*rsat)))
                    else 
                        rain2vap(ix,iy,iz) = 0.
                    endif

                enddo
            enddo
        enddo

    end subroutine calc_rainevap
    

end module microphysics
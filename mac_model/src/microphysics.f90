module microphysics
    !This module contains subroutines for doing the microphysical calculations

    implicit none

    contains

    subroutine check_negs
        use run_constants, only: minmix,nx,ny,nz,mincld,minrain
        use model_vars, only: rvb,rvp,rcp,rrp

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate


        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    if ((rvb(iz)+rvp(ix,iy,iz,3))<0.) then
                        rvp(ix,iy,iz,3) = -rvb(iz)
                    endif 

                    if (rcp(ix,iy,iz,3) < 0.) then
                        rcp(ix,iy,iz,3) = 0.
                    endif 

                    if (rrp(ix,iy,iz,3) < 0.) then
                        rrp(ix,iy,iz,3) = 0.
                    endif 
                enddo
            enddo
        enddo
    end subroutine check_negs

    subroutine sat_adjust
        use model_vars, only: vap2cld, pib,thb,rvb,pip,thp,rvp,rcp         
        use thermo_functions, only: calc_rsat
        use constants, only: p00, cp, rd,lv
        use run_constants, only: nx,ny,nz,rd2t

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        real :: rsat,th,pi,phi,rv,t,p

        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    th = thp(ix,iy,iz,3) + thb(iz)
                    pi = pip(ix,iy,iz,3) + pib(iz)
                    t = th*pi
                    p = p00 * (pi**(cp/rd))
                    rv = rvp(ix,iy,iz,3) + rvb(iz)

                    rsat = calc_rsat(t, p)

                    phi = rsat * (17.27 * 237. * lv)/(cp*(t-36.)**2)

                    vap2cld(ix,iy,iz) = (rv - rsat)/(1+phi)

                    ! if microphysics wants to condense more vapor than is available
                    ! set condensed vapor to be maximum amount of vapor
                    if (vap2cld(ix,iy,iz) > rv) then
                        vap2cld(ix,iy,iz) = rv
                    ! if it wants to evaporate more cloud than is available
                    ! set negative condensed water to be maximum amount of cloud 
                    else if (-vap2cld(ix,iy,iz) > rcp(ix,iy,iz,3)) then
                        vap2cld(ix,iy,iz) = -rcp(ix,iy,iz,3)
                    endif

                    thp(ix,iy,iz,3) = thp(ix,iy,iz,3) + (lv/(cp*pib(iz))) * vap2cld(ix,iy,iz)
                    rcp(ix,iy,iz,3) = rcp(ix,iy,iz,3) + vap2cld(ix,iy,iz)
                    rvp(ix,iy,iz,3) = rvp(ix,iy,iz,3) - vap2cld(ix,iy,iz)
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
        use run_constants, only: nz,nx,ny,accrrate,mincld,minrain
        use model_vars, only:cld2rain_accr,rcp,rrp,rhoub

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    if (rcp(ix,iy,iz,2)>0.) then
                        cld2rain_accr(ix,iy,iz) = accrrate * rhoub(iz) * rcp(ix,iy,iz,2) * (rrp(ix,iy,iz,2)**(7./8.))

                        if (cld2rain_accr(ix,iy,iz)<0.) then
                            cld2rain_accr(ix,iy,iz)=0.
                        endif
                    else 
                        cld2rain_accr(ix,iy,iz)=0.
                        rcp(ix,iy,iz,2) = 0.
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

        real :: fvent,rsat,th,pi,rv,pb,t,p

        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    fvent = rain_vent(rhoub(iz),rrp(ix,iy,iz,2))

                    th = thp(ix,iy,iz,2)+thb(iz)
                    pi = pip(ix,iy,iz,2)+pib(iz)
                    t = th*pi
                    p = p00 * (pi**(cp/rd))
                    pb = p00*pib(iz) **(cp/rd)
                    rv = rvp(ix,iy,iz,2)+rvb(iz)
                    
                    rsat = calc_rsat(t, p)

                    !condition here in case it's supersaturated so we aren't condensing onto raindrops
                    if (rv < rsat) then
                        rain2vap(ix,iy,iz) = (1/rhoub(iz)) * ((1-(rv/rsat))*fvent*(rhoub(iz)*rrp(ix,iy,iz,2))**0.525) &
                                                        / (2.03e4 + (9.58e6/(pb*rsat)))
                    else 
                        rain2vap(ix,iy,iz) = 0.
                    endif
                
                enddo
            enddo
        enddo

    end subroutine calc_rainevap
    
    subroutine apply_micro_tends
        use model_vars, only: rcp_tend_total, rrp_tend_total,pib &
                            , rain2vap,cld2rain_accr,cld2rain_auto,vap2cld, rcp, rrp,rvp,rvb,thp
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
                    ! check that the amount being removed through accretion and autoconversion is less than the existing amount of cloud
                    !print*,rcp_tend_total(ix,iy,iz),rrp_tend_total(ix,iy,iz)
                    
                    cldavail = (rcp(ix,iy,iz,1)*rd2t) + rcp_tend_total(ix,iy,iz)
                    cldsink = cld2rain_accr(ix,iy,iz) + cld2rain_auto(ix,iy,iz) 
                    
                    if (cldavail>0.) then
                        if (cldsink > cldavail) then
                            cldrat = cld2rain_accr(ix,iy,iz)/cldsink

                            cld2rain_accr(ix,iy,iz) = cldrat * cldavail
                            cld2rain_auto(ix,iy,iz) = (1-cldrat) * cldavail
                        endif
                    else
                        cld2rain_accr(ix,iy,iz) = 0.
                        cld2rain_auto(ix,iy,iz) = 0.
                    endif

                    if (cld2rain_accr(ix,iy,iz)<0.) then
                        cld2rain_accr(ix,iy,iz) = 0.
                    endif

                    if (cld2rain_auto(ix,iy,iz)<0.) then
                        cld2rain_auto(ix,iy,iz) = 0.
                    endif

                    rcp(ix,iy,iz,3) = rcp(ix,iy,iz,3) - dt*(cld2rain_accr(ix,iy,iz)+cld2rain_auto(ix,iy,iz))

                    rainavail = rrp(ix,iy,iz,1)*rd2t + rrp_tend_total(ix,iy,iz)

                    if (rainavail>0) then 
                        if (rain2vap(ix,iy,iz) > rainavail) then
                            rain2vap(ix,iy,iz) = rainavail
                        endif
                    else
                        rain2vap(ix,iy,iz) = 0.
                    endif
                    
                    !more checks on negative values -- not sure why this should be needed but if not some of the rain2vap becomes negative sometimes? maybe floating pt math
                    if (rain2vap(ix,iy,iz)<0.) then
                        rain2vap(ix,iy,iz) = 0.
                    endif

                    rrp(ix,iy,iz,3) = rrp(ix,iy,iz,3) +  dt*(cld2rain_accr(ix,iy,iz)+cld2rain_auto(ix,iy,iz)-rain2vap(ix,iy,iz))
                    thp(ix,iy,iz,3) = thp(ix,iy,iz,3) -  dt*(lv/(cp*pib(iz))) * (rain2vap(ix,iy,iz))
                    rvp(ix,iy,iz,3) = rvp(ix,iy,iz,3) +  dt*rain2vap(ix,iy,iz)
                    
                enddo
            enddo
        enddo
    end subroutine apply_micro_tends



end module microphysics
module microphysics

    implicit none
    
    contains

    subroutine check_micro_zeros
        use model_vars, only: rvb,rvp, rcp, rrp
        use run_constants, only: nx, nz
        
        implicit none

        integer :: ix, iz
        do iz=1,nz
            do ix=1,nx
                if ((rvp(ix,iz,3)+rvb(iz))<0.) then
                    rvp(ix,iz,3) = -rvb(iz)
                endif 

                if (rcp(ix,iz,3)<0.) then
                    rcp(ix,iz,3) = 0.
                endif

                if (rrp(ix,iz,3)<0.) then
                    rrp(ix,iz,3) = 0.
                endif
            enddo
        enddo
        
    end subroutine check_micro_zeros

    subroutine sat_adjust
        use model_vars, only: vap2cld, rvp, rcp, pip, thp             &   
                        , rvb, thb, pib
        use run_constants, only: nx, nz
        use thermo_functions, only: calc_rsat
        use constants, only: lv, cp, p00, rd

        implicit none

        integer :: ix,iz
        real    :: phi,rsat,th,pi,t,p,rv

        do iz=1,nz
            do ix=1,nx
                th = thp(ix,iz,3) + thb(iz)
                pi = pip(ix,iz,3) + pib(iz)
                t = th*pi
                p = p00 * pi**(cp/rd)
                rv = rvp(ix,iz,3) + rvb(iz)

                rsat = calc_rsat(t, p)

                phi = rsat * (17.27*237.*lv)/(cp*(t-36.)**2)
                vap2cld(ix,iz) = (rv - rsat)/(1+phi)

                if (vap2cld(ix,iz) > rv) then
                    vap2cld(ix,iz) = rv
                else if (-vap2cld(ix,iz) > rcp(ix,iz,3)) then 
                    vap2cld(ix,iz) = -rcp(ix,iz,3)
                endif

                thp(ix,iz,3) = thp(ix,iz,3) + (lv/(cp*pib(iz)))*vap2cld(ix,iz)
                rvp(ix,iz,3) = rvp(ix,iz,3) - vap2cld(ix,iz)
                rcp(ix,iz,3) = rcp(ix,iz,3) + vap2cld(ix,iz)
            enddo
        enddo
    end subroutine sat_adjust

    subroutine autoconversion
        use model_vars, only: cld2rain_auto, rcp
        use run_constants, only: nx, nz,cldautothresh,autorate

        implicit none

        integer :: ix,iz

        do iz=1,nz
            do ix=1,nx
                if (rcp(ix,iz,2) > cldautothresh) then
                    cld2rain_auto(ix,iz) = autorate * (rcp(ix,iz,2) - cldautothresh)
                else
                    cld2rain_auto(ix,iz) = 0.
                endif
            enddo
        enddo
    end subroutine autoconversion

    subroutine accretion
        use run_constants, only: nz,nx,accrrate,mincld,minrain
        use model_vars, only:cld2rain_accr,rcp,rrp,rhoub

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate

        do iz = 1, nz
            do ix = 1, nx
                if (rcp(ix,iz,2)>0.) then
                    cld2rain_accr(ix,iz) = accrrate * rhoub(iz) * rcp(ix,iz,2) * (rrp(ix,iz,2)**(7./8.))

                    if (cld2rain_accr(ix,iz)<0.) then
                        cld2rain_accr(ix,iz)=0.
                    endif
                else 
                    cld2rain_accr(ix,iz)=0.
                    rcp(ix,iz,2) = 0.
                endif 

            enddo
        enddo
        
    end subroutine accretion


    subroutine rainevap
        use constants, only: p00, cp, rd
        use run_constants, only: nz,nx,accrrate
        use model_vars, only:rain2vap,rrp,rvp,rvb,rhoub,thb,thp,pib,pip
        use thermo_functions, only: calc_rsat,rain_vent

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate

        real :: fvent,rsat,th,pi,rv,pb,t,p

        do iz = 1, nz
            do ix = 1, nx
                fvent = rain_vent(rhoub(iz),rrp(ix,iz,2))

                th = thp(ix,iz,2) + thb(iz)
                pi = pip(ix,iz,2) + pib(iz)
                t = th*pi
                p = p00 * pi**(cp/rd)
                pb = p00*pib(iz) **(cp/rd)
                rv = rvp(ix,iz,2) + rvb(iz)

                rsat = calc_rsat(t, p)

                !condition here in case it's supersaturated so we aren't condensing onto raindrops
                if (rv < rsat) then
                    rain2vap(ix,iz) = (1/rhoub(iz)) * ((1-(rv/rsat))*fvent*(rhoub(iz)*rrp(ix,iz,2))**0.525) &
                                                    / (2.03e4 + (9.58e6/(pb*rsat)))
                else 
                    rain2vap(ix,iz) = 0.
                endif

            enddo
        enddo
    end subroutine rainevap

    subroutine apply_micro_tends
        use model_vars, only: cld2rain_auto, cld2rain_accr, rain2vap, rcp, rrp, rvp, thp &
                        ,pib,rvp_tend_total,rcp_tend_total,rrp_tend_total,thp_tend_total
        use run_constants, only: nx, nz,rd2t, d2t,dt
        use constants, only: lv, cp

        implicit none

        integer :: ix,iz

        real :: cldavail, cldsink, cldrat, rainavail,cldex,rainex

        do iz = 1, nz
            do ix = 1, nx
                ! not yet balanced

                ! check that the amount being removed through accretion and autoconversion is less than the existing amount of cloud
                !print*,rcp_tend_total(ix,iz),rrp_tend_total(ix,iz)
                
                cldavail = (rcp(ix,iz,1)*rd2t) + rcp_tend_total(ix,iz)
                cldsink = cld2rain_accr(ix,iz) + cld2rain_auto(ix,iz) 
                
                if (cldavail>0) then
                    if (cldsink > cldavail) then
                        cldrat = cld2rain_accr(ix,iz)/cldsink

                        cld2rain_accr(ix,iz) = cldrat * cldavail
                        cld2rain_auto(ix,iz) = (1-cldrat) * cldavail
                    endif
                else
                    cld2rain_accr(ix,iz) = 0.
                    cld2rain_auto(ix,iz) = 0.
                endif

                if (cld2rain_accr(ix,iz)<0) then
                    cld2rain_accr(ix,iz) = 0.
                endif

                if (cld2rain_auto(ix,iz)<0) then
                    cld2rain_auto(ix,iz) = 0.
                endif

                rcp(ix,iz,3) = rcp(ix,iz,3) - dt*(cld2rain_accr(ix,iz) + cld2rain_auto(ix,iz) )

                rainavail = rrp(ix,iz,1)*rd2t + rrp_tend_total(ix,iz)

                if (rainavail>0) then 
                    if (rain2vap(ix,iz) > rainavail) then
                        rain2vap(ix,iz) = rainavail
                    endif
                else
                    rain2vap(ix,iz) = 0.
                endif
                
                !more checks on negative values -- not sure why this should be needed but if not some of the rain2vap becomes negative sometimes? maybe floating pt math
                if (rain2vap(ix,iz)<0) then
                    rain2vap(ix,iz) = 0.
                endif

                rrp(ix,iz,3) = rrp(ix,iz,3) +  dt*(cld2rain_accr(ix,iz) + cld2rain_auto(ix,iz) -rain2vap(ix,iz))
                thp(ix,iz,3) = thp(ix,iz,3) -  dt*(lv/(cp*pib(iz)))*rain2vap(ix,iz)
                rvp(ix,iz,3) = rvp(ix,iz,3) + dt*rain2vap(ix,iz)
            enddo
        enddo


    end subroutine apply_micro_tends

end module microphysics
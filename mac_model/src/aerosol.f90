module aerosol
    implicit none

    contains

    real function calc_mass(dp, rho)
        use constants, only: trigpi

        implicit none
        real :: rho,dp

        calc_mass = rho * (4./3.) * trigpi * (dp/2.)**3.
    end function calc_mass

    real function calc_dp(m, rho)
        use constants, only: trigpi

        implicit none
        real :: rho,m

        calc_dp = 2 *((m/(rho*(4./3.)*trigpi)))**(1/3.)
    end function calc_dp
 

    subroutine init_aerosol
        use model_vars, only: np, mp
        use run_constants, only: nx, nz, npb, ndb
        use constants, only: trigpi

        ! very basic aerosol intialization for now
        ! constant everywhere and just two lognormal modes

        real :: kappa=0.5       !hygroscopicity parameter, assume 0.5
        real :: rhop = 1400.    !density of particle, assume ammonium sulfate
        real :: dpi=1.e-9       ! smallest particle diameter (1nm)
        real :: ntot = 1000.       ! total aerosol number conc [units of #/kg]
        real :: dpg = 100.e-9
        real :: sigma = 1.8

        ! define mass and diameters of bins
        ! only going to use this in this function for set-up
        real :: mp_c(npb), dp_c(npb)
        real :: mp_e(npb+1),dp_e(npb+1)

        real :: ratio = 2. ! mass bin ratio, by default use mass doubling bins
        integer :: ipb ! counter for aerosol bins 
        real :: a = 0.

        ! set up mass doubling bins, define both bin centers and edges
        dp_e(1) = dpi
        mp_e(1) = calc_mass(dpi,rhop)
        mp_c(1) = mp_e(1) * ratio**0.5
        dp_c(1) = calc_dp(mp_c(1),rhop)

        do ipb=2,npb
            mp_e(ipb) = mp_e(ipb-1)*ratio
            dp_e(ipb) = calc_dp(mp_e(ipb),rhop)
            
            mp_c(ipb) = mp_c(ipb-1) * ratio
            dp_c(ipb) = calc_dp(mp_c(ipb),rhop)
        enddo

        mp_e(npb+1) = mp_e(npb)*ratio
        dp_e(npb+1) = calc_dp(mp_e(npb+1),rhop)

        print*,'bins ok'

        ! set up lognormal distribution
        do ipb=1,npb
            np(:,:,ipb,:) = ((dp_c(ipb) * 2.303 * ntot/((2*trigpi)**0.5 * LOG(sigma) * dp_c(ipb)) &
                        * EXP(-(LOG(dp_c(ipb))-LOG(dpg))**2/(2*LOG(sigma)**2)))) * (LOG10(dp_e(ipb+1))-LOG10(dp_e(ipb)))
            mp(:,:,ipb,:) = mp_c(ipb) * np(:,:,ipb,:)
        enddo           

        print*,'dist ok'
    end subroutine init_aerosol


    subroutine microphysics_driv
        !this is the main microphysics driver
        use constants, only: cp, rd, p00, lv
        use run_constants, only: nx,nz,npb,ndb
        use model_vars, only: thp, pip, rvp, thb, pib, rvb, np, nc
        use thermo_functions, only: calc_rsat

        implicit none

        integer :: ix,iz,ipb,idb ! counters
        real :: th, pi, rv, t, p, rvsat ! temporary variables for saturation calc

        !first loop over all spatial points
        do iz=1,nz
            do ix=1,nx
                th = thp(ix,iz,2) + thb(iz)
                pi = pip(ix,iz,2) + pib(iz)
                t = th*pi
                p = p00 * (pi**(cp/rd))
                rv = rvp(ix,iz,2) + rvb(iz)
                rvsat = calc_rsat(t,p)

                ! if the grid point is supersaturated and
                ! there are fewer cloud droplets than available CCN locally
                if (rv>rvsat) then
                    ! calculate local CCN concentration
                    print*,SUM(nc(ix,iz,:,2)),SUM(np(ix,iz,:,2))
                    !if (nc(ix,iz,:,2)>np(ix,iz,:,2)) then
                        ! calculate local total cloud drop concentration

                    !endif
                endif 
            enddo
        enddo


        !calculate saturation
        !if saturated & more than past supersaturation (think about how to do this criteria), do activation
        !if there is cloud, do condensation onto clouds
        !account for temp, rv change due to condensation
        !if there is cloud, do coagulation 
        !do rainout -- if no time, could scrap this

    end subroutine microphysics_driv

end module aerosol
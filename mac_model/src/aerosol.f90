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
        integer :: iab ! counter for aerosol bins 

        ! set up mass doubling bins, define both bin centers and edges
        dp_e(1) = dpi
        mp_e(1) = calc_mass(dpi,rhop)
        mp_c(1) = mp_e(1) * ratio**0.5
        dp_c(1) = calc_dp(mp_c(1),rhop)

        do iab=2,npb
            mp_e(iab) = mp_e(iab-1)*ratio
            dp_e(iab) = calc_dp(mp_e(iab),rhop)
            
            mp_c(iab) = mp_c(iab-1) * ratio
            dp_c(iab) = calc_dp(mp_c(iab),rhop)
        enddo

        mp_e(npb+1) = mp_e(npb)*ratio
        dp_e(npb+1) = calc_dp(mp_e(npb+1),rhop)

        ! set up lognormal distribution
        do iab=1,npb
            np(:,:,iab,:) = ((dp_c(iab) * 2.303 * ntot/((2*trigpi)**0.5 * LOG(sigma) * dp_c(iab)) &
                        * EXP(-(LOG(dp_c(iab))-LOG(dpg))**2/(2*LOG(sigma)**2)))) &
                     * (LOG10(dp_e(iab+1))-LOG10(dp_e(iab)))
            mp(:,:,iab,:) = mp_c(iab) * np(:,:,iab,:)
        enddo           
    end subroutine init_aerosol

    subroutine microphysics
        !this is the main microphysics driver


        !calculate saturation
        !if saturated & more than past supersaturation (think about how to do this criteria), do activation
        !if there is cloud, do condensation onto clouds
        !account for temp, rv change due to condensation
        !if there is cloud, do coagulation 
        !do rainout -- if no time, could scrap this

    end subroutine microphysics

end module aerosol
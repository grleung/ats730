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

        calc_dp = 2. *((m/(rho*(4./3.)*trigpi)))**(1./3.)
    end function calc_dp
 
    real function calc_scrit(dp, T, kappa)
        use constants, only: MW_w, sigma_w,rhol, rd

        implicit none

        real :: dp, T, kappa
        real :: A, B

        A = 4*MW_w*sigma_w/(rd*T*rhol)
        B = kappa * dp**3
        calc_scrit=EXP(((4*A**3)/(27*B))**0.5)

    end function calc_scrit


    real function calc_lambd(p,t)
        use constants, only: trigpi, R, MW_air,mu
        implicit none

        real:: p, t

        calc_lambd = (2*mu)/(p*(8*MW_air)/(trigpi*R*t))
    end function calc_lambd

    real function calc_dahneke(dp, lambd,alpha)
        implicit none

        real:: dp, lambd, alpha, kn

        !knudsen number
        kn = 2*lambd/dp
        calc_dahneke = (1+kn)/(1+(2*kn*(1+kn)/alpha))
    end function calc_dahneke

    subroutine init_aerosol
        use model_vars, only: np, mp, mdropbin_lims, mpartbin_lims, nc, mlc, mpc 
        use run_constants, only: nx, nz, npartbin, ndropbin   &
                            , kappa,rhop,ntot,dpg,sigma
        use constants, only: trigpi,rhol

        ! very basic aerosol intialization for now
        ! constant everywhere and just two lognormal modes

        real :: dpi=1.e-9       ! smallest particle diameter (1nm)
        
        ! define mass and diameters of bins
        ! only going to use this in this function for set-up
        real :: mp_c(npartbin), dp_c(npartbin), dp_e(npartbin+1)

        real :: ratio = 2. ! mass bin ratio, by default use mass doubling bins
        integer :: ipb,idb,ix,iz,it ! counter for aerosol bins 
        real :: a = 0.

        ! set up mass doubling bins, define both bin centers (_c) and edges (_e)
        dp_e(1) = dpi
        mpartbin_lims(1) = calc_mass(dpi,rhop)
        mp_c(1) = mpartbin_lims(1) * ratio**0.5
        dp_c(1) = calc_dp(mp_c(1),rhop)

        do ipb=2,npartbin
            mpartbin_lims(ipb) = mpartbin_lims(ipb-1)*ratio
            dp_e(ipb) = calc_dp(mpartbin_lims(ipb),rhop)
            
            mp_c(ipb) = mp_c(ipb-1) * ratio
            dp_c(ipb) = calc_dp(mp_c(ipb),rhop)
        enddo

        ! define extra point for largest bin edge
        mpartbin_lims(npartbin+1) = mpartbin_lims(npartbin)*ratio
        dp_e(npartbin+1) = calc_dp(mpartbin_lims(npartbin+1),rhop)

        print*,'bins ok'

        ! set up lognormal distribution
        ! for now this is based on unimodal distribution with a peak at dpg (hard coded), but should be easy to change this to be Namelist param
        do it =1,3
            do ix=2,nx-1
                do iz=2,nz-1
                    do ipb=1,npartbin
                        np(ix,iz,ipb,it) = ((dp_c(ipb) * 2.303 * ntot/((2*trigpi)**0.5 * LOG(sigma) * dp_c(ipb)) &
                                    * EXP(-(LOG(dp_c(ipb))-LOG(dpg))**2/(2*LOG(sigma)**2)))) * (LOG10(dp_e(ipb+1))-LOG10(dp_e(ipb)))
                        mp(ix,iz,ipb,it) = mp_c(ipb) * np(ix,iz,ipb,it)

                        do idb=1,ndropbin
                            nc(ix,iz,ipb,idb,it) = 0. 
                            mlc(ix,iz,ipb,idb,it) = 0. 
                            mpc(ix,iz,ipb,idb,it) = 0. 
                        enddo
                    enddo
                enddo
            enddo
        enddo           

        print*,'dist ok'

        ! set up droplet bins
        mdropbin_lims(1) = calc_mass(1.e-8,rhol) ! first bin has radius of 10nm
        print*,calc_dp(mdropbin_lims(1),rhol)
        do idb=2,ndropbin+1
            mdropbin_lims(idb) = mdropbin_lims(idb-1) * ratio !mass doubling bins
        enddo

        print*,'drops',ndropbin,calc_dp(mdropbin_lims(1),rhol),calc_dp(mdropbin_lims(ndropbin+1),rhol)
        print*,'particles',npartbin,calc_dp(mpartbin_lims(1),rhop),calc_dp(mpartbin_lims(npartbin+1),rhop)

    end subroutine init_aerosol


    

end module aerosol
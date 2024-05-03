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

        calc_lambd = (2*mu)/(p*((8*MW_air)/(trigpi*R*t))**0.5)
    end function calc_lambd

    real function calc_dahneke(dp, lambd,alpha)
        implicit none

        real:: dp, lambd, alpha, kn

        !knudsen number
        kn = 2*lambd/dp
        calc_dahneke = (1+kn)/(1+(2*kn*(1+kn)/alpha))
    end function calc_dahneke

    real function calc_stokes(dp, p, t)
        use constants, only: rhol,g,mu

        implicit none

        real :: dp, p, t, lambd,cc

        lambd = calc_lambd(p,t)
        ! stokes estimate of terminal velocity
        cc = 1. + (2.*lambd/dp)*(1.257+(0.4*EXP((-1.1*dp)/(2.*lambd))))
        calc_stokes = (cc*rhol*g*(dp**2))/(18*mu)

    end function calc_stokes

    real function calc_drag(dp,v,p,t)
        use constants, only: mu, g, rhol
        use thermo_functions, only: calc_rhoair

        implicit none

        real :: dp, v, p, t, rho_air, re, cc, cd,lambd

        lambd = calc_lambd(p,t)

        ! slip correction factor 
        cc =  1. + (2.*lambd/dp)*(1.257+(0.4*EXP((-1.1*dp)/(2.*lambd))))
        rho_air = calc_rhoair(p,t)

        !reynolds number
        re = rho_air*v*dp/mu

        ! drag coefficient
        cd = (24./re) * (1+((3./16.)*re*0.43))

        ! set drag to Fg=Fdrag
        calc_drag = ((4*rhol*dp*g*cc)/(3*rho_air*cd))**0.5
    end function calc_drag

    real function calc_vt(dp,p,t)
        implicit none

        real :: uo, un,dp,p,t
        integer :: i=0
        logical :: repeat=.True.

        uo = calc_stokes(dp,p,t)

        do while (repeat)
            un = calc_drag(dp,uo,p,t)

            ! do until values converge

            if (ABS(un-uo)<=5.e-2) then
                repeat = .False.
            else if (i>=20) then
                repeat = .False.
            endif

            uo = un
            i = i+1
        enddo
    end function calc_vt

    subroutine init_aerosol
        use model_vars, only: np, mp, mdropbin_lims, mpartbin_lims, nd, mld, mpd 
        use run_constants, only: nx, nz, npartbin, ndropbin,ratio   &
                            , kappa,rhop,ntot,dpg,sigma
        use constants, only: trigpi,rhol

        ! very basic aerosol intialization for now
        ! constant everywhere and just two lognormal modes

        real :: dpi=5.e-9       ! smallest particle diameter (5nm)
        
        ! define mass and diameters of bins
        ! only going to use this in this function for set-up
        real :: mp_c(npartbin), dp_c(npartbin), dp_e(npartbin+1)

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
                            nd(ix,iz,ipb,idb,it) = 0. 
                            mld(ix,iz,ipb,idb,it) = 0. 
                            mpd(ix,iz,ipb,idb,it) = 0. 
                        enddo
                    enddo
                enddo
            enddo
        enddo           

        print*,'dist ok'

        print*,'total particles',SUM(np(2,2,:,2))

        ! set up droplet bins
        mdropbin_lims(1) = calc_mass(1.e-7,rhol) ! first bin has radius of 100nm
        print*,calc_dp(mdropbin_lims(1),rhol)
        do idb=2,ndropbin+1
            mdropbin_lims(idb) = mdropbin_lims(idb-1) * ratio !mass doubling bins
        enddo

        print*,'drops (mm)',ndropbin,calc_dp(mdropbin_lims(1),rhol)*1e3,calc_dp(mdropbin_lims(ndropbin+1),rhol)*1e3
        print*,'particles (micron)',npartbin,calc_dp(mpartbin_lims(1),rhop)*1e6,calc_dp(mpartbin_lims(npartbin+1),rhop)*1e6

    end subroutine init_aerosol


    

end module aerosol
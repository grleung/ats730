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
 
    real function calc_scrit(dp, T, kappa)
        use constants, only: MW_w, sigma_w,rhol, rd

        implicit none

        real :: dp, T, kappa
        real :: A, B

        A = 4*MW_w*sigma_w/(rd*T*rhol)
        B = kappa * dp**3
        calc_scrit=EXP(((4*A**3)/(27*B))**0.5)

    end function calc_scrit

    real function calc_esat(T)
        implicit none

        real :: T

        calc_esat = 100.*EXP(53.68-(6743./T) - (4.845*LOG(T)))
    end function calc_esat

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
        use model_vars, only: np, mp, mdb, nc, mc, mpc, nr, mr, mpr
        use run_constants, only: nx, nz, npb, ndb
        use constants, only: trigpi,rhol

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
        integer :: ipb,idb,ix,iz,it ! counter for aerosol bins 
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
        do it =1,3
            do ix=2,nx-1
                do iz=2,nz-1
                    do ipb=1,npb
                        np(ix,iz,ipb,it) = ((dp_c(ipb) * 2.303 * ntot/((2*trigpi)**0.5 * LOG(sigma) * dp_c(ipb)) &
                                    * EXP(-(LOG(dp_c(ipb))-LOG(dpg))**2/(2*LOG(sigma)**2)))) * (LOG10(dp_e(ipb+1))-LOG10(dp_e(ipb)))
                        mp(ix,iz,ipb,it) = mp_c(ipb) * np(ix,iz,ipb,it)

                        do idb=1,ndb
                            nc(ix,iz,ipb,idb,it) = 0. 
                            mc(ix,iz,ipb,idb,it) = 0. 
                            mpc(ix,iz,ipb,idb,it) = 0. 
                            nr(ix,iz,ipb,idb,it) = 0. 
                            mr(ix,iz,ipb,idb,it) = 0. 
                            mpr(ix,iz,ipb,idb,it) = 0. 
                        enddo
                    enddo
                enddo
            enddo
        enddo           

        print*,'dist ok'

        ! set up droplet bins
        mdb(1) = calc_mass(1.e-7,rhol) ! first bin has radius of 1micron
        print*,calc_dp(mdb(1),rhol)
        do idb=2,ndb
            mdb(idb) = mdb(idb-1) * ratio !mass doubling bins
        enddo
    

    end subroutine init_aerosol


    subroutine microphysics_driv
        !this is the main microphysics driver
        implicit none
    
        call check_negs
        print*,'check negs'

        call activation
        print*,'activation'

        call condensation
        print*,'condensation'

        call check_negs
        print*,'check negs'

        !calculate saturation
        !if saturated & more than past supersaturation (think about how to do this criteria), do activation
        !if there is cloud, do condensation onto clouds
        !account for temp, rv change due to condensation
        !if there is cloud, do coagulation 
        !do rainout -- if no time, could scrap this

    end subroutine microphysics_driv

    subroutine activation
        use constants, only: cp, rd, p00, lv
        use run_constants, only: nx,nz,npb,ndb
        use model_vars, only: mdb, thp, pip, rvp, thb, pib, rvb, np, mp,nc,mc,mpc
        use thermo_functions, only: calc_rsat

        implicit none

        real :: kappa=0.5       !hygroscopicity parameter, assume 0.5
        real :: rhop = 1400.    !density of particle, assume ammonium sulfate
        
        integer :: ix,iz,ipb,idb ! counters
        real :: th, pi, rv, t, p, rvsat, Samb ! temporary variables for saturation calc
        real :: ma,dp,Sc ! temporary variables for activation calc

        !first loop over all spatial points
        do iz=2,nz-1
            do ix=2,nx-1
                th = thp(ix,iz,3) + thb(iz)
                pi = pip(ix,iz,3) + pib(iz)
                t = th*pi
                p = p00 * (pi**(cp/rd))
                rv = rvp(ix,iz,3) + rvb(iz)
                rvsat = calc_rsat(t,p)
                Samb = rv/rvsat

                ! if the grid point is supersaturated and
                ! there are fewer cloud droplets than available CCN locally
                if (Samb>1.) then
                    ! calculate local CCN concentration
                    if (SUM(nc(ix,iz,:,:,3)<SUM(np(ix,iz,:,3)))) then
                        ! calculate local total cloud drop concentration
                        do ipb=1,npb
                            if ((np(ix,iz,ipb,3) > 0.) .and. (mp(ix,iz,ipb,3) > 0.)) then
                                ma = mp(ix,iz,ipb,3)/np(ix,iz,ipb,3) ! mean mass of aerosol particle in bin
                                dp = calc_dp(ma, rhop)
                                Sc = calc_scrit(dp, t, kappa)

                                if (Samb>Sc) then !if ambinet supersaturation is above critical, do activation
                                    ! determine which cloud bin the activate droplets should go into
                                    if ((5*ma) < mdb(1)) then
                                        idb = 1
                                    else
                                        idb = COUNT(mdb(:)<(5*ma))
                                    endif
                                    
                                    ! add cloud droplets & mass & aerosol mass
                                    nc(ix,iz,ipb,idb,3) = nc(ix,iz,ipb,idb,3) + np(ix,iz,ipb,3) 
                                    mc(ix,iz,ipb,idb,3) = mc(ix,iz,ipb,idb,3) + (mdb(idb)-ma)
                                    mpc(ix,iz,ipb,idb,3) = mpc(ix,iz,ipb,idb,3) + (np(ix,iz,ipb,3)*ma)

                                    ! add the equivalent amount of latent heating to THP
                                    thp(ix,iz,3) = thp(ix,iz,3) + (mdb(idb)-ma)*(lv/(cp*pib(iz)))
                                    ! remove equivalent amount of water from RVP
                                    rvp(ix,iz,3) = rvp(ix,iz,3) - (mdb(idb)-ma)

                                    ! remove aerosol from unprocessed aerosol bins
                                    np(ix,iz,ipb,3) = 0. 
                                    mp(ix,iz,ipb,3) = 0.
                                endif 
                            endif
                        enddo
                    endif
                endif 
            enddo
        enddo
    end subroutine activation

    subroutine condensation
        use constants, only: cp, rd, p00, lv,rhol,MW_w,sigma_w,trigpi,dg,kair,R
        use run_constants, only: nx,nz,npb,ndb,dt,d2t
        use model_vars, only: mdb, thp, pip, rvp, thb, pib, rvb, np, mp,nc,mc,mpc
        use thermo_functions, only: calc_rsat

        implicit none

        real :: kappa=0.5       !hygroscopicity parameter, assume 0.5
        real :: rhop = 1400.    !density of particle, assume ammonium sulfate
        
        integer :: ix,iz,ipb,idb ! counters
        real :: th, pi, rv, t, p, rvsat, Samb,Sc ! temporary variables for saturation calc
        real :: md_each,mw_each,mp_each,mdf_each! temporary variables for activation cal
        real :: dpd_each,dpw_each,dpp_each,lambd, Seq, G, Im,A,beta, condensed_water
        real :: t1=0., t2=0., t3=0.
        integer :: jdb

        real :: ncf(npb,ndb),mcf(npb,ndb),mpcf(npb,ndb)

        !first loop over all spatial points
        do iz=2,nz-1
            do ix=2,nx-1
                th = thp(ix,iz,3) + thb(iz)
                pi = pip(ix,iz,3) + pib(iz)
                t = th*pi
                p = p00 * (pi**(cp/rd))
                rv = rvp(ix,iz,3) + rvb(iz)
                rvsat = calc_rsat(t,p)
                Samb = rv/rvsat

                A=4*MW_w*sigma_w/(rd*t*rhol)
                lambd = calc_lambd(p,t)

                do ipb=1,npb
                    do idb=1,ndb
                        ncf(ipb,idb)  = 0.
                        mcf(ipb,idb)  = 0.
                        mpcf(ipb,idb)  = 0.
                    enddo
                enddo


                do ipb=1,npb
                    do idb=1,ndb
                        if ((nc(ix,iz,ipb,idb,3) > 0.) .and. (mc(ix,iz,ipb,idb,3) > 0.)) then
                            !print*,'mc',mc(ix,iz,ipb,idb,3)
                            !print*,'nc',nc(ix,iz,ipb,idb,3)

                            mw_each = mc(ix,iz,ipb,idb,3)/nc(ix,iz,ipb,idb,3) !water mass per particle
                            if (mw_each<0.) then
                                print*,mw_each ! id on't understand why but if i comment this out it doens't work and gets fp error
                            endif 
                            mp_each = mpc(ix,iz,ipb,idb,3)/nc(ix,iz,ipb,idb,3) !dry mass per particle
                            !print*,mp_each
                            md_each = mw_each + mp_each !total cloud droplet mass per particle
                            !print*,'mass',mw_each,mp_each,md_each

                            dpw_each = calc_dp(mw_each,rhol)
                            dpp_each = calc_dp(mp_each,rhop)
                            dpd_each = ((dpw_each**3.) + (dpp_each**3.))**(1/3.)

                            !print*,'dp',dpw_each,dpp_each,dpd_each

                            beta = calc_dahneke(dpd_each, lambd, 1.)
                            !print*,'beta',beta
                            
                            Seq = (dpd_each**3 - dpp_each**3)/(dpd_each**3 - (1-kappa)*dpp_each**3) * EXP(A/dpd_each)
                            !print*,'Seq',Seq

                            t1 = (R*t)/(calc_esat(t)*dg*MW_w)
                            t2 = (lv/(kair*t))
                            t3 = (lv*MW_w/(R*t)-1)
                            if (t1<0.) then
                                print*,t1,t2,t3
                            endif

                            G = t1 + (t2*t3)
                            !print*,'G',G
                            
                            Im = 2*trigpi*dpd_each*beta*(Samb-Seq)*(1/G)

                            mdf_each = md_each + (Im*d2t)
                            !print*,Im,md_each,mdf_each
                        
                            ! if the new mass per droplet is more than the smallest cloud bin
                            if (mdf_each>mdb(1)) then
                                condensed_water = (Im*d2t) * nc(ix,iz,ipb,idb,3) 

                                ! what bin should newly grown (or shrunk) drop go to 
                                jdb = COUNT(mdb(:)<mdf_each)

                                ncf(ipb,jdb) = ncf(ipb,jdb)  +  nc(ix,iz,ipb,idb,3) 
                                mcf(ipb,jdb) = mcf(ipb,jdb)  + condensed_water
                                mpcf(ipb,jdb) = mpcf(ipb,jdb)  +  (mp_each * nc(ix,iz,ipb,idb,3))
                                !print*,'condensed', condensed_water
                            else
                                !otherwise the evaporation has made droplet smaller than smallest possible droplet, return aerosol to environment
                                condensed_water = md_each * nc(ix,iz,ipb,idb,3) 
                                
                                np(ix,iz,ipb,3) = np(ix,iz,ipb,3) + nc(ix,iz,ipb,idb,3)
                                mp(ix,iz,ipb,3) = mp(ix,iz,ipb,3) + (nc(ix,iz,ipb,idb,3)*mp_each)
                            endif

                            ! add the equivalent amount of latent heating to THP
                            thp(ix,iz,3) = thp(ix,iz,3) + condensed_water*(lv/(cp*pib(iz)))
                            ! remove equivalent amount of water from RVP
                            rvp(ix,iz,3) = rvp(ix,iz,3) - condensed_water
                        endif   
                    enddo ! end droplet bin loop
                enddo ! end particle bin loop

                do ipb=1,npb
                    do idb=1,ndb
                        nc(ix,iz,ipb,idb,3)  = ncf(ipb,idb)
                        mc(ix,iz,ipb,idb,3)  = mcf(ipb,idb)
                        mpc(ix,iz,ipb,idb,3)  = mpcf(ipb,idb)
                    enddo
                enddo

            enddo
        enddo
    
    end subroutine condensation

    subroutine check_negs
        use run_constants, only: nx,nz,npb,ndb
        use model_vars, only: np, mp,nc,mc,mpc,nr,mr,mpr
        
        implicit none

        integer :: ix,iz,ipb,idb ! counters
        
        do ix=2,nx-1
            do iz=2,nz-1
                do ipb=1,npb
                    if ((np(ix,iz,ipb,3)<0.) .or. (mp(ix,iz,ipb,3)<0.)) then
                        np(ix,iz,ipb,3)=0.
                        mp(ix,iz,ipb,3)=0.
                    endif

                    do idb=1,ndb
                        if ((nc(ix,iz,ipb,idb,3)<0.) .or. (mc(ix,iz,ipb,idb,3)<0.) .or. (mpc(ix,iz,ipb,idb,3)<0.)) then
                            nc(ix,iz,ipb,idb,3)=0.
                            mc(ix,iz,ipb,idb,3)=0.
                            mpc(ix,iz,ipb,idb,3) = 0.
                        endif

                        if ((nr(ix,iz,ipb,idb,3)<0.) .or. (mr(ix,iz,ipb,idb,3)<0.) .or. (mpr(ix,iz,ipb,idb,3)<0.)) then
                            nr(ix,iz,ipb,idb,3)=0.
                            mr(ix,iz,ipb,idb,3)=0.
                            mpr(ix,iz,ipb,idb,3) = 0.
                        endif
                    enddo
                enddo
            enddo
        enddo

    end subroutine check_negs

end module aerosol
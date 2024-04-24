module microphysics
    implicit none

    contains
    
    subroutine microphysics_driv
        !this is the main microphysics driver
        implicit none
    
        call check_negs
        print*,'check negs'

        call activation
        print*,'activation'

        call check_negs
        print*,'check negs'

        call condensation
        print*,'condensation'

        !call check_negs
        !print*,'check negs'

        call collisioncoalescence 
        print*,'collision-coalescence'
         
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
        use run_constants, only: nx,nz,npartbin,ndropbin, kappa, rhop
        use model_vars, only: mdropbin_lims, thp, pip, rvp, thb, pib, rvb, np, mp,nc,mlc,mpc
        use thermo_functions, only: calc_rsat,calc_satfrac
        use aerosol, only: calc_mass, calc_dp, calc_scrit

        implicit none
        
        integer :: ix,iz,ipartbin,idropbin ! counters
        real :: th, pi, rv, t, p, rvsat, Samb ! temporary variables for saturation calc
        real :: mp_each,dp_each,Sc_each ! temporary variables for activation calc

        !first loop over all spatial points
        do iz=2,nz-1
            do ix=2,nx-1
                ! calculate actual values instead of perturbation values of theta, pressure, rv
                th = thp(ix,iz,3) + thb(iz)
                pi = pip(ix,iz,3) + pib(iz)
                t = th*pi
                p = p00 * (pi**(cp/rd))
                rv = rvp(ix,iz,3) + rvb(iz)
                rvsat = calc_rsat(t,p)
                Samb = calc_satfrac(rv,rvsat)

                ! if the grid point is supersaturated and
                if (Samb>1.) then
                    ! and there are fewer cloud droplets than available CCN locally
                    ! this is what HUCM does -- I don't think it is that well physically-justified, but ok for now
                    if (SUM(nc(ix,iz,:,:,3)<SUM(np(ix,iz,:,3)))) then
                        ! then we should activate particles! 
                        do ipartbin=1,npartbin
                            ! for each particle bin, check if there are any particles at all (saves us some time looping over bins that have no aerosol)
                            if ((np(ix,iz,ipartbin,3) > 0.) .and. (mp(ix,iz,ipartbin,3) > 0.)) then
                                mp_each = mp(ix,iz,ipartbin,3)/np(ix,iz,ipartbin,3) ! mean mass of aerosol particle in bin
                                dp_each = calc_dp(mp_each, rhop) ! mean particle diameter in bin
                                Sc_each = calc_scrit(dp_each, t, kappa) !corresponding mean critical saturation in bin

                                if (Samb>Sc_each) then !if ambient supersaturation is above critical, do activation
                                    ! determine which cloud bin the activate droplets should go into
                                    if ((5*mp_each) < mdropbin_lims(1)) then
                                        ! if the mean aerosol size is way smaller than the first bin, just shove them in first droplet bin
                                        idropbin = 1
                                    else
                                        ! else, find the appropriate bin where droplet size is just smaller than our activated particle
                                        idropbin = COUNT(mdropbin_lims(:)<(5*mp_each))
                                    endif
                                    
                                    ! add liquid number & mass & aerosol mass
                                    ! we assume all the particles within a bin activate at the same time if they hit the mean critical supersat for the bin
                                    nc(ix,iz,ipartbin,idropbin,3) = nc(ix,iz,ipartbin,idropbin,3) + np(ix,iz,ipartbin,3) 
                                    mlc(ix,iz,ipartbin,idropbin,3) = mlc(ix,iz,ipartbin,idropbin,3) + (mdropbin_lims(idropbin)-mp_each)
                                    mpc(ix,iz,ipartbin,idropbin,3) = mpc(ix,iz,ipartbin,idropbin,3) + (np(ix,iz,ipartbin,3)*mp_each)

                                    ! add the equivalent amount of latent heating to THP
                                    thp(ix,iz,3) = thp(ix,iz,3) + (mdropbin_lims(idropbin)-mp_each)*(lv/(cp*pib(iz)))
                                    ! remove equivalent amount of water from RVP
                                    rvp(ix,iz,3) = rvp(ix,iz,3) - (mdropbin_lims(idropbin)-mp_each)

                                    ! remove aerosol from unprocessed aerosol bins
                                    np(ix,iz,ipartbin,3) = 0. 
                                    mp(ix,iz,ipartbin,3) = 0.
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
        use run_constants, only: nx,nz,npartbin,ndropbin,dt,d2t,kappa,rhop
        use model_vars, only: mdropbin_lims, thp, pip, rvp, thb, pib, rvb, np, mp,nc,mlc,mpc
        use thermo_functions, only: calc_rsat, calc_esat
        use aerosol, only: calc_mass, calc_dp, calc_lambd, calc_dahneke

        implicit none

        integer :: ix,iz,ipartbin,idropbin ! counters
        integer :: jdropbin ! index of new bin after condensation

        real :: th, pi, rv, t, p, rvsat, Samb,Sc ! temporary variables for saturation calc
        real :: md_each_i,ml_each_i,mp_each_i ! per droplet in bin, mass total, liquid mass, and particle mass
        real :: md_each_j ! final mass per droplet after condensation
        real :: dpd_each,dpw_each,dpp_each,lambd, Seq, G, Im,A,beta, condensed_water
        real :: t1=0., t2=0., t3=0. !more temporary variables
        real :: minmass = 1.e-30

        real :: nc_final(npartbin,ndropbin),mlc_final(npartbin,ndropbin),mpc_final(npartbin,ndropbin)

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

                ! calculate some terms for Kohler curve later -- these only depend on temp and pressure, so do them per spatial point
                A=4*MW_w*sigma_w/(rd*t*rhol)
                lambd = calc_lambd(p,t)

                ! initialize the output arrays with zero
                do ipartbin=1,npartbin
                    do idropbin=1,ndropbin
                        nc_final(ipartbin,idropbin)  = 0.
                        mlc_final(ipartbin,idropbin)  = 0.
                        mpc_final(ipartbin,idropbin)  = 0.
                    enddo
                enddo

                do ipartbin=1,npartbin
                    do idropbin=1,ndropbin
                        if ((nc(ix,iz,ipartbin,idropbin,3) > 0.) .and. (mlc(ix,iz,ipartbin,idropbin,3) > 0.)) then
                            print*,mlc(ix,iz,ipartbin,idropbin,3),nc(ix,iz,ipartbin,idropbin,3)
                            ml_each_i = mlc(ix,iz,ipartbin,idropbin,3)/nc(ix,iz,ipartbin,idropbin,3) !water mass per particle

                            if (ml_each_i<0.) then
                                print*,ml_each_i ! i don't understand why but if i comment this out it doens't work and gets fp error
                            endif 

                            print*,mpc(ix,iz,ipartbin,idropbin,3),nc(ix,iz,ipartbin,idropbin,3)
                            mp_each_i = mpc(ix,iz,ipartbin,idropbin,3)/nc(ix,iz,ipartbin,idropbin,3) !dry mass per particle

                            print*,ml_each_i,mp_each_i
                                
                            if ((ml_each_i<minmass) .or. (mp_each_i<minmass)) then
                                print*,ml_each_i 
                            else
                                md_each_i = ml_each_i + mp_each_i !total cloud droplet mass per particle

                                ! calculate the equivalent diameter for the mass of water, mass of particle
                                ! then calculate the resulting droplet diameter as cubed sum

                                dpw_each = calc_dp(ml_each_i,rhol)
                                dpp_each = calc_dp(mp_each_i,rhop)
                                dpd_each = ((dpw_each**3.) + (dpp_each**3.))**(1/3.)
                                
                                ! now that we know droplet size, do full Kohler calculation for Seq 
                                Seq = (dpd_each**3 - dpp_each**3)/(dpd_each**3 - (1-kappa)*dpp_each**3) * EXP(A/dpd_each)
                                
                                ! calculating the amount of condensation following Pruppracher & Klemt
                                ! this is not elegant code but if I don't have this print statement then I get non-normal values
                                t1 = (R*t)/(calc_esat(t)*dg*MW_w)
                                t2 = (lv/(kair*t))
                                t3 = (lv*MW_w/(R*t)-1)
                                if (t1<0.) then
                                    print*,t1,t2,t3
                                endif

                                G = t1 + (t2*t3)
                                
                                ! calculate dahneke parameter
                                beta = calc_dahneke(dpd_each, lambd, 1.)

                                ! amount of water condensed this timestep
                                Im = 2*trigpi*dpd_each*beta*(Samb-Seq)*(1/G)

                                print*,'rate',Im
                                ! the final amount of droplet mass will be the original mass + Im * timestep
                                md_each_j = md_each_i + (Im*d2t)
                                
                                ! if the new mass per droplet is more than the smallest cloud bin
                                if (md_each_j>mdropbin_lims(1)) then
                                    ! then we need to find which bin to put the new droplet in 
                                    condensed_water = (Im*d2t) * nc(ix,iz,ipartbin,idropbin,3) 

                                    ! what bin should newly grown (or shrunk) drop go to 
                                    jdropbin = COUNT(mdropbin_lims(:)<md_each_j)

                                    print*,'new bin',jdropbin

                                    ! assign these post-condensation values to temporary arrays
                                    nc_final(ipartbin,jdropbin) = nc_final(ipartbin,jdropbin)  +  nc(ix,iz,ipartbin,idropbin,3) 
                                    mlc_final(ipartbin,jdropbin) = mlc_final(ipartbin,jdropbin)  + condensed_water
                                    mpc_final(ipartbin,jdropbin) = mpc_final(ipartbin,jdropbin)  +  (mp_each_i * nc(ix,iz,ipartbin,idropbin,3))
                                    
                                else
                                    !otherwise the evaporation has made droplet smaller than smallest possible droplet, return aerosol to environment
                                    condensed_water = md_each_i * nc(ix,iz,ipartbin,idropbin,3) 
                                    
                                    np(ix,iz,ipartbin,3) = np(ix,iz,ipartbin,3) + nc(ix,iz,ipartbin,idropbin,3)
                                    mp(ix,iz,ipartbin,3) = mp(ix,iz,ipartbin,3) + (nc(ix,iz,ipartbin,idropbin,3)*mp_each_i)
                                endif

                                print*,'moved to new bin'

                                ! add the equivalent amount of latent heating to THP
                                thp(ix,iz,3) = thp(ix,iz,3) + condensed_water*(lv/(cp*pib(iz)))
                                ! remove equivalent amount of water from RVP
                                rvp(ix,iz,3) = rvp(ix,iz,3) - condensed_water

                                print*,'adjust th and rv'
                            endif
                        endif   
                    enddo ! end droplet bin loop
                enddo ! end particle bin loop

                ! move from temporary arrays to final
                do ipartbin=1,npartbin
                    do idropbin=1,ndropbin
                        nc(ix,iz,ipartbin,idropbin,3)  = nc_final(ipartbin,idropbin)
                        mlc(ix,iz,ipartbin,idropbin,3)  = mlc_final(ipartbin,idropbin)
                        mpc(ix,iz,ipartbin,idropbin,3)  = mpc_final(ipartbin,idropbin)
                    enddo
                enddo
            
            enddo
        enddo

        print*,'done!'
        do ix=2,nx-1
            do iz=2,nz-1
                do ipartbin=1,npartbin
                    do idropbin=1,ndropbin
                        if ((nc(ix,iz,ipartbin,idropbin,3)<0.) .or. (mlc(ix,iz,ipartbin,idropbin,3)<0.) .or. (mpc(ix,iz,ipartbin,idropbin,3)<0.)) then
                            print*,nc(ix,iz,ipartbin,idropbin,3) ,mlc(ix,iz,ipartbin,idropbin,3) ,mpc(ix,iz,ipartbin,idropbin,3) 
                            nc(ix,iz,ipartbin,idropbin,3) = 0.
                            mlc(ix,iz,ipartbin,idropbin,3) = 0.
                            mpc(ix,iz,ipartbin,idropbin,3) = 0.
                            print*,nc(ix,iz,ipartbin,idropbin,3) ,mlc(ix,iz,ipartbin,idropbin,3) ,mpc(ix,iz,ipartbin,idropbin,3) 
                        endif
                    enddo
                enddo
            enddo
        enddo

        print*,'negs ok'
    end subroutine condensation

    subroutine collisioncoalescence 
        use constants, only: cp, rd, p00, lv,rhol,MW_w,sigma_w,trigpi,dg,kair,R
        use run_constants, only: nx,nz,npartbin,ndropbin,dt,d2t,rhop
        use model_vars, only: mdropbin_lims, mpartbin_lims,thp, pip, thb, pib, np, mp,nc,mlc,mpc
        use thermo_functions, only: calc_rhoair
        use aerosol, only: calc_mass, calc_dp

        implicit none

        integer :: ix,iz,ipb,jpb,kpb,idb,jdb,kdb ! counters
        
        real :: rho_air, th, pi, t, p , rv, rvsat ! extra variables for calculations

        real :: ncf(npartbin,ndropbin),mcf(npartbin,ndropbin),mpcf(npartbin,ndropbin) !final cloud arrays
        real :: nct(ndropbin), mct(ndropbin), mpct(ndropbin) ! total cloud values 
        real :: mw_each_i,mp_each_i,md_each_i,mw_each_j,mp_each_j,md_each_j &
                ,frac_i(npartbin),frac_j(npartbin),vol_each_i,vol_each_j,dpd_each_j   &
                ,pj,K,Jij,md_each_k,mp_each_k,mpc_i(npartbin),mpc_j(npartbin) !temporary variables for collision coalescence calc


        real :: minmass = 1.e-30

        print*,'start coagulation'
        !first loop over all spatial points
        do iz=2,nz-1
            do ix=2,nx-1
                th = thp(ix,iz,3) + thb(iz)
                pi = pip(ix,iz,3) + pib(iz)
                t = th*pi
                p = p00 * (pi**(cp/rd))
                rho_air = calc_rhoair(t,p)
                print*,'initial temps calc'

                ! calculate total cloud number, water mass, aerosol mass in each cloud drop bin
                do idb =1,ndropbin
                    nct(idb) = SUM(nc(ix,iz,:,idb,3))
                    mct(idb) = SUM(mlc(ix,iz,:,idb,3))
                    mpct(idb) = SUM(mpc(ix,iz,:,idb,3))
                enddo 
                print*,'total calc'

                ! initially set the final variables to be same as initial (post-advection and condensation, but pre-collision coalescence)
                ncf(:,:) = nc(ix,iz,:,:,3)
                mcf(:,:) = mlc(ix,iz,:,:,3)
                mpcf(:,:) = mpc(ix,iz,:,:,3)

                print*,'set pre-coagulation values'
                
                do idb = 1,ndropbin
                    if (((mct(idb)+mpct(idb))> 0.) .and. (nct(idb)>0.)) then 
                        do jdb = idb,ndropbin
                            if (((mct(jdb)+mpct(jdb))>0.) .and. (nct(jdb)>0.)) then
                                print*,nct(idb),nct(jdb)

                                mw_each_i = mct(idb)/nct(idb) !water mass per drop in bin i
                                if (mw_each_i<0.) then
                                    print*,mw_each_i ! i don't understand why but if i comment this out it doens't work and gets fp error
                                endif 
                                mp_each_i = mpct(idb)/nct(idb) !dry mass per particle
                                md_each_i = mw_each_i + mp_each_i !total cloud droplet mass per particle
                                
                                mw_each_j = mct(jdb)/nct(jdb) !water mass per drop in bin j
                                if (mw_each_j<0.) then
                                    print*,mw_each_j ! i don't understand why but if i comment this out it doens't work and gets fp error
                                endif 
                                mp_each_j = mpct(jdb)/nct(jdb) !dry mass per particle
                                md_each_j = mw_each_j + mp_each_j !total cloud droplet mass per particle

                                print*,'initial droplet sizes calc'

                                if ((mw_each_i<minmass) .or. (mp_each_i<minmass) .or. (mw_each_j<minmass) .or. (mp_each_j<minmass)) then
                                    print*,mw_each_i
                                else
                                    !calculate total volume of each particle
                                    vol_each_i = ((mw_each_i/rhol) + (mp_each_i/rhop)) * 100.**3
                                    vol_each_j = ((mw_each_j/rhol) + (mp_each_j/rhop)) * 100.**3


                                    ! calculate larger droplet diameter
                                    dpd_each_j = (calc_dp(mw_each_j,rhol)**3 + calc_dp(mp_each_j,rhop)**3)**(1/3.)

                                    print*,'each size calc'

                                    if (dpd_each_j > 100.e-6) then
                                        pj = 5.78e3 * (vol_each_i + vol_each_j)
                                    else 
                                        pj = 9.44e9 * (vol_each_i**2 + vol_each_j**2)
                                    end if

                                    K = pj * rho_air * (1/100.)**3

                                    ! calculate coagulation rate
                                    Jij = K * nct(idb) * nct(jdb)

                                    print*,'rates calc'

                                    ! calculate fraction of cloud in each aerosol bin for distributing later
                                    frac_i(:) = nc(ix,iz,:,idb,3)/nct(idb)
                                    frac_j(:) = nc(ix,iz,:,jdb,3)/nct(jdb)

                                    do ipb=1,npartbin
                                        if (nc(ix,iz,ipb,idb,3)>0.) then
                                            mpc_i(ipb) = mpc(ix,iz,ipb,idb,3)/nc(ix,iz,ipb,idb,3)
                                        else
                                            mpc_i(ipb) = mpartbin_lims(ipb)
                                        endif
                                    enddo 

                                    do jpb=1,npartbin
                                        if (nc(ix,iz,jpb,idb,3)>0.) then
                                            mpc_j(jpb) = mpc(ix,iz,jpb,idb,3)/nc(ix,iz,jpb,idb,3)
                                        else
                                            mpc_j(jpb) = mpartbin_lims(jpb)
                                        endif
                                    enddo 

                                    ! redistribute cloud number/mass based on collision-coalescence rate
                                    ncf(:,idb) = ncf(:,idb) - (Jij*d2t*frac_i(:))
                                    ncf(:,jdb) = ncf(:,idb) - (Jij*d2t*frac_j(:))

                                    mcf(:,idb) = mcf(:,idb) - (Jij*d2t*frac_i(:)*(md_each_i-mpc_i(:)))
                                    mcf(:,jdb) = mcf(:,jdb) - (Jij*d2t*frac_j(:)*(md_each_j-mpc_j(:)))

                                    mpcf(:,idb) = mpcf(:,idb) - (Jij*d2t*frac_i(:)*mpc_i(:))
                                    mpcf(:,jdb) = mpcf(:,jdb) - (Jij*d2t*frac_j(:)*mpc_j(:))

                                    print*,'subtract from old bins'

                                    ! find new droplet mass bin
                                    md_each_k = md_each_j + md_each_j
                                    kdb = COUNT(mdropbin_lims(:)<md_each_k) 
                                    
                                    if (md_each_k>mdropbin_lims(ndropbin)) then 
                                        print*,'too big droplet!'
                                    endif

                                    print*,'calc new bin'

                                    do ipb=1,npartbin
                                        if (frac_i(ipb)>0.) then
                                            do jpb=1,npartbin 
                                                if (frac_j(jpb)>0.) then
                                                    mp_each_k = mpc_i(ipb) + mpc_j(jpb)

                                                    ! will need to make particle bin dimension
                                                    kpb = COUNT(mpartbin_lims(:)<mp_each_k)

                                                    if (mp_each_k > md_each_k) then
                                                        print*,'more aerosol than water! smth is wrong'
                                                    endif

                                                    ncf(kpb,kdb) = ncf(kpb,kdb) + (Jij*d2t*frac_i(ipb)*frac_j(jpb))
                                                    mcf(kpb,kdb) = mcf(kpb,kdb) + (Jij*d2t*frac_i(ipb)*frac_j(jpb)*(md_each_k-mp_each_k))
                                                    mpcf(kpb,kdb) = mpcf(kpb,kdb) + (Jij*d2t*frac_i(ipb)*frac_j(jpb)*mp_each_k)
                                                endif 
                                            enddo
                                        endif
                                    enddo
                                endif
                            endif 
                        enddo
                    endif 
                enddo

                print*,'start negs check'

                

                do ipb=1,npartbin
                    do idb=1,ndropbin
                        !print*,ncf(ipb,idb),mcf(ipb,idb),mpcf(ipb,idb)
                        !negative values check
                        if ((ncf(ipb,idb)<0.) .or. (mcf(ipb,idb)<0.) .or. (mpcf(ipb,idb)<0.)) then
                            print*,ncf(ipb,idb) ,mcf(ipb,idb) ,mpcf(ipb,idb) 
                        
                            ncf(ipb,idb) = 0.
                            mcf(ipb,idb) = 0.
                            mpcf(ipb,idb) = 0.
                            print*,ncf(ipb,idb) ,mcf(ipb,idb) ,mpcf(ipb,idb) 
                        endif

                        print*,'put in main arrays'
                        ! move from temporary arrays to main data arrays
                        nc(ix,iz,ipb,idb,3)  = ncf(ipb,idb)
                        mlc(ix,iz,ipb,idb,3)  = mcf(ipb,idb)
                        mpc(ix,iz,ipb,idb,3)  = mpcf(ipb,idb)

                        print*,'main arrays ok'

                        if ((nc(ix,iz,ipb,idb,3)<0.) .or. (mlc(ix,iz,ipb,idb,3)<0.) .or. (mpc(ix,iz,ipb,idb,3)<0.)) then
                            print*,nc(ix,iz,ipb,idb,3) ,mlc(ix,iz,ipb,idb,3) ,mpc(ix,iz,ipb,idb,3) 
                            nc(ix,iz,ipb,idb,3) = 0.
                            mlc(ix,iz,ipb,idb,3) = 0.
                            mpc(ix,iz,ipb,idb,3) = 0.
                            print*,nc(ix,iz,ipb,idb,3) ,mlc(ix,iz,ipb,idb,3) ,mpc(ix,iz,ipb,idb,3) 
                        endif

                        !print*,nc(ix,iz,ipb,idb,3),mc(ix,iz,ipb,idb,3),mpc(ix,iz,ipb,idb,3)
                    enddo
                enddo

                print*,'end negs check'
            enddo
        enddo

        print*,'all collision done'

    end subroutine collisioncoalescence
    

    subroutine check_negs
        use run_constants, only: nx,nz,npartbin,ndropbin
        use model_vars, only: np, mp,nc,mlc,mpc
        
        implicit none

        integer :: ix,iz,ipb,idb ! counters
        real:: minmass=1.e-30
        
        do ix=2,nx-1
            do iz=2,nz-1
                do ipb=1,npartbin
                    if ((np(ix,iz,ipb,3)<0.) .or. (mp(ix,iz,ipb,3)<minmass)) then
                        np(ix,iz,ipb,3)=0.
                        mp(ix,iz,ipb,3)=0.
                    endif

                    do idb=1,ndropbin
                        if ((nc(ix,iz,ipb,idb,3)<0.) .or. (mlc(ix,iz,ipb,idb,3)<minmass) .or. (mpc(ix,iz,ipb,idb,3)<minmass)) then
                            nc(ix,iz,ipb,idb,3)=0.
                            mlc(ix,iz,ipb,idb,3)=0.
                            mpc(ix,iz,ipb,idb,3) = 0.
                        endif
                    enddo
                enddo
            enddo
        enddo

    end subroutine check_negs
end module microphysics 
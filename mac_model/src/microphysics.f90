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

        call check_negs
        print*,'check negs'

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
        use model_vars, only: mdropbin_lims, thp, pip, rvp, thb, pib, rvb, np, mp, nd, mld, mpd
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
                    if (SUM(nd(ix,iz,:,:,3)<SUM(np(ix,iz,:,3)))) then
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
                                        idropbin = COUNT(mdropbin_lims(:ndropbin)<=(5*mp_each))
                                    endif

                                    ! add liquid number & mass & aerosol mass
                                    ! we assume all the particles within a bin activate at the same time if they hit the mean critical supersat for the bin
                                    nd(ix,iz,ipartbin,idropbin,3) = nd(ix,iz,ipartbin,idropbin,3) + np(ix,iz,ipartbin,3) 
                                    mld(ix,iz,ipartbin,idropbin,3) = mld(ix,iz,ipartbin,idropbin,3) + (mdropbin_lims(idropbin)-mp_each)*np(ix,iz,ipartbin,3)
                                    mpd(ix,iz,ipartbin,idropbin,3) = mpd(ix,iz,ipartbin,idropbin,3) + (np(ix,iz,ipartbin,3)*mp_each)

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
        use model_vars, only: mdropbin_lims,mpartbin_lims, thp, pip, rvp, thb, pib, rvb, np, mp,nd,mld,mpd
        use thermo_functions, only: calc_rsat, calc_esat, calc_satfrac
        use aerosol, only: calc_mass, calc_dp, calc_lambd, calc_dahneke

        implicit none

        integer :: ix,iz,ipartbin,idropbin ! counters
        integer :: jdropbin ! index of new bin after condensation

        real :: th, pi, rv, t, p, rvsat, Samb,Sc ! temporary variables for saturation calc
        real :: md_each_i,mld_each_i,mpd_each_i ! per droplet in bin, mass total, liquid mass, and particle mass
        real :: md_each_j ! final mass per droplet after condensation
        real :: dpd_each_i,dpld_each_i,dppd_each_i,lambd, Seq, G, Im,A,beta, condensed_water
        real :: t1=0., t2=0., t3=0. !more temporary variables
        real :: minmass = 1.e-30

        real :: nd_final(npartbin,ndropbin),mld_final(npartbin,ndropbin),mpd_final(npartbin,ndropbin)
        real :: temp
        !first loop over all spatial points
        do iz=2,nz-1
            do ix=2,nx-1
                th = thp(ix,iz,3) + thb(iz)
                pi = pip(ix,iz,3) + pib(iz)
                t = th*pi
                p = p00 * (pi**(cp/rd))
                rv = rvp(ix,iz,3) + rvb(iz)
                rvsat = calc_rsat(t,p)
                Samb = calc_satfrac(rv,rvsat)

                ! calculate some terms for Kohler curve later -- these only depend on temp and pressure, so do them per spatial point
                A=4*MW_w*sigma_w/(rd*t*rhol)
                lambd = calc_lambd(p,t)

                ! initialize the output arrays with zero
                do ipartbin=1,npartbin
                    do idropbin=1,ndropbin
                        nd_final(ipartbin,idropbin)  = 0.
                        mld_final(ipartbin,idropbin)  = 0.
                        mpd_final(ipartbin,idropbin)  = 0.
                    enddo
                enddo

                do ipartbin=1,npartbin
                    do idropbin=1,ndropbin
                        ! check if there are any droplets in this bin at all, if not we can skip it and save calculation time
                        if ((nd(ix,iz,ipartbin,idropbin,3) > 0.) .and. ((mld(ix,iz,ipartbin,idropbin,3)+mpd(ix,iz,ipartbin,idropbin,3)) > 0.)) then
                            ! mass of liquid in each particle/drop bin
                            mld_each_i = mld(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3) ! mass of liquid in droplet per particle/drop bin
                            if (mld_each_i<0.) then
                                print*,mld_each_i ! i don't understand why but if i comment this out it doesn't work and gets fp error
                            endif 

                            ! dry particle mass in each particle/drop bin
                            mpd_each_i = mpd(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3) 

                            if (mpd_each_i<0.) then
                                print*,mpd_each_i ! i don't understand why but if i comment this out it doesn't work and gets fp error
                            endif 
  
                            if ((mld_each_i<minmass) .or. (mpd_each_i<minmass)) then
                                ! check that the mass is above some threshold, if not the calc_dp function returns 0 and there is a problem with saturation calc
                                print*,'too small mass, skip this bin'
                            else
                                ! total cloud mass in each particle/drop bin
                                md_each_i = mld_each_i + mpd_each_i 

                                ! calculate the equivalent diameter for the mass of water, mass of particle
                                ! then calculate the resulting droplet diameter as cubed sum
                                dpld_each_i = calc_dp(mld_each_i,rhol)
                                dppd_each_i = calc_dp(mpd_each_i,rhop)
                                dpd_each_i = ((dpld_each_i**3.) + (dppd_each_i**3.))**(1/3.)
                                
                                ! now that we know droplet size, do full Kohler calculation for Seq 
                                Seq = (dpd_each_i**3 - dppd_each_i**3)/(dpd_each_i**3 - ((1-kappa)*dppd_each_i**3)) * EXP(A/dpd_each_i)
                                
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
                                beta = calc_dahneke(dpd_each_i, lambd, 1.)

                                ! amount of water condensed this timestep
                                Im = 2.*trigpi*dpd_each_i*beta*(Samb-Seq)*(1/G)

                                ! the final amount of droplet mass will be the original mass + Im * 2*dt
                                md_each_j = md_each_i + (Im*d2t)
                                
                                ! TO DO:  add a check to say we can't condense more than available water vapor

                                ! if the new mass per droplet is more than the smallest cloud bin
                                if (md_each_j>mdropbin_lims(1)) then
                                    ! total amount of condensed water
                                    condensed_water = (Im*d2t) * nd(ix,iz,ipartbin,idropbin,3) 

                                    ! what bin should newly grown (or shrunk) drop go to 
                                    jdropbin = COUNT(mdropbin_lims(:ndropbin)<=md_each_j)

                                    
                                    if (jdropbin<ndropbin) then
                                        ! assign these post-condensation values to their new bin
                                        nd_final(ipartbin,jdropbin) = nd_final(ipartbin,jdropbin)  +  nd(ix,iz,ipartbin,idropbin,3) 
                                        mld_final(ipartbin,jdropbin) = mld_final(ipartbin,jdropbin)  + mld(ix,iz,ipartbin,idropbin,3) + condensed_water
                                        ! this works because condensation doesn't change the mass of dry particle, so we're just moving around the particle mass
                                        mpd_final(ipartbin,jdropbin) = mpd_final(ipartbin,jdropbin)  +  (mpd_each_i * nd(ix,iz,ipartbin,idropbin,3))
                                    endif
                                else
                                    !otherwise the evaporation has made droplet smaller than smallest possible droplet, return aerosol to environment
                                    condensed_water = -(md_each_i-mpd_each_i) * nd(ix,iz,ipartbin,idropbin,3) 
                                    
                                    np(ix,iz,ipartbin,3) = np(ix,iz,ipartbin,3) + nd(ix,iz,ipartbin,idropbin,3)
                                    mp(ix,iz,ipartbin,3) = mp(ix,iz,ipartbin,3) + (nd(ix,iz,ipartbin,idropbin,3)*mpd_each_i)
                                endif

                                ! add the equivalent amount of latent heating to THP
                                thp(ix,iz,3) = thp(ix,iz,3) + condensed_water*(lv/(cp*pib(iz)))
                                ! remove equivalent amount of water from RVP
                                rvp(ix,iz,3) = rvp(ix,iz,3) - condensed_water

                            endif
                        endif   
                    enddo ! end droplet bin loop
                enddo ! end particle bin loop

                ! move from temporary arrays to final
                do ipartbin=1,npartbin
                    do idropbin=1,ndropbin
                        nd(ix,iz,ipartbin,idropbin,3)  = nd_final(ipartbin,idropbin)
                        mld(ix,iz,ipartbin,idropbin,3)  = mld_final(ipartbin,idropbin)
                        mpd(ix,iz,ipartbin,idropbin,3)  = mpd_final(ipartbin,idropbin)

                        
                    enddo ! end droplet bin loop
                enddo ! end particle bin loop
    
            enddo ! end x loop
        enddo ! end z loop

    end subroutine condensation

    subroutine collisioncoalescence 
        use constants, only: cp, rd, p00, lv,rhol,MW_w,sigma_w,trigpi,dg,kair,R
        use run_constants, only: nx,nz,npartbin,ndropbin,dt,d2t,rhop
        use model_vars, only: mdropbin_lims, mpartbin_lims,thp, pip, thb, pib, np, mp,nd,mld,mpd
        use thermo_functions, only: calc_rhoair
        use aerosol, only: calc_mass, calc_dp

        implicit none

        integer :: ix,iz,idropbin,jdropbin,kdropbin,ipartbin,jpartbin,kpartbin ! counters
        ! i and j are the bins doing the collision-coalescence (where j is the bigger drop)
        ! k is the final bin they end up in after collision
        
        real :: rho_air, th, pi, t, p , rv, rvsat ! extra variables for calculations

        real :: nd_final(npartbin,ndropbin),mld_final(npartbin,ndropbin),mpd_final(npartbin,ndropbin) !temporary cloud arrays
        real :: nd_total(ndropbin), mld_total(ndropbin), mpd_total(ndropbin) ! total cloud values in each droplet bin
        real :: mld_each_i,mpd_each_i,md_each_i,vol_each_i                  &
                ,mld_each_j,mpd_each_j,md_each_j,vol_each_j,dpd_each_j      &
                ,mld_each_k,mpd_each_k,md_each_k                            &
                ,frac_n_i(npartbin),frac_n_j(npartbin),frac_mp_i(npartbin),frac_mp_j(npartbin)                        & ! fraction of droplets of this droplet size that are in each particle bin
                ,pj,Kij,Jij                                               !temporary variables for collision coalescence calc


        real :: minmass = 1.e-30
        real :: temp

        !first loop over all spatial points
        do iz=2,nz-1
            do ix=2,nx-1
                th = thp(ix,iz,3) + thb(iz)
                pi = pip(ix,iz,3) + pib(iz)
                t = th*pi
                p = p00 * (pi**(cp/rd))
                rho_air = calc_rhoair(t,p)
                !print*,'initial temps calc'

                ! calculate total cloud number, water mass, aerosol mass in each cloud drop bin
                do idropbin =1,ndropbin
                    nd_total(idropbin) = SUM(nd(ix,iz,:,idropbin,3))
                    mld_total(idropbin) = SUM(mld(ix,iz,:,idropbin,3))
                    mpd_total(idropbin) = SUM(mpd(ix,iz,:,idropbin,3))
                enddo 
                !print*,'total calc'

                ! initially set the final variables to values post-advection and condensation, but pre-collision coalescence
                nd_final(:,:) = nd(ix,iz,:,:,3)
                mld_final(:,:) = mld(ix,iz,:,:,3)
                mpd_final(:,:) = mpd(ix,iz,:,:,3)
                !print*,'set pre-coagulation values'

                
                !loop over first droplet (smaller) bins
                do idropbin = 1,ndropbin
                    ! do a check to see if there's any mass/nuber in this bin -- if not we can skip and save computations
                    if (((mld_total(idropbin)+mpd_total(idropbin))> 0.) .and. (nd_total(idropbin)>0.)) then 
                        mld_each_i = mld_total(idropbin)/nd_total(idropbin) !water mass per drop in bin i
                        
                        if (mld_each_i<0.) then
                            print*,mld_each_i ! i don't understand why but if i comment this out it doesn't work and gets fp error
                        endif 

                        mpd_each_i = mpd_total(idropbin)/nd_total(idropbin) !dry mass per drop in bin i
                        md_each_i = mld_each_i + mpd_each_i !total droplet mass per drop in bin i

                        ! loop over second (bigger) droplet bin

                        do jdropbin = idropbin,ndropbin
                            ! do a check to see if there's any mass/nuber in this bin -- if not we can skip and save computations
                            if (((mld_total(jdropbin)+mpd_total(jdropbin))>0.) .and. (nd_total(jdropbin)>0.)) then
                                
                                mld_each_j = mld_total(jdropbin)/nd_total(jdropbin) !water mass per drop in bin j
                                if (mld_each_j<0.) then
                                    print*,mld_each_j ! i don't understand why but if i comment this out it doesn't work and gets fp error
                                endif 
                                mpd_each_j = mpd_total(jdropbin)/nd_total(jdropbin) !dry mass per drop in bin j
                                md_each_j = mld_each_j + mpd_each_j !total droplet mass per drop in bin j

                                ! if any of the masses are too small, ignore these collisions
                                if ((mld_each_i<minmass) .or. (mpd_each_i<minmass) .or. (mld_each_j<minmass) .or. (mpd_each_j<minmass)) then
                                    print*,mld_each_i
                                else
                                    ! calculate larger droplet diameter
                                    ! collision-coalescence rate is parameterized different for small and big collector drops (below/above 100micron) 
                                    dpd_each_j = (calc_dp(mld_each_j,rhol)**3 + calc_dp(mpd_each_j,rhop)**3)**(1/3.)

                                    !calculate total volume of each particle
                                    vol_each_i = ((mld_each_i/rhol) + (mpd_each_i/rhop)) * 100.**3
                                    vol_each_j = ((mld_each_j/rhol) + (mpd_each_j/rhop)) * 100.**3

                                    if (dpd_each_j > 100.e-6) then
                                        pj = 5.78e3 * (vol_each_i + vol_each_j)
                                    else 
                                        pj = 9.44e9 * (vol_each_i**2 + vol_each_j**2)
                                    end if

                                    ! collision-coalescence rate constant
                                    Kij = pj * rho_air * (1/100.)**3

                                    ! calculate coagulation rate between i and j bins
                                    Jij = Kij * nd_total(idropbin) * nd_total(jdropbin)

                                    !print*,'rates calc'
                                    
                                    do ipartbin=1,npartbin
                                        ! calculate fraction of cloud number/particle mass that are in each particle bin (for a given droplet bin) for distributing later
                                        frac_n_i(ipartbin) = nd(ix,iz,ipartbin,idropbin,3)/nd_total(idropbin)

                                        if ((mpd(ix,iz,ipartbin,idropbin,3)>0.) .and. (nd(ix,iz,ipartbin,idropbin,3)>0)) then
                                            frac_mp_i(ipartbin) = mpd(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3)
                                        else
                                            frac_mp_i(ipartbin) = mpartbin_lims(ipartbin)
                                        endif

                                        if (frac_mp_i(ipartbin)>mpartbin_lims(ipartbin+1)) then
                                            print*,'too big i',ipartbin,frac_mp_i(ipartbin),mpartbin_lims(ipartbin:ipartbin+1),mpd(ix,iz,ipartbin,idropbin,3),nd(ix,iz,ipartbin,idropbin,3)
                                            stop
                                        else if (frac_mp_i(ipartbin)<mpartbin_lims(ipartbin)) then
                                            print*,'too small i',ipartbin,frac_mp_i(ipartbin),mpartbin_lims(ipartbin:ipartbin+1),mpd(ix,iz,ipartbin,idropbin,3),nd(ix,iz,ipartbin,idropbin,3)
                                            stop
                                        endif
                                
                                        frac_n_j(ipartbin) = nd(ix,iz,ipartbin,jdropbin,3)/nd_total(jdropbin)

                                        if ((mpd(ix,iz,ipartbin,jdropbin,3)>0.)  .and. (nd(ix,iz,ipartbin,idropbin,3)>0)) then
                                            frac_mp_j(ipartbin) = mpd(ix,iz,ipartbin,jdropbin,3)/nd(ix,iz,ipartbin,jdropbin,3)
                                        else
                                            frac_mp_j(ipartbin) = mpartbin_lims(ipartbin)
                                        endif

                                        if (frac_mp_j(ipartbin)>mpartbin_lims(ipartbin+1)) then
                                            print*,'too big j',ipartbin,frac_mp_j(ipartbin),mpartbin_lims(ipartbin:ipartbin+1),mpd(ix,iz,ipartbin,idropbin,3),nd(ix,iz,ipartbin,idropbin,3)
                                            stop
                                        else if (frac_mp_j(ipartbin)<mpartbin_lims(ipartbin)) then
                                            print*,'too small j',ipartbin,frac_mp_j(ipartbin),mpartbin_lims(ipartbin:ipartbin+1),mpd(ix,iz,ipartbin,idropbin,3),nd(ix,iz,ipartbin,idropbin,3)
                                            stop
                                        endif
                                        
                                        ! redistribute cloud number/mass based on collision-coalescence rate
                                        nd_final(ipartbin,idropbin) = nd_final(ipartbin,idropbin) - (Jij*d2t*frac_n_i(ipartbin))
                                        nd_final(ipartbin,jdropbin) = nd_final(ipartbin,idropbin) - (Jij*d2t*frac_n_j(ipartbin))

                                        ! all our droplets in a droplet bin are assumed to be the same size, but total size = particle + liquid water
                                        ! so the bigger the dry particle, the less liquid mass there is
   
                                        mld_final(ipartbin,idropbin) = mld_final(ipartbin,idropbin) - (Jij*d2t*frac_n_i(ipartbin)*(md_each_i-frac_mp_i(ipartbin)))
                                        mld_final(ipartbin,jdropbin) = mld_final(ipartbin,jdropbin) - (Jij*d2t*frac_n_j(ipartbin)*(md_each_j-frac_mp_j(ipartbin)))

                                        mpd_final(ipartbin,idropbin) = mpd_final(ipartbin,idropbin) - (Jij*d2t*frac_n_i(ipartbin)*frac_mp_i(ipartbin))
                                        mpd_final(ipartbin,jdropbin) = mpd_final(ipartbin,jdropbin) - (Jij*d2t*frac_n_j(ipartbin)*frac_mp_j(ipartbin))

   
                                    enddo 
                                    !print*,'subtract from old bins'

                                    ! find new droplet mass bin
                                    md_each_k = md_each_i + md_each_j
                                    kdropbin = COUNT(mdropbin_lims(:ndropbin)<=md_each_k) 
                            
                                    if (kdropbin>ndropbin) then
                                        print*,'collision makes too big drop',kdropbin,ndropbin,md_each_k
                                    else if ((kdropbin<idropbin) .or. (kdropbin<jdropbin)) then
                                        print*,'something is wrong, new bin is smaller than original bin drops',idropbin,jdropbin,kdropbin,md_each_i,md_each_j,md_each_k
                                    else
                                        ! check to make sure the new droplet has not become bigger than the biggest droplet bin we have
                                        if (md_each_k>mdropbin_lims(ndropbin+1)) then 
                                            print*,'too big droplet!'
                                        else
                                            do ipartbin=1,npartbin
                                                ! check that there are actually droplets in this bin
                                                if ((frac_n_i(ipartbin)>0.) .and. (frac_mp_i(ipartbin)>0.)) then
                                                    do jpartbin=1,npartbin 
                                                        ! check that there are actually droplets in this bin
                                                        if ((frac_n_j(jpartbin)>0.) .and. (frac_mp_j(jpartbin)>0.)) then
                                                    
                                                            mpd_each_k = frac_mp_i(ipartbin) + frac_mp_j(jpartbin)

                                                            ! the assingment to bin is actually correct,but smth is wrong when i calculate the initial frac_mp_i and frac_mp_j bins 
                                                            ! will need to make particle bin dimension
                                                            kpartbin = COUNT(mpartbin_lims(:npartbin)<=mpd_each_k)

                                                            if ((kpartbin < ipartbin) .or. (kpartbin < jpartbin)) then
                                                                print*,'something wrong, part bin smaller than original bins',ipartbin,jpartbin,kpartbin,frac_mp_i(ipartbin) ,frac_mp_j(jpartbin),mpd_each_k
                                                            endif
                                                            
                                                            if (mpd_each_k > md_each_k) then
                                                                print*,'more aerosol than water! smth is wrong',mpd_each_k,md_each_k
                                                            endif

                                                            nd_final(kpartbin,kdropbin) = nd_final(kpartbin,kdropbin) + (Jij*d2t*frac_n_i(ipartbin)*frac_n_j(jpartbin))
                                                            mld_final(kpartbin,kdropbin) = mld_final(kpartbin,kdropbin) + (Jij*d2t*frac_n_i(ipartbin)*frac_n_j(jpartbin)*(md_each_k-mpd_each_k))
                                                            mpd_final(kpartbin,kdropbin) = mpd_final(kpartbin,kdropbin) + (Jij*d2t*frac_n_i(ipartbin)*frac_n_j(jpartbin)*mpd_each_k)
                                                            
                                                            

                                                        endif 
                                                    enddo ! end j particle
                                                endif
                                            enddo ! end i particle
                                        endif
                                    endif
                                endif
                            endif 
                        enddo
                    endif 
                enddo

                !print*,'start negs check'

                ! move from temporary arrays to final
                do ipartbin=1,npartbin
                    do idropbin=1,ndropbin
                        if ((nd_final(ipartbin,idropbin)<=0.) .or. (nd_final(ipartbin,idropbin)<=0.) .or. (nd_final(ipartbin,idropbin)<=0.)) then
                            nd(ix,iz,ipartbin,idropbin,3) = 0.
                            mld(ix,iz,ipartbin,idropbin,3) = 0.
                            mpd(ix,iz,ipartbin,idropbin,3) = 0.
                        else
                            nd(ix,iz,ipartbin,idropbin,3)  = nd_final(ipartbin,idropbin)
                            mld(ix,iz,ipartbin,idropbin,3)  = mld_final(ipartbin,idropbin)
                            mpd(ix,iz,ipartbin,idropbin,3)  = mpd_final(ipartbin,idropbin)

                            
                        endif
                    enddo ! end droplet bin loop
                enddo ! end particle bin loop

                !print*,'values in main arrays now'
                
            enddo ! end x loop
        enddo ! end z loop
    end subroutine collisioncoalescence
    

    subroutine check_negs
        use run_constants, only: nx,nz,npartbin,ndropbin
        use model_vars, only: np, mp,nd,mld,mpd
        
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
                        if ((nd(ix,iz,ipb,idb,3)<0.) .or. (mld(ix,iz,ipb,idb,3)<minmass) .or. (mpd(ix,iz,ipb,idb,3)<minmass)) then
                            nd(ix,iz,ipb,idb,3)=0.
                            mld(ix,iz,ipb,idb,3)=0.
                            mpd(ix,iz,ipb,idb,3) = 0.
                        endif
                    enddo
                enddo
            enddo
        enddo

    end subroutine check_negs
end module microphysics 
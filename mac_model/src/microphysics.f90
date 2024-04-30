module microphysics
    implicit none

    contains
    
    subroutine microphysics_driv
        !this is the main microphysics driver
        use constants, only : p00, cp, rd
        use run_constants, only: nx, nz
        use model_vars, only: pib, pip, thb, thp, rvb, rvp,samb,t,p
        use thermo_functions, only: calc_rsat, calc_satfrac

        implicit none

        real :: th, pi, rv, rvsat
        integer :: ix, iz

        do ix=2,nx-1
            do iz=2,nz-1
                th = thp(ix,iz,3) + thb(iz)
                pi = pip(ix,iz,3) + pib(iz)
                t(ix,iz) = th*pi
                p(ix,iz) = p00 * (pi**(cp/rd))
                rv = rvp(ix,iz,3) + rvb(iz)
                rvsat = calc_rsat(t(ix,iz),p(ix,iz))
                Samb(ix,iz) = calc_satfrac(rv,rvsat)
            enddo
        enddo
    
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
        use model_vars, only: mdropbin_lims, mpartbin_lims  &
                        ,  np, mp, nd, mld, mpd, thp, rvp   &
                        ,  t, samb, pib
        use thermo_functions, only: calc_rsat,calc_satfrac
        use aerosol, only: calc_mass, calc_dp, calc_scrit

        implicit none
        
        integer :: ix,iz,ipartbin,idropbin ! counters
        real :: th, pi, rv, rvsat ! temporary variables for saturation calc
        real :: mp_each,dp_each,Sc_each ! temporary variables for activation calc
        real :: np_final(npartbin),mp_final(npartbin),nd_final(npartbin,ndropbin),mld_final(npartbin,ndropbin),mpd_final(npartbin,ndropbin)
        real :: temp1,temp2,temp3
        integer :: jpartbin,jdropbin
        !first loop over all spatial points
        do iz=2,nz-1
            do ix=2,nx-1
                ! initialize the output arrays with zero
                do ipartbin=1,npartbin
                    np_final(ipartbin) = np(ix,iz,ipartbin,3)
                    mp_final(ipartbin) = mp(ix,iz,ipartbin,3)

                    do idropbin=1,ndropbin
                        nd_final(ipartbin,idropbin)  = nd(ix,iz,ipartbin,idropbin,3)
                        mld_final(ipartbin,idropbin)  = mld(ix,iz,ipartbin,idropbin,3)
                        mpd_final(ipartbin,idropbin)  = mpd(ix,iz,ipartbin,idropbin,3)
                    enddo
                enddo
                
                ! if the grid point is supersaturated and
                if (samb(ix,iz)>1.) then
                    !print*,'ambient ss',Samb,rv,rvsat
                    ! and there are fewer cloud droplets than available CCN locally
                    ! this is what HUCM does -- I don't think it is that well physically-justified, but ok for now
                    if (SUM(nd(ix,iz,:,:,3)<SUM(np(ix,iz,:,3)))) then
                        ! then we should activate particles! 
                        do ipartbin=1,npartbin
                            ! for each particle bin, check if there are any particles at all (saves us some time looping over bins that have no aerosol)
                            if ((np(ix,iz,ipartbin,3) > 0.) .and. (mp(ix,iz,ipartbin,3) > 0.)) then
                                mp_each = mp(ix,iz,ipartbin,3)/np(ix,iz,ipartbin,3) ! mean mass of aerosol particle in bin
                                dp_each = calc_dp(mp_each, rhop) ! mean particle diameter in bin
                                Sc_each = calc_scrit(dp_each, t(ix,iz), kappa) !corresponding mean critical saturation in bin

                                if (samb(ix,iz)>=Sc_each) then !if ambient supersaturation is above critical, do activation
                                    
                                    ! determine which cloud bin the activate droplets should go into
                                    if ((2*mp_each) < mdropbin_lims(1)) then
                                        ! if the mean aerosol size is way smaller than the first bin, just shove them in first droplet bin
                                        idropbin = 1
                                    else
                                        ! else, find the appropriate bin where droplet size is just smaller than our activated particle
                                        idropbin = COUNT(mdropbin_lims(:ndropbin)<(2*mp_each))
                                    endif

                                    ! add liquid number & mass & aerosol mass
                                    ! we assume all the particles within a bin activate at the same time if they hit the mean critical supersat for the bin
                                    nd_final(ipartbin,idropbin) = nd_final(ipartbin,idropbin) + np(ix,iz,ipartbin,3) 
                                    mld_final(ipartbin,idropbin) = mld_final(ipartbin,idropbin) + ((0.5*(mdropbin_lims(idropbin)+mdropbin_lims(idropbin+1)))-mp_each)*np(ix,iz,ipartbin,3)
                                    mpd_final(ipartbin,idropbin) = mpd_final(ipartbin,idropbin) + (np(ix,iz,ipartbin,3)*mp_each)

                                    ! add the equivalent amount of latent heating to THP
                                    thp(ix,iz,3) = thp(ix,iz,3) + ((0.5*(mdropbin_lims(idropbin)+mdropbin_lims(idropbin+1)))*(lv/(cp*pib(iz))))
                                    ! remove equivalent amount of water from RVP
                                    rvp(ix,iz,3) = rvp(ix,iz,3) - ((0.5*(mdropbin_lims(idropbin)+mdropbin_lims(idropbin+1)))-mp_each)

                                    ! remove aerosol from unprocessed aerosol bins
                                    np_final(ipartbin) = np_final(ipartbin) - np(ix,iz,ipartbin,3)
                                    mp_final(ipartbin) = mp_final(ipartbin) - (np(ix,iz,ipartbin,3)*mp_each)
                                !else
                                !    print*,'no activate',ipartbin,Sc_each, Samb,np(ix,iz,ipartbin,3)
                                endif 
                            endif
                        enddo ! end part loop
                    endif
                endif 

                do ipartbin=1,npartbin
                    np(ix,iz,ipartbin,3) = np_final(ipartbin) 
                    mp(ix,iz,ipartbin,3) = mp_final(ipartbin) 

                    if (np(ix,iz,ipartbin,3)>0.) then 
                        temp1 = mp(ix,iz,ipartbin,3)/np(ix,iz,ipartbin,3) 

                        if ((temp1<mpartbin_lims(ipartbin)) .or. (temp1>mpartbin_lims(ipartbin+1))) then
                            jpartbin = COUNT(mpartbin_lims(:npartbin)<=temp1)
                            print*,'actv part size no match',ipartbin,jpartbin
                            mp(ix,iz,jpartbin,3) = mp(ix,iz,jpartbin,3) + mp(ix,iz,ipartbin,3)
                            mp(ix,iz,ipartbin,3) = 0.
                            np(ix,iz,jpartbin,3) = np(ix,iz,jpartbin,3) + np(ix,iz,ipartbin,3)
                            np(ix,iz,ipartbin,3) = 0.
                        endif
                    endif 

                    do idropbin=1,ndropbin
                        nd(ix,iz,ipartbin,idropbin,3) = nd_final(ipartbin,idropbin) 
                        mld(ix,iz,ipartbin,idropbin,3) = mld_final(ipartbin,idropbin)  
                        mpd(ix,iz,ipartbin,idropbin,3) = mpd_final(ipartbin,idropbin)  

                        if (nd(ix,iz,ipartbin,idropbin,3)>0.) then 

                            temp2 = mld(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3) 
                            if (temp2<0.) then
                                print*,temp2
                            endif 
                            temp3 = mpd(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3) 
                            if (temp3<0.) then
                                print*,temp3
                            endif 
                        

                            if (((temp2+temp3)<mdropbin_lims(idropbin)) .or. ((temp2+temp3)>mdropbin_lims(idropbin+1)) .or. (temp3<mpartbin_lims(ipartbin)) .or. (temp3>mpartbin_lims(ipartbin+1))) then
                                
                                jdropbin = COUNT(mdropbin_lims(:ndropbin)<=(temp2+temp3))
                                jpartbin = COUNT(mpartbin_lims(:npartbin)<=temp3)

                                print*,'actv drop size no match',ix,iz,idropbin,jdropbin,temp2,temp3,temp2+temp3,ipartbin,jpartbin,temp3
                        
                                mpd(ix,iz,jpartbin,jdropbin,3) = mpd(ix,iz,jpartbin,jdropbin,3) + mpd(ix,iz,ipartbin,idropbin,3)
                                mpd(ix,iz,ipartbin,idropbin,3) = 0.
                                
                                mld(ix,iz,jpartbin,jdropbin,3) = mld(ix,iz,jpartbin,jdropbin,3) + mld(ix,iz,ipartbin,idropbin,3)
                                mld(ix,iz,ipartbin,idropbin,3) = 0.
                        
                                nd(ix,iz,jpartbin,jdropbin,3) = nd(ix,iz,jpartbin,jdropbin,3) + nd(ix,iz,ipartbin,idropbin,3)
                                nd(ix,iz,ipartbin,idropbin,3) = 0.
                            endif
                    endif 

                    enddo
                enddo
            enddo
        enddo
    end subroutine activation

    subroutine condensation
        use constants, only: cp, rd, p00, lv,rhol,MW_w,sigma_w,trigpi,dg,kair,R
        use run_constants, only: nx,nz,npartbin,ndropbin,dt,d2t,kappa,rhop
        use model_vars, only: mdropbin_lims,mpartbin_lims, thp, rvp             &
                        , np, mp,nd,mld,mpd                                     &
                        , samb, t, p, pib
        use thermo_functions, only: calc_esat
        use aerosol, only: calc_mass, calc_dp, calc_lambd, calc_dahneke

        implicit none

        integer :: ix,iz,ipartbin,idropbin ! counters
        integer :: jdropbin ! index of new bin after condensation

        real :: md_each_i,mld_each_i,mpd_each_i ! per droplet in bin, mass total, liquid mass, and particle mass
        real :: md_each_j ! final mass per droplet after condensation
        real :: dpd_each_i,dpld_each_i,dppd_each_i,lambd, Seq, G, Im,A,beta, condensed_water
        real :: t1=0., t2=0., t3=0. !more temporary variables
        real :: minmass = 1.e-30

        real :: nd_final(npartbin,ndropbin),mld_final(npartbin,ndropbin),mpd_final(npartbin,ndropbin)
        real :: temp2,temp3
        integer :: jpartbin

        !first loop over all spatial points
        do iz=2,nz-1
            do ix=2,nx-1
                ! calculate some terms for Kohler curve later -- these only depend on temp and pressure, so do them per spatial point
                A=4*MW_w*sigma_w/(rd*t(ix,iz)*rhol)
                lambd = calc_lambd(p(ix,iz),t(ix,iz))

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
                                t1 = (R*t(ix,iz))/(calc_esat(t(ix,iz))*dg*MW_w)
                                t2 = (lv/(kair*t(ix,iz)))
                                t3 = (lv*MW_w/(R*t(ix,iz))-1)
                                if (t1<0.) then
                                    print*,t1,t2,t3
                                endif

                                G = t1 + (t2*t3)
                                
                                ! calculate dahneke parameter
                                beta = calc_dahneke(dpd_each_i, lambd, 1.)

                                ! amount of water condensed this timestep
                                Im = 2.*trigpi*dpd_each_i*beta*(samb(ix,iz)-Seq)*(1/G)

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

                do ipartbin=1,npartbin
                    do idropbin=1,ndropbin
                        nd(ix,iz,ipartbin,idropbin,3) = nd_final(ipartbin,idropbin) 
                        mld(ix,iz,ipartbin,idropbin,3) = mld_final(ipartbin,idropbin)  
                        mpd(ix,iz,ipartbin,idropbin,3) = mpd_final(ipartbin,idropbin)  

                        if (nd(ix,iz,ipartbin,idropbin,3)>0.) then 

                            temp2 = mld(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3) 
                            if (temp2<0.) then
                                print*,temp2
                            endif 
                            temp3 = mpd(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3) 
                            if (temp3<0.) then
                                print*,temp3
                            endif 
                            if (((temp2+temp3)<mdropbin_lims(idropbin)) .or. ((temp2+temp3)>mdropbin_lims(idropbin+1)) .or. (temp3<mpartbin_lims(ipartbin)) .or. (temp3>mpartbin_lims(ipartbin+1))) then
                                
                                jdropbin = COUNT(mdropbin_lims(:ndropbin)<=(temp2+temp3))
                                jpartbin = COUNT(mpartbin_lims(:npartbin)<=temp3)

                                print*,'cond drop size no match',idropbin,jdropbin,ipartbin,jpartbin
                        
                                mpd(ix,iz,jpartbin,jdropbin,3) = mpd(ix,iz,jpartbin,jdropbin,3) + mpd(ix,iz,ipartbin,idropbin,3)
                                mpd(ix,iz,ipartbin,idropbin,3) = 0.
                                
                                mld(ix,iz,jpartbin,jdropbin,3) = mld(ix,iz,jpartbin,jdropbin,3) + mld(ix,iz,ipartbin,idropbin,3)
                                mld(ix,iz,ipartbin,idropbin,3) = 0.
                        
                                nd(ix,iz,jpartbin,jdropbin,3) = nd(ix,iz,jpartbin,jdropbin,3) + nd(ix,iz,ipartbin,idropbin,3)
                                nd(ix,iz,ipartbin,idropbin,3) = 0.
                            endif
                        endif 

                    enddo
                enddo
                
        
            enddo ! end x loop
        enddo ! end z loop

    end subroutine condensation

    subroutine collisioncoalescence
        use constants, only: cp, rd, p00, lv,rhol,MW_w,sigma_w,trigpi,dg,kair,R
        use run_constants, only: nx,nz,npartbin,ndropbin,dt,d2t,rhop
        use model_vars, only: mdropbin_lims,mpartbin_lims,thp                   &
                        , np, mp,nd,mld,mpd                                     &
                        , samb, t, p, pib
        use thermo_functions, only: calc_rhoair
        use aerosol, only: calc_mass, calc_dp

        implicit none

        integer :: ix,iz,idropbin,jdropbin,kdropbin,ipartbin,jpartbin,kpartbin ! counters
        real :: rho_air
        real :: nd_tot_dbin(ndropbin),mld_tot_dbin(ndropbin), mpd_tot_dbin(ndropbin) ! total number in eachd droplet bin
        real :: nd_final(npartbin,ndropbin),mld_final(npartbin,ndropbin),mpd_final(npartbin,ndropbin) 
        real :: mld_each_i, mpd_each_i,md_each_i,mld_each_j,mpd_each_j,md_each_j
        real :: mpd_pbin_i(npartbin),mpd_pbin_j(npartbin)
        real :: mld_each_k, mpd_each_k, md_each_k
        real :: dpd_each_j,vol_each_i,vol_each_j,pj,Kij,Jij
        real :: fracn_dbin(npartbin,ndropbin)
        real :: temp2, temp3

        real :: minnum = 1.e-10, minmass = 1.e-30

        do iz = 2,nz-1
            do ix = 2,nx-1
                rho_air = calc_rhoair(t(ix,iz),p(ix,iz))

                do idropbin=1,ndropbin
                    nd_tot_dbin(idropbin) = 0.
                    mld_tot_dbin(idropbin) = 0.
                    mpd_tot_dbin(idropbin) = 0.
                enddo

                do idropbin=1,ndropbin 
                    do ipartbin=1,npartbin
                        ! calculate the total droplet number in each droplet bin
                        nd_tot_dbin(idropbin) = nd_tot_dbin(idropbin) + nd(ix,iz,ipartbin,idropbin,3)
                        mld_tot_dbin(idropbin) = mld_tot_dbin(idropbin) + mld(ix,iz,ipartbin,idropbin,3)
                        mpd_tot_dbin(idropbin) = mpd_tot_dbin(idropbin) + mpd(ix,iz,ipartbin,idropbin,3)

                        ! initialize the output variables temporarily
                        nd_final(ipartbin,idropbin) = nd(ix,iz,ipartbin,idropbin,3)
                        mld_final(ipartbin,idropbin) = mld(ix,iz,ipartbin,idropbin,3)
                        mpd_final(ipartbin,idropbin) = mpd(ix,iz,ipartbin,idropbin,3)
                       
                    enddo

                    do ipartbin=1,npartbin
                        !find how much number mixing ratio of each drop bin is in each particle bin
                        if (nd_tot_dbin(idropbin) > 1.e-2) then
                            fracn_dbin(ipartbin,idropbin) = nd(ix,iz,ipartbin,idropbin,3)/nd_tot_dbin(idropbin)
                        else
                            fracn_dbin(ipartbin,idropbin) = 0.
                        endif 
                    enddo
                enddo


                ! let's only calculate the collisions over droplets
                do idropbin=1,ndropbin 
                    if ((nd_tot_dbin(idropbin)>minnum) .and. ((mld_tot_dbin(idropbin)+mpd_tot_dbin(idropbin))>minmass)) then
                        mld_each_i = mld_tot_dbin(idropbin)/nd_tot_dbin(idropbin)

                        if (mld_each_i<0.) then
                            print*,mld_each_i !this doesn't do anything but stops an fp error
                        endif

                        mpd_each_i = mpd_tot_dbin(idropbin)/nd_tot_dbin(idropbin)
                        md_each_i = mld_each_i + mpd_each_i


                        !loop over larger dropbin
                        do jdropbin=idropbin,ndropbin
                            if ((nd_tot_dbin(jdropbin)>minnum) .and. ((mld_tot_dbin(jdropbin)+mpd_tot_dbin(jdropbin))>minmass)) then
                                mld_each_j = mld_tot_dbin(jdropbin)/nd_tot_dbin(jdropbin)
        
                                if (mld_each_j<0.) then
                                    print*,mld_each_j !this doesn't do anything but stops an fp error
                                endif
        
                                mpd_each_j = mpd_tot_dbin(jdropbin)/nd_tot_dbin(jdropbin)
                                md_each_j = mld_each_j + mpd_each_j
                                
                                dpd_each_j = ((calc_dp(mld_each_j,rhol)**3.) + (calc_dp(mpd_each_j,rhop)**3.))**(1/3.)

                                vol_each_i = ((mld_each_i/rhol) + (mpd_each_i/rhop)) * (100.**3.)
                                vol_each_j = ((mld_each_j/rhol) + (mpd_each_j/rhop)) * (100.**3.)

                                if (dpd_each_j>100.e-6) then
                                    pj = 5.78e3 * (vol_each_i+vol_each_j)
                                else
                                    pj = 9.44e9 * (vol_each_i**2 + vol_each_j**2)
                                endif

                                Kij = pj * rho_air * (100.**(-3.))
                                Jij = Kij * nd_tot_dbin(idropbin) * nd_tot_dbin(jdropbin)

                                do ipartbin=1,npartbin
                                    if (nd(ix,iz,ipartbin,idropbin,3)>0.) then
                                        mpd_pbin_i(ipartbin) = mpd(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3)
                                    else
                                        mpd_pbin_i(ipartbin) = mpartbin_lims(ipartbin)
                                    endif 

                                    if (nd(ix,iz,ipartbin,jdropbin,3)>0.) then
                                        mpd_pbin_j(ipartbin) = mpd(ix,iz,ipartbin,jdropbin,3)/nd(ix,iz,ipartbin,jdropbin,3)
                                    else
                                        mpd_pbin_j(ipartbin) = mpartbin_lims(ipartbin)
                                    endif 
                                    

                                    nd_final(ipartbin,idropbin) = nd_final(ipartbin,idropbin) - (Jij*d2t*fracn_dbin(ipartbin,idropbin))
                                    nd_final(ipartbin,jdropbin) = nd_final(ipartbin,jdropbin) - (Jij*d2t*fracn_dbin(ipartbin,jdropbin))

                                    mld_final(ipartbin,idropbin) = mld_final(ipartbin,idropbin) - (Jij*d2t*fracn_dbin(ipartbin,idropbin)*(md_each_i-mpd_pbin_i(ipartbin)))
                                    mld_final(ipartbin,jdropbin) = mld_final(ipartbin,jdropbin) - (Jij*d2t*fracn_dbin(ipartbin,jdropbin)*(md_each_j-mpd_pbin_j(ipartbin)))

                                    mpd_final(ipartbin,idropbin) = mpd_final(ipartbin,idropbin) - (Jij*d2t*fracn_dbin(ipartbin,idropbin)*mpd_pbin_i(ipartbin))
                                    mpd_final(ipartbin,jdropbin) = mpd_final(ipartbin,jdropbin) - (Jij*d2t*fracn_dbin(ipartbin,jdropbin)*mpd_pbin_j(ipartbin))
                                enddo ! end particle bins

                                !find new droplet bin
                                md_each_k = md_each_i + md_each_j
                                kdropbin = COUNT(mdropbin_lims(:ndropbin)<=md_each_k)

                                if (kdropbin>ndropbin) then
                                    print*,'collision makes too big drop',kdropbin,ndropbin,md_each_k
                                else if ((kdropbin<idropbin) .or. (kdropbin<jdropbin)) then
                                    print*,'something is wrong, new bin is smaller than original bin drops',idropbin,jdropbin,kdropbin,md_each_i,md_each_j,md_each_k
                                else
                                    do ipartbin=1,npartbin
                                        if (fracn_dbin(ipartbin,idropbin)>0.) then
                                            do jpartbin=1,npartbin
                                                if (fracn_dbin(jpartbin,jdropbin)>0.) then
                                                    mpd_each_k = mpd_pbin_i(ipartbin) + mpd_pbin_j(jpartbin)
                                                    kpartbin = COUNT(mpartbin_lims(:npartbin)<mpd_each_k)

                                                    if ((kpartbin < ipartbin) .or. (kpartbin < jpartbin)) then
                                                        print*,'something wrong, part bin smaller than original bins',ipartbin,jpartbin,kpartbin
                                                    endif
                                                    
                                                    if (mpd_each_k > md_each_k) then
                                                        print*,'more aerosol than water! smth is wrong',mpd_each_k,md_each_k
                                                    endif

                                                    nd_final(kpartbin,kdropbin) = nd_final(kpartbin,kdropbin) + (Jij * d2t * fracn_dbin(ipartbin,idropbin) * fracn_dbin(ipartbin,jdropbin))
                                                    mld_final(kpartbin,kdropbin) = mld_final(kpartbin,kdropbin) + (Jij * d2t * fracn_dbin(ipartbin,idropbin) * fracn_dbin(ipartbin,jdropbin) * (md_each_k-mpd_each_k))
                                                    mpd_final(kpartbin,kdropbin) = mpd_final(kpartbin,kdropbin) + (Jij * d2t * fracn_dbin(ipartbin,idropbin) * fracn_dbin(ipartbin,jdropbin) * mpd_each_k)
                                                endif
                                            enddo ! end j particle bin
                                        endif
                                    enddo ! end  i particle bins
                                endif
                            endif
                        enddo ! end j drop bin
                    endif !end check if there is anything in this bin
                enddo ! end i drop bin

                
                do idropbin=1,ndropbin
                    do ipartbin=1,npartbin
                        if ((nd_final(ipartbin,idropbin)<minnum) .and. (mld_final(ipartbin,idropbin)<minmass) .and. (mpd_final(ipartbin,idropbin)<minmass)) then
                            nd(ix,iz,ipartbin,idropbin,3) = 0.
                            mld(ix,iz,ipartbin,idropbin,3) = 0.
                            mpd(ix,iz,ipartbin,idropbin,3) = 0.
                        else
                            nd(ix,iz,ipartbin,idropbin,3) = nd_final(ipartbin,idropbin)
                            mld(ix,iz,ipartbin,idropbin,3) = mld_final(ipartbin,idropbin)
                            mpd(ix,iz,ipartbin,idropbin,3) = mpd_final(ipartbin,idropbin)
                        endif

                        if (nd(ix,iz,ipartbin,idropbin,3)>0.) then 

                            temp2 = mld(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3) 
                            if (temp2<0.) then
                                print*,temp2
                            endif 
                            temp3 = mpd(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3) 
                            if (temp3<0.) then
                                print*,temp3
                            endif 
                            if (((temp2+temp3)<mdropbin_lims(idropbin)) .or. ((temp2+temp3)>mdropbin_lims(idropbin+1)) .or. (temp3<mpartbin_lims(ipartbin)) .or. (temp3>mpartbin_lims(ipartbin+1))) then
                                
                                jdropbin = COUNT(mdropbin_lims(:ndropbin)<=(temp2+temp3))
                                if (jdropbin==0.) then
                                    jdropbin = 1
                                endif
                                jpartbin = COUNT(mpartbin_lims(:npartbin)<=temp3)

                                if (jpartbin==0.) then
                                    jpartbin = 1
                                endif

                                print*,'coll drop size no match',idropbin,jdropbin,ipartbin,jpartbin
                        
                                mpd(ix,iz,jpartbin,jdropbin,3) = mpd(ix,iz,jpartbin,jdropbin,3) + mpd(ix,iz,ipartbin,idropbin,3)
                                mpd(ix,iz,ipartbin,idropbin,3) = 0.
                                
                                mld(ix,iz,jpartbin,jdropbin,3) = mld(ix,iz,jpartbin,jdropbin,3) + mld(ix,iz,ipartbin,idropbin,3)
                                mld(ix,iz,ipartbin,idropbin,3) = 0.
                        
                                nd(ix,iz,jpartbin,jdropbin,3) = nd(ix,iz,jpartbin,jdropbin,3) + nd(ix,iz,ipartbin,idropbin,3)
                                nd(ix,iz,ipartbin,idropbin,3) = 0.
                            endif
                        endif 
                    enddo ! end part bin
                enddo ! end drop bin


            enddo ! end x loop
        enddo  ! end z loop




    end subroutine collisioncoalescence


    subroutine check_negs
        use run_constants, only: nx,nz,npartbin,ndropbin
        use model_vars, only: np,mp,nd,mld,mpd,mdropbin_lims,mpartbin_lims
        
        implicit none

        integer :: ix,iz,ipartbin,idropbin ! counters
        real:: minmass=1.e-30, minnum=1.e-10
        real :: mp_each_p, mp_each_d, ml_each_d
        integer :: jpartbin,jdropbin
        
        do ix=2,nx-1
            do iz=2,nz-1
                do ipartbin=1,npartbin
                    if ((np(ix,iz,ipartbin,3)<minnum) .or. (mp(ix,iz,ipartbin,3)<minmass)) then
                        np(ix,iz,ipartbin,3)=0.
                        mp(ix,iz,ipartbin,3)=0.
                    !else if (np(ix,iz,ipartbin,3)>0.) then
                        !print*,ix,iz,ipartbin,mp(ix,iz,ipartbin,3),np(ix,iz,ipartbin,3)
                    !    mp_each_p = mp(ix,iz,ipartbin,3)/np(ix,iz,ipartbin,3)
                        
                        !if the mass of each particle does not match with the bin, move it to the right bin
                    !    if ((mp_each_p<mpartbin_lims(ipartbin)) .or. (mp_each_p>mpartbin_lims(ipartbin+1))) then

                    !        print*,'part mass no match'
                    !        jpartbin = COUNT(mpartbin_lims(:npartbin)<=mp_each_p)

                    !        mp(ix,iz,jpartbin,3) = mp(ix,iz,jpartbin,3) + mp(ix,iz,ipartbin,3)
                    !        mp(ix,iz,ipartbin,3) = 0.
                    !        np(ix,iz,jpartbin,3) = np(ix,iz,jpartbin,3) + np(ix,iz,ipartbin,3)
                    !        np(ix,iz,ipartbin,3) = 0.
                    !    endif  
                    endif

                    do idropbin=1,ndropbin
                        if ((nd(ix,iz,ipartbin,idropbin,3)<minnum) .or. (mld(ix,iz,ipartbin,idropbin,3)<minmass) .or. (mpd(ix,iz,ipartbin,idropbin,3)<minmass)) then
                            nd(ix,iz,ipartbin,idropbin,3)=0.
                            mld(ix,iz,ipartbin,idropbin,3)=0.
                            mpd(ix,iz,ipartbin,idropbin,3) = 0.
                        !else if (nd(ix,iz,ipartbin,idropbin,3)>0.) then
                            !print*,mpd(ix,iz,ipartbin,idropbin,3),nd(ix,iz,ipartbin,idropbin,3)
                        !    mp_each_d = mpd(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3)
                        !    if (mp_each_d<0.) then
                        !        print*,mp_each_d
                        !    endif
                        !    ml_each_d = mld(ix,iz,ipartbin,idropbin,3)/nd(ix,iz,ipartbin,idropbin,3)
                           
                        !    if (ml_each_d<0.) then
                        !        print*,ml_each_d
                        !    endif
                        
                        !    if (((mp_each_d+ml_each_d)<mdropbin_lims(idropbin)) .or. ((mp_each_d+ml_each_d)>mdropbin_lims(idropbin+1))) then
                        !        print*,'droplet mass no match'
                        !        jdropbin = COUNT(mdropbin_lims(:ndropbin)<=(mp_each_d+ml_each_d))
                        !        jpartbin = COUNT(mpartbin_lims(:npartbin)<=mp_each_d)
                        
                        !        mpd(ix,iz,jpartbin,jdropbin,3) = mpd(ix,iz,jpartbin,jdropbin,3) + mpd(ix,iz,ipartbin,idropbin,3)
                        !        mpd(ix,iz,ipartbin,idropbin,3) = 0.
                                
                        !        mld(ix,iz,jpartbin,jdropbin,3) = mld(ix,iz,jpartbin,jdropbin,3) + mld(ix,iz,ipartbin,idropbin,3)
                        !        mld(ix,iz,ipartbin,idropbin,3) = 0.
                        
                        !        nd(ix,iz,jpartbin,jdropbin,3) = nd(ix,iz,jpartbin,jdropbin,3) + nd(ix,iz,ipartbin,idropbin,3)
                        !        nd(ix,iz,ipartbin,idropbin,3) = 0.
                        !    endif
                        endif
                    enddo
                enddo
            enddo
        enddo

        ! also check that they are in the right bins depending on mld/nd, mpd/nd, etc.
        ! (mld+mpd)/nd should be in mdbin_lims(idropbin:idropbin+1) and mpd/nd should be in mpdbins(ipartbin:ipartbin+1)

    end subroutine check_negs
end module microphysics 
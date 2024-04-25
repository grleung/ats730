module solve_prog
    !This module contains the full tendency equations for u,w,theta_prime,PI_prime from HW4

    implicit none

    contains
    
    subroutine tendencies
        implicit none

        call zero_tends

        call calc_u_tend
        print*,'u tend'
        call calc_w_tend
        print*,'w tend'
        call calc_pip_tend
        print*,'pip tend'
        call calc_thp_tend
        print*,'thp tend'
        call calc_rvp_tend
        print*,'rvp tend'
        call calc_aero_tend
        print*,'aero tend'
        call calc_liquid_tend
        print*,'liquid tend'

        call apply_tends

        
        !add advection/diffusion for aerosol,cloud,rain

    end subroutine tendencies

    subroutine zero_tends
        use run_constants, only: nz,nx,npartbin,ndropbin
        use model_vars, only:u_xadv,u_zadv,u_pgf,u_xdiff,u_zdiff,u_tend_total                        &
                            ,w_xadv,w_zadv,w_pgf,w_buoy,w_xdiff,w_zdiff,w_tend_total                 &
                            ,thp_xadv,thp_zadv,thp_meanadv,thp_xdiff,thp_zdiff,thp_tend_total    &
                            ,pip_xadv,pip_zadv,pip_xdiff,pip_zdiff,pip_tend_total                &
                            ,rvp_xadv,rvp_zadv,rvp_meanadv,rvp_xdiff,rvp_zdiff,rvp_tend_total    &
                            ,np_tend_total,mp_tend_total,nd_tend_total,mld_tend_total,mpd_tend_total  

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate
        integer :: iab, idb

        do iz = 1, nz
            do ix = 1, nx
                u_xadv(ix,iz)=0.
                u_zadv(ix,iz)=0.
                u_pgf(ix,iz)=0.
                u_xdiff(ix,iz)=0.
                u_zdiff(ix,iz)=0.
                u_tend_total(ix,iz)=0.
                w_xadv(ix,iz)=0.
                w_zadv(ix,iz)=0.
                w_pgf(ix,iz)=0.
                w_buoy(ix,iz)=0.
                w_xdiff(ix,iz)=0.
                w_zdiff(ix,iz)=0.
                w_tend_total(ix,iz)=0.
                thp_xadv(ix,iz)=0.
                thp_zadv(ix,iz)=0.
                thp_meanadv(ix,iz)=0.
                thp_xdiff(ix,iz)=0.
                thp_zdiff(ix,iz)=0.
                thp_tend_total(ix,iz)=0.
                pip_xadv(ix,iz)=0.
                pip_zadv(ix,iz)=0.
                pip_xdiff(ix,iz)=0.
                pip_zdiff(ix,iz)=0.
                pip_tend_total(ix,iz)=0.
                rvp_xadv(ix,iz)=0.
                rvp_zadv(ix,iz)=0.
                rvp_meanadv(ix,iz)=0.
                rvp_xdiff(ix,iz)=0.
                rvp_zdiff(ix,iz)=0.
                rvp_tend_total(ix,iz)=0.

                do iab=1,npartbin
                    np_tend_total(ix,iz,iab) = 0.
                    mp_tend_total(ix,iz,iab) = 0.

                    do idb=1,ndropbin
                        nd_tend_total(ix,iz,iab,idb) = 0.
                        mld_tend_total(ix,iz,iab,idb) = 0.
                        mpd_tend_total(ix,iz,iab,idb) = 0.
                    enddo
                enddo
                
            enddo
        enddo 
    end subroutine zero_tends

    subroutine calc_u_tend
        use run_constants, only: nz,nx,kmx,kmz,rdx,rdz,d2t
        use constants, only: cp
        use model_vars, only:thvb,rhoub,rhowb,pip,up,wp &
                            ,u_xadv,u_zadv,u_pgf,u_xdiff,u_zdiff,u_tend_total    

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate

        
        ! calculate tendency in u
        ! this is equation 6 in HW4

        ! loop over real/unique points
        
        do ix = 2, nx-1
            do iz = 2, nz-1
                ! first term in u-tendency equation: horizontal advection term = -d(uu)/dx
                u_xadv(ix,iz) = - rdx * (((0.5*(up(ix+1,iz,2)+up(ix,iz,2)))**2) &
                                        - ((0.5*(up(ix,iz,2)+up(ix-1,iz,2)))**2))

                ! second term in u-tendency equation: vertical advection term = -1/rho * d(rho*u*w)/dz
                u_zadv(ix,iz) = - (rdz/(rhoub(iz)))                                                 &
                                    * 0.25 * ((rhowb(iz+1) *    (wp(ix-1,iz+1,2)    +   wp(ix,iz+1,2))  * (up(ix,iz,2)    +   up(ix,iz+1,2))) &
                                    -         (rhowb(iz)   *    (wp(ix-1,iz,2)      +   wp(ix,iz,2))    * (up(ix,iz-1,2)    +   up(ix,iz,2))))

                ! third term in u-tendency equation: pressure gradient term = -cp * theta_base * d(pi_pert)/ dx
                u_pgf(ix,iz) = -cp*thvb(iz)*rdx                                                   & 
                                    * (pip(ix,iz,2)-pip(ix-1,iz,2))

                ! fourth term in u-tendency equation: horizontal diffusion
                u_xdiff(ix,iz) = kmx * rdx * rdx * (up(ix-1,iz,1) - (2*up(ix,iz,2)) + up(ix+1,iz,1))

                ! fifth term in u-tendency equation: vertical diffusion
                u_zdiff(ix,iz) = kmz * rdz * rdz * (up(ix,iz-1,1) - (2*up(ix,iz,1)) + up(ix,iz+1,1))
                
                u_tend_total(ix,iz) = u_xadv(ix,iz) + u_zadv(ix,iz) + u_pgf(ix,iz) + u_xdiff(ix,iz) + u_zdiff(ix,iz)
            enddo ! end x loop
        enddo ! end z loop
    end subroutine calc_u_tend

    subroutine calc_w_tend
        use run_constants, only: nz,nx,kmx,kmz,rdx,rdz,d2t
        use constants, only: cp,g
        use model_vars, only:thvb,thb,rhoub,rhowb,pip,up,wp,thp,rvp,mld &
                            ,w_xadv,w_zadv,w_pgf,w_buoy,w_xdiff,w_zdiff,w_tend_total    

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate tendency in w
        ! this is equation 7 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iz = 3, nz-1
                ! first term in w-tendency equation: horizontal advection term 
                w_xadv(ix,iz) = - rdx                                                                              &
                                    * 0.25 * (((up(ix+1,iz,2)+up(ix+1,iz-1,2))*(wp(ix+1,iz,2)+wp(ix,iz,2)))         &                                        
                                             -((up(ix,iz,2)+up(ix,iz-1,2))*(wp(ix,iz,2)+wp(ix-1,iz,2))) )

                ! second term in w-tendency equation: vertical advection term 
                w_zadv(ix,iz) =  - rdz/rhoub(iz) &
                                        * ((rhowb(iz+1)*(0.5*(wp(ix,iz+1,2)+wp(ix,iz,2)))**2) &
                                        - (rhowb(iz)*(0.5*(wp(ix,iz,2)+wp(ix,iz-1,2)))**2))

                ! third term in w-tendency equation: pressure gradient term 
                w_pgf(ix,iz) = -cp * rdz * 0.25 * (thvb(iz)+thvb(iz-1)) * (pip(ix,iz,2)-pip(ix,iz-1,2))

                ! fourth term in w-tendency equation: pressure gradient term 
                w_buoy(ix,iz) = g * ((thp(ix,iz,2)+thp(ix,iz-1,2))/(thb(iz)+thb(iz-1))                  &
                                    + 0.61*0.5*(rvp(ix,iz,2)+rvp(ix,iz-1,2))                            &
                                    - 0.5*(SUM(mld(ix,iz,:,:,2))+SUM(mld(ix,iz-1,:,:,2)))     &
                                    )

                ! fifth term in w-tendency equation: horizontal diffusion
                w_xdiff(ix,iz) = kmx * rdx * rdx * (wp(ix-1,iz,1) - (2*wp(ix,iz,1)) + wp(ix+1,iz,1))

                ! sixth term in w-tendency equation: vertical diffusion
                w_zdiff(ix,iz) = kmz * rdz * rdz * (wp(ix,iz-1,1) - (2*wp(ix,iz,1)) + wp(ix,iz+1,1))
                
                w_tend_total(ix,iz) = w_xadv(ix,iz) +w_zadv(ix,iz) +w_pgf(ix,iz) + w_buoy(ix,iz) + w_xdiff(ix,iz) + w_zdiff(ix,iz)
            enddo ! end x loop
        enddo ! end z loop

    end subroutine calc_w_tend

    subroutine calc_thp_tend
        use run_constants, only: nz,nx,khx,khz,rdx,rdz,d2t
        use constants, only: cp
        use model_vars, only:thb,rhoub,rhowb,up,wp,thp &
                            ,thp_xadv,thp_zadv,thp_meanadv,thp_pgf,thp_xdiff,thp_zdiff,thp_tend_total  

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate tendency in theta
        ! this is equation 8 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iz = 2, nz-1
                ! first term in thp-tendency equation: horizontal advection term = -d(u * thp)/dx
                thp_xadv(ix,iz) = - rdx                                                                   &
                                    * 0.5 * ((up(ix+1,iz,2) * (thp(ix+1,iz,2) +  thp(ix,iz,2)))     &
                                    -        (up(ix,iz,2)   * (thp(ix,iz,2)   +  thp(ix-1,iz,2))))

                ! second term in thp-tendency equation: vertical advection term = -1/rho * d(rho*w * thp)/dz
                thp_zadv(ix,iz) = - rdz/(rhoub(iz))                                                 &
                                    * 0.5 *  ((rhowb(iz+1) *    wp(ix,iz+1,2)     * (thp(ix,iz+1,2)    +   thp(ix,iz,2))) &
                                    -         (rhowb(iz)   *    wp(ix,iz,2)     * (thp(ix,iz,2)    +   thp(ix,iz-1,2))))

                ! term in thp-tendency equation: mean state advection
                thp_meanadv(ix,iz) = -0.5 * (rdz/rhoub(iz))                                           &
                                                * ((rhowb(iz) * wp(ix,iz,2) * (thb(iz)-thb(iz-1)))      &
                                                +  (rhowb(iz+1) * wp(ix,iz+1,2) * (thb(iz+1)-thb(iz))))

                ! fourth term in thp-tendency equation: horizontal diffusion
                thp_xdiff(ix,iz) = khx * rdx * rdx * (thp(ix-1,iz,1) - (2*thp(ix,iz,1)) + thp(ix+1,iz,1))

                ! fifth term in thp-tendency equation: vertical diffusion
                thp_zdiff(ix,iz) = khz * rdz * rdz * (thp(ix,iz-1,1) - (2*thp(ix,iz,1)) + thp(ix,iz+1,1))

                thp_tend_total(ix,iz) = thp_xadv(ix,iz) + thp_zadv(ix,iz) + thp_meanadv(ix,iz) + thp_xdiff(ix,iz) + thp_zdiff(ix,iz)
            enddo ! end x loop
        enddo ! end z loop

    end subroutine calc_thp_tend

    subroutine calc_pip_tend
        use run_constants, only: nz,nx,khx,khz,rdx,rdz,cs
        use constants, only: cp
        use model_vars, only:thvb,rhoub,rhowb,up,wp,pip &
                     ,pip_xadv,pip_zadv,pip_xdiff,pip_zdiff,pip_tend_total          


        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate tendency in perturbation exner function
        ! this is equation 9 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iz = 2, nz-1
                ! first term in pip-tendency equation: horizontal advection term
                pip_xadv(ix,iz) = -((cs**2)*rdx/(cp*thvb(iz))) * (up(ix+1,iz,2)-up(ix,iz,2))

                ! second term in pip-tendency equation: vertical advection term 
                pip_zadv(ix,iz) = -(((cs**2)*rdz*0.5/(rhoub(iz)*cp*(thvb(iz)**2))))   &
                                *  ((rhowb(iz+1)*wp(ix,iz+1,2)*(thvb(iz+1)+thvb(iz))) &
                                  -(rhowb(iz)*wp(ix,iz,2)*(thvb(iz)+thvb(iz-1))))

                ! third term in thp-tendency equation: horizontal diffusion
                pip_xdiff(ix,iz) = khx * rdx * rdx * (pip(ix-1,iz,1) - (2*pip(ix,iz,1)) + pip(ix+1,iz,1))

                ! fourth term in thp-tendency equation: vertical diffusion
                pip_zdiff(ix,iz) = khz * rdz * rdz * (pip(ix,iz-1,1) - (2*pip(ix,iz,1)) + pip(ix,iz+1,1))

                pip_tend_total(ix,iz) = pip_xadv(ix,iz)+pip_zadv(ix,iz)+pip_xdiff(ix,iz)+pip_zdiff(ix,iz)
            enddo ! end x loop
        enddo

        
    end subroutine calc_pip_tend

    subroutine calc_rvp_tend
        use run_constants, only: nz,nx,khx,khz,rdx,rdz,cs
        use constants, only: cp
        use model_vars, only:rvb,rhoub,rhowb,up,wp,rvp &
                        ,rvp_xadv,rvp_zadv,rvp_meanadv,rvp_xdiff,rvp_zdiff,rvp_tend_total

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate tendency in rvp
        ! loop over real/unique points
        do ix = 2, nx-1
            do iz = 2, nz-1
                rvp_xadv(ix,iz) = - rdx                                                                   &
                                    * 0.5 * ((up(ix+1,iz,2) * (rvp(ix+1,iz,2) +  rvp(ix,iz,2)))     &
                                    -        (up(ix,iz,2)   * (rvp(ix,iz,2)   +  rvp(ix-1,iz,2))))

                rvp_zadv(ix,iz) = - rdz/(rhoub(iz))                                                 &
                                    * 0.5 *  ((rhowb(iz+1) *    wp(ix,iz+1,2)     * (rvp(ix,iz+1,2)    +   rvp(ix,iz,2))) &
                                    -         (rhowb(iz)   *    wp(ix,iz,2)     * (rvp(ix,iz,2)    +   rvp(ix,iz-1,2))))

                rvp_meanadv(ix,iz) = -0.5 * (1/rhoub(iz)) *   rdz                                     &
                                                * ((rhowb(iz) * wp(ix,iz,2) * (rvb(iz)-rvb(iz-1)))      &
                                                +  (rhowb(iz+1) * wp(ix,iz+1,2) * (rvb(iz+1)-rvb(iz))))

                rvp_xdiff(ix,iz) = khx * rdx * rdx * (rvp(ix-1,iz,1) - (2*rvp(ix,iz,1)) + rvp(ix+1,iz,1))

                rvp_zdiff(ix,iz) = khz * rdz * rdz * (rvp(ix,iz-1,1) - (2*rvp(ix,iz,1)) + rvp(ix,iz+1,1))

                rvp_tend_total(ix,iz) = rvp_xadv(ix,iz) + rvp_zadv(ix,iz) + rvp_meanadv(ix,iz) + rvp_xdiff(ix,iz) + rvp_zdiff(ix,iz)
            enddo ! end x loop
        enddo ! end z loop

    end subroutine calc_rvp_tend

    subroutine calc_aero_tend
        use run_constants, only: nz,nx,npartbin,khx,khz,rdx,rdz,cs
        use constants, only: cp
        use model_vars, only:rvb,rhoub,rhowb,up,wp,np,mp &
                        ,np_tend_total,mp_tend_total

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate
        integer :: ipb ! counter for particle bins

        real :: mp_each


        do ix = 2, nx-1
            do iz = 2, nz-1
                do ipb=1,npartbin

                    mp_tend_total(ix,iz,ipb) = -rdx * 0.5 * ((up(ix+1,iz,2)*(mp(ix+1,iz,ipb,2)+mp(ix,iz,ipb,2)))         &
                                                            -(up(ix,iz,2)*(mp(ix,iz,ipb,2)+mp(ix-1,iz,ipb,2)))) &
                                            - (rdz/rhoub(iz)) * 0.5 * ((rhowb(iz+1)*wp(ix,iz+1,2)*(mp(ix,iz+1,ipb,2)+mp(ix,iz,ipb,2))) &
                                                                        - (rhowb(iz)*wp(ix,iz,2)*(mp(ix,iz,ipb,2)+mp(ix,iz-1,ipb,2)))) &
                                            + (khx * rdx * rdx * (mp(ix-1,iz,ipb,1) - 2*mp(ix,iz,ipb,1) + mp(ix+1,iz,ipb,1))) &
                                            + (khz * rdz * rdz * (mp(ix,iz-1,ipb,1) - 2*mp(ix,iz,ipb,1) + mp(ix,iz+1,ipb,1))) 


                    np_tend_total(ix,iz,ipb) = -rdx * 0.5 * ((up(ix+1,iz,2)*(np(ix+1,iz,ipb,2)+np(ix,iz,ipb,2)))         &
                                                            -(up(ix,iz,2)*(np(ix,iz,ipb,2)+np(ix-1,iz,ipb,2)))) &
                                            - (rdz/rhoub(iz)) * 0.5 * ((rhowb(iz+1)*wp(ix,iz+1,2)*(np(ix,iz+1,ipb,2)+np(ix,iz,ipb,2))) &
                                                                        - (rhowb(iz)*wp(ix,iz,2)*(np(ix,iz,ipb,2)+np(ix,iz-1,ipb,2)))) &
                                            + (khx * rdx * rdx * (np(ix-1,iz,ipb,1) - 2*np(ix,iz,ipb,1) + np(ix+1,iz,ipb,1))) &
                                            + (khz * rdz * rdz * (np(ix,iz-1,ipb,1) - 2*np(ix,iz,ipb,1) + np(ix,iz+1,ipb,1))) 

                    
                enddo
            enddo ! end x loop
        enddo ! end z loop

    end subroutine calc_aero_tend

    

    subroutine calc_liquid_tend
        use run_constants, only: nz,nx,npartbin,ndropbin,khx,khz,rdx,rdz,cs
        use constants, only: cp
        use model_vars, only:rvb,rhoub,rhowb,up,wp,nd,mld,mpd &
                        ,nd_tend_total,mld_tend_total,mpd_tend_total

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate
        integer :: ipb ! counter for particle bins
        integer :: idb !counter for droplet bins

        real :: mp_each


        do ix = 2, nx-1
            do iz = 2, nz-1
                do ipb=1,npartbin
                    do idb=1,ndropbin
                        mld_tend_total(ix,iz,ipb,idb) = -rdx * 0.5 * ((up(ix+1,iz,2)*(mld(ix+1,iz,ipb,idb,2)+mld(ix,iz,ipb,idb,2)))         &
                                                                -(up(ix,iz,2)*(mld(ix,iz,ipb,idb,2)+mld(ix-1,iz,ipb,idb,2)))) &
                                                - (rdz/rhoub(iz)) * 0.5 * ((rhowb(iz+1)*wp(ix,iz+1,2)*(mld(ix,iz+1,ipb,idb,2)+mld(ix,iz,ipb,idb,2))) &
                                                                            - (rhowb(iz)*wp(ix,iz,2)*(mld(ix,iz,ipb,idb,2)+mld(ix,iz-1,ipb,idb,2)))) &
                                                + (khx * rdx * rdx * (mld(ix-1,iz,ipb,idb,1) - 2*mld(ix,iz,ipb,idb,1) + mld(ix+1,iz,ipb,idb,1))) &
                                                + (khz * rdz * rdz * (mld(ix,iz-1,ipb,idb,1) - 2*mld(ix,iz,ipb,idb,1) + mld(ix,iz+1,ipb,idb,1))) 

                        nd_tend_total(ix,iz,ipb,idb) = -rdx * 0.5 * ((up(ix+1,iz,2)*(nd(ix+1,iz,ipb,idb,2)+nd(ix,iz,ipb,idb,2)))         &
                                                                -(up(ix,iz,2)*(nd(ix,iz,ipb,idb,2)+nd(ix-1,iz,ipb,idb,2)))) &
                                                - (rdz/rhoub(iz)) * 0.5 * ((rhowb(iz+1)*wp(ix,iz+1,2)*(nd(ix,iz+1,ipb,idb,2)+nd(ix,iz,ipb,idb,2))) &
                                                                            - (rhowb(iz)*wp(ix,iz,2)*(nd(ix,iz,ipb,idb,2)+nd(ix,iz-1,ipb,idb,2)))) &
                                                + (khx * rdx * rdx * (nd(ix-1,iz,ipb,idb,1) - 2*nd(ix,iz,ipb,idb,1) + nd(ix+1,iz,ipb,idb,1))) &
                                                + (khz * rdz * rdz * (nd(ix,iz-1,ipb,idb,1) - 2*nd(ix,iz,ipb,idb,1) + nd(ix,iz+1,ipb,idb,1))) 

                        mpd_tend_total(ix,iz,ipb,idb) = -rdx * 0.5 * ((up(ix+1,iz,2)*(mpd(ix+1,iz,ipb,idb,2)+mpd(ix,iz,ipb,idb,2)))         &
                                                                -(up(ix,iz,2)*(mpd(ix,iz,ipb,idb,2)+mpd(ix-1,iz,ipb,idb,2)))) &
                                                - (rdz/rhoub(iz)) * 0.5 * ((rhowb(iz+1)*wp(ix,iz+1,2)*(mpd(ix,iz+1,ipb,idb,2)+mpd(ix,iz,ipb,idb,2))) &
                                                                            - (rhowb(iz)*wp(ix,iz,2)*(mpd(ix,iz,ipb,idb,2)+mpd(ix,iz-1,ipb,idb,2)))) &
                                                + (khx * rdx * rdx * (mpd(ix-1,iz,ipb,idb,1) - 2*mpd(ix,iz,ipb,idb,1) + mpd(ix+1,iz,ipb,idb,1))) &
                                                + (khz * rdz * rdz * (mpd(ix,iz-1,ipb,idb,1) - 2*mpd(ix,iz,ipb,idb,1) + mpd(ix,iz+1,ipb,idb,1))) 

                    enddo ! end drop loop
                enddo ! end particle loop
            enddo ! end x loop
        enddo ! end z loop

    end subroutine calc_liquid_tend

    


    
    subroutine apply_tends
        use run_constants, only: nz,nx,dt,npartbin,ndropbin
        use model_vars, only:it,thp,pip,up,wp,rvp,np,mp,nd,mld,mpd                                &
                            ,u_tend_total,w_tend_total,thp_tend_total,pip_tend_total,rvp_tend_total &
                            ,np_tend_total,mp_tend_total,nd_tend_total,mld_tend_total        &
                            ,mpd_tend_total

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate
        integer :: iab,idb

        real :: rdx,rdz,d2t
        
        d2t       = (dt+dt)         ! 2*dt
        
        ! calculate actual future values
        ! loop over real/unique points
        do ix = 2, nx-1
            do iz = 2, nz-1
                if (it==1) then !for very frist timestep
                    up(ix,iz,3) = up(ix,iz,2) + (dt * u_tend_total(ix,iz))
                    wp(ix,iz,3) = wp(ix,iz,2) + (dt * w_tend_total(ix,iz))
                    thp(ix,iz,3) = thp(ix,iz,2) + (dt * thp_tend_total(ix,iz))
                    pip(ix,iz,3) = pip(ix,iz,2) + (dt * pip_tend_total(ix,iz))
                    rvp(ix,iz,3) = rvp(ix,iz,2) + (dt * rvp_tend_total(ix,iz))

                    do iab=1,npartbin
                        np(ix,iz,iab,3) = np(ix,iz,iab,2) + (dt * np_tend_total(ix,iz,iab)) 
                        mp(ix,iz,iab,3) = mp(ix,iz,iab,2) + (dt * mp_tend_total(ix,iz,iab))

                        do idb=1,ndropbin
                            nd(ix,iz,iab,idb,3) = nd(ix,iz,iab,idb,2) + (dt * nd_tend_total(ix,iz,iab,idb)) 
                            mld(ix,iz,iab,idb,3) = mld(ix,iz,iab,idb,2) + (dt * mld_tend_total(ix,iz,iab,idb))                        
                            mpd(ix,iz,iab,idb,3) = mpd(ix,iz,iab,idb,2) + (dt * mpd_tend_total(ix,iz,iab,idb))
                            
                        enddo
                    enddo 
                else
                    up(ix,iz,3) = up(ix,iz,1) + (d2t * u_tend_total(ix,iz))
                    wp(ix,iz,3) = wp(ix,iz,1) + (d2t * w_tend_total(ix,iz))
                    thp(ix,iz,3) = thp(ix,iz,1) + (d2t * thp_tend_total(ix,iz))
                    pip(ix,iz,3) = pip(ix,iz,1) + (d2t * pip_tend_total(ix,iz))
                    rvp(ix,iz,3) = rvp(ix,iz,1) + (d2t * rvp_tend_total(ix,iz))

                    do iab=1,npartbin
                        np(ix,iz,iab,3) = np(ix,iz,iab,1) + (d2t * np_tend_total(ix,iz,iab)) 
                        mp(ix,iz,iab,3) = mp(ix,iz,iab,1) + (d2t * mp_tend_total(ix,iz,iab))

                        do idb=1,ndropbin
                            nd(ix,iz,iab,idb,3) = nd(ix,iz,iab,idb,1) + (d2t * nd_tend_total(ix,iz,iab,idb)) 
                            mld(ix,iz,iab,idb,3) = mld(ix,iz,iab,idb,1) + (d2t * mld_tend_total(ix,iz,iab,idb))
                            mpd(ix,iz,iab,idb,3) = mpd(ix,iz,iab,idb,1) + (d2t * mpd_tend_total(ix,iz,iab,idb))
                        enddo
                    enddo 
                endif
            enddo ! end x loop
        enddo ! end z loop
    end subroutine apply_tends

end module solve_prog
module solve_prog
    !This module contains the full tendency equations for u,w,theta_prime,PI_prime, RVP_prime, RCP, RRP

    implicit none

    contains
    
    subroutine tendencies
        implicit none

        ! first zero out tendencies 
        call zero_tends

        ! calculate each of the tendency terms
        call calc_tend_u
        call calc_tend_v
        call calc_tend_w
        call calc_tend_thp
        call calc_tend_pip
        call calc_tend_rvp
        call calc_tend_rcp
        call calc_tend_rrp

        ! do actual time integration to apply tendencies we calculated above
        call apply_tends

    end subroutine tendencies

    subroutine zero_tends
        use run_constants, only: nz,nx,ny
        use model_vars, only:vap2cld,rain2vap,cld2rain_accr,cld2rain_auto                                           &
                            ,u_xadv,u_yadv,u_zadv,u_pgf,u_xdiff,u_ydiff,u_zdiff,u_tend_total                        &
                            ,v_xadv,v_yadv,v_zadv,v_pgf,v_xdiff,v_ydiff,v_zdiff,v_tend_total                        &
                            ,w_xadv,w_yadv,w_zadv,w_pgf,w_buoy,w_xdiff,w_ydiff,w_zdiff,w_tend_total                 &
                            ,thp_xadv,thp_yadv,thp_zadv,thp_meanadv,thp_xdiff,thp_ydiff,thp_zdiff,thp_tend_total    &
                            ,pip_xadv,pip_yadv,pip_zadv,pip_xdiff,pip_ydiff,pip_zdiff,pip_tend_total                &
                            ,rvp_xadv,rvp_yadv,rvp_zadv,rvp_meanadv,rvp_xdiff,rvp_ydiff,rvp_zdiff,rvp_tend_total    &
                            ,rcp_xadv,rcp_yadv,rcp_zadv,rcp_xdiff,rcp_ydiff,rcp_zdiff,rcp_tend_total                &
                            ,rrp_xadv,rrp_yadv,rrp_zadv,rrp_xdiff,rrp_ydiff,rrp_zdiff,rrp_tend_total   

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        ! first, reset all tendencies to zero
        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    vap2cld(ix,iy,iz)=0.
                    rain2vap(ix,iy,iz)=0.
                    cld2rain_accr(ix,iy,iz)=0.
                    cld2rain_auto(ix,iy,iz)=0.
                    u_xadv(ix,iy,iz)=0.
                    u_yadv(ix,iy,iz)=0.
                    u_zadv(ix,iy,iz)=0.
                    u_pgf(ix,iy,iz)=0.
                    u_xdiff(ix,iy,iz)=0.
                    u_ydiff(ix,iy,iz)=0.
                    u_zdiff(ix,iy,iz)=0.
                    u_tend_total(ix,iy,iz)=0.
                    v_xadv(ix,iy,iz)=0.
                    v_yadv(ix,iy,iz)=0.
                    v_zadv(ix,iy,iz)=0.
                    v_pgf(ix,iy,iz)=0.
                    v_xdiff(ix,iy,iz)=0.
                    v_ydiff(ix,iy,iz)=0.
                    v_zdiff(ix,iy,iz)=0.
                    v_tend_total(ix,iy,iz)=0.
                    w_xadv(ix,iy,iz)=0.
                    w_yadv(ix,iy,iz)=0.
                    w_zadv(ix,iy,iz)=0.
                    w_pgf(ix,iy,iz)=0.
                    w_buoy(ix,iy,iz)=0.
                    w_xdiff(ix,iy,iz)=0.
                    w_ydiff(ix,iy,iz)=0.
                    w_zdiff(ix,iy,iz)=0.
                    w_tend_total(ix,iy,iz)=0.
                    thp_xadv(ix,iy,iz)=0.
                    thp_yadv(ix,iy,iz)=0.
                    thp_zadv(ix,iy,iz)=0.
                    thp_meanadv(ix,iy,iz)=0.
                    thp_xdiff(ix,iy,iz)=0.
                    thp_ydiff(ix,iy,iz)=0.
                    thp_zdiff(ix,iy,iz)=0.
                    thp_tend_total(ix,iy,iz)=0.
                    pip_xadv(ix,iy,iz)=0.
                    pip_yadv(ix,iy,iz)=0.
                    pip_zadv(ix,iy,iz)=0.
                    pip_xdiff(ix,iy,iz)=0.
                    pip_ydiff(ix,iy,iz)=0.
                    pip_zdiff(ix,iy,iz)=0.
                    pip_tend_total(ix,iy,iz)=0.
                    rvp_xadv(ix,iy,iz)=0.
                    rvp_yadv(ix,iy,iz)=0.
                    rvp_zadv(ix,iy,iz)=0.
                    rvp_meanadv(ix,iy,iz)=0.
                    rvp_xdiff(ix,iy,iz)=0.
                    rvp_ydiff(ix,iy,iz)=0.
                    rvp_zdiff(ix,iy,iz)=0.
                    rvp_tend_total(ix,iy,iz)=0.
                    rcp_xadv(ix,iy,iz)=0.
                    rcp_yadv(ix,iy,iz)=0.
                    rcp_zadv(ix,iy,iz)=0.
                    rcp_xdiff(ix,iy,iz)=0.
                    rcp_ydiff(ix,iy,iz)=0.
                    rcp_zdiff(ix,iy,iz)=0.
                    rcp_tend_total(ix,iy,iz)=0.
                    rrp_xadv(ix,iy,iz)=0.
                    rrp_yadv(ix,iy,iz)=0.
                    rrp_zadv(ix,iy,iz)=0.
                    rrp_xdiff(ix,iy,iz)=0.
                    rrp_ydiff(ix,iy,iz)=0.
                    rrp_zdiff(ix,iy,iz)=0.
                    rrp_tend_total(ix,iy,iz)=0.
                enddo
            enddo
        enddo 
    end subroutine zero_tends

    subroutine calc_tend_u
        use run_constants, only: nz,nx,ny,kmx,kmy,kmz,rdx,rdy,rdz
        use constants, only: cp
        use model_vars, only:up,vp,wp,pip,rhoub,rhowb,thvb                                       &
                            ,u_xadv,u_yadv,u_zadv,u_pgf,u_xdiff,u_ydiff,u_zdiff,u_tend_total                        

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate tendency in u
        ! this is equation 6 in HW4

        ! loop over real/unique points
        
        do ix = 2, nx-1
            do iy = 2, ny-1
                do iz = 2, nz-1 
                    ! term1 in u-tendency equation: x advection term = -d(uu)/dx
                    u_xadv(ix,iy,iz) = - rdx * (((0.5*(up(ix+1,iy,iz,2)+up(ix,iy,iz,2)))**2) &
                                            - ((0.5*(up(ix,iy,iz,2)+up(ix-1,iy,iz,2)))**2))

                    ! term2 in u-tendency equation: y advection term = -d(uv)/dy
                    u_yadv(ix,iy,iz) = - rdy                                                &
                                        * 0.25 * (((vp(ix-1,iy+1,iz,2)    +   vp(ix,iy+1,iz,2))  * (up(ix,iy,iz,2)      +   up(ix,iy+1,iz,2))) &
                                        -         ((vp(ix-1,iy,iz,2)      +   vp(ix,iy,iz,2))    * (up(ix,iy-1,iz,2)    +   up(ix,iy,iz,2))))
                    
                    ! term3 in u-tendency equation: vertical advection term = -1/rho * d(rho*u*w)/dz
                    u_zadv(ix,iy,iz) = - (rdz/(rhoub(iz)))                                                 &
                                        * 0.25 * ((rhowb(iz+1) *    (wp(ix-1,iy,iz+1,2)    +   wp(ix,iy,iz+1,2))  * (up(ix,iy,iz,2)    +   up(ix,iy,iz+1,2))) &
                                        -         (rhowb(iz)   *    (wp(ix-1,iy,iz,2)      +   wp(ix,iy,iz,2))    * (up(ix,iy,iz-1,2)    +   up(ix,iy,iz,2))))

                    ! term4 in u-tendency equation: pressure gradient term = -cp * thetav_base * d(pi_pert)/ dx
                    u_pgf(ix,iy,iz) = -cp*thvb(iz)*rdx                                                   & 
                                        * (pip(ix,iy,iz,2)-pip(ix-1,iy,iz,2))

                    ! term5 in u-tendency equation: x diffusion
                    u_xdiff(ix,iy,iz) = kmx * rdx * rdx * (up(ix-1,iy,iz,1) - (2*up(ix,iy,iz,1)) + up(ix+1,iy,iz,1))

                    ! term6 in u-tendency equation: y diffusion
                    u_ydiff(ix,iy,iz) = kmy * rdy * rdy * (up(ix,iy-1,iz,1) - (2*up(ix,iy,iz,1)) + up(ix,iy+1,iz,1))

                    ! term7 in u-tendency equation: vertical diffusion
                    u_zdiff(ix,iy,iz) = kmz * rdz * rdz * (up(ix,iy,iz-1,1) - (2*up(ix,iy,iz,1)) + up(ix,iy,iz+1,1))
                    
                    u_tend_total(ix,iy,iz) = u_xadv(ix,iy,iz) + u_yadv(ix,iy,iz) + u_zadv(ix,iy,iz) &
                                            + u_pgf(ix,iy,iz) + u_xdiff(ix,iy,iz) + u_ydiff(ix,iy,iz) + u_zdiff(ix,iy,iz)
                enddo ! end x loop
            enddo
        enddo ! end z loop

    end subroutine calc_tend_u

    subroutine calc_tend_v
        use run_constants, only: nz,nx,ny,kmx,kmy,kmz,rdx,rdy,rdz
        use constants, only: cp
        use model_vars, only:up,vp,wp,pip,rhoub,rhowb,thvb                                       &
                            ,v_xadv,v_yadv,v_zadv,v_pgf,v_xdiff,v_ydiff,v_zdiff,v_tend_total                        

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        
        ! calculate tendency in v
        ! this is equation 6 in HW4 but for v instead of u

        ! loop over real/unique points
        
        do ix = 2, nx-1
            do iy = 2, ny-1
                do iz = 2, nz-1 
                    ! term in v-tendency equation: x advection term = -d(uv)/dx
                    v_xadv(ix,iy,iz) = -rdx * 0.25 * (((up(ix+1,iy-1,iz,2)+up(ix+1,iy,iz,2))*(vp(ix,iy,iz,2)+vp(ix+1,iy,iz,2)))    &
                                                    - ((up(ix,iy-1,iz,2)+up(ix,iy,iz,2))*(vp(ix-1,iy,iz,2)+vp(ix,iy,iz,2))))
                                            
                    ! term in v-tendency equation: y advection term = -d(vv)/dy
                    v_yadv(ix,iy,iz) = - rdy * (((0.5*(vp(ix,iy+1,iz,2)+vp(ix,iy,iz,2)))**2) &
                                            - ((0.5*(vp(ix,iy,iz,2)+vp(ix,iy-1 ,iz,2)))**2))

                    ! term in v-tendency equation: vertical advection term = -1/rho * d(rho*v*w)/dz
                    v_zadv(ix,iy,iz) = - (rdz/(rhoub(iz)))                                                 &
                                        * 0.25 * ((rhowb(iz+1) *    (wp(ix,iy-1,iz+1,2)    +   wp(ix,iy,iz+1,2))  * (vp(ix,iy,iz,2)    +   vp(ix,iy,iz+1,2))) &
                                        -         (rhowb(iz)   *    (wp(ix,iy-1,iz,2)      +   wp(ix,iy,iz,2))    * (vp(ix,iy,iz-1,2)    +   vp(ix,iy,iz,2))))

                    ! term in v-tendency equation: pressure gradient term = -cp * thetav_base * d(pi_pert)/ dy
                    v_pgf(ix,iy,iz) = -cp*thvb(iz)*rdy                                                  & 
                                        * (pip(ix,iy,iz,2)-pip(ix,iy-1,iz,2))

                    ! term in v-tendency equation: x diffusion
                    v_xdiff(ix,iy,iz) = kmx * rdx * rdx * (vp(ix-1,iy,iz,1) - (2*vp(ix,iy,iz,1)) + vp(ix+1,iy,iz,1))

                    ! term in v-tendency equation: y diffusion
                    v_ydiff(ix,iy,iz) = kmy * rdy * rdy * (vp(ix,iy-1,iz,1) - (2*vp(ix,iy,iz,1)) + vp(ix,iy+1,iz,1))

                    ! term in v-tendency equation: vertical diffusion
                    v_zdiff(ix,iy,iz) = kmz * rdz * rdz * (vp(ix,iy,iz-1,1) - (2*vp(ix,iy,iz,1)) + vp(ix,iy,iz+1,1))
                    
                    v_tend_total(ix,iy,iz) = v_xadv(ix,iy,iz) + v_yadv(ix,iy,iz) + v_zadv(ix,iy,iz) &
                                        + v_pgf(ix,iy,iz) +v_xdiff(ix,iy,iz) + v_ydiff(ix,iy,iz) + v_zdiff(ix,iy,iz)
                enddo ! end x loop
            enddo
        enddo ! end z loop

    end subroutine calc_tend_v

    subroutine calc_tend_w
        use run_constants, only: nz,nx,ny,kmx,kmy,kmz,rdx,rdy,rdz
        use constants, only: cp, g
        use model_vars, only:up,vp,wp,pip,thp,rhoub,rhowb,thvb,thb,rvp,rcp,rrp                               &
                            ,w_xadv,w_yadv,w_zadv,w_pgf,w_buoy,w_xdiff,w_ydiff,w_zdiff,w_tend_total                        

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate tendency in w
        ! this is equation 7 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iy = 2, ny-1
                do iz = 2, nz-1
                    ! term in w-tendency equation: x advection term 
                    w_xadv(ix,iy,iz) = - rdx                                                                              &
                                        * 0.25 * (((up(ix+1,iy,iz,2)+up(ix+1,iy,iz-1,2))*(wp(ix+1,iy,iz,2)+wp(ix,iy,iz,2)))         &                                        
                                                -((up(ix,iy,iz,2)+up(ix,iy,iz-1,2))*(wp(ix,iy,iz,2)+wp(ix-1,iy,iz,2))) )

                    ! term in w-tendency equation: y advection term 
                    w_yadv(ix,iy,iz) = - rdy                                                                              &
                                        * 0.25 * (((vp(ix,iy+1,iz,2)+vp(ix,iy+1,iz-1,2))*(wp(ix,iy+1,iz,2)+wp(ix,iy,iz,2)))         &                                        
                                                -((vp(ix,iy,iz,2)+vp(ix,iy,iz-1,2))*(wp(ix,iy,iz,2)+wp(ix,iy-1,iz,2))) )

                    ! term in w-tendency equation: vertical advection term 
                    w_zadv(ix,iy,iz) =  - rdz/rhoub(iz) &
                                            * ((rhowb(iz+1)*(0.5*(wp(ix,iy,iz+1,2)+wp(ix,iy,iz,2)))**2) &
                                            - (rhowb(iz)*(0.5*(wp(ix,iy,iz,2)+wp(ix,iy,iz-1,2)))**2))

                    ! term in w-tendency equation: pressure gradient term 
                    w_pgf(ix,iy,iz) = -cp * rdz * 0.25 * (thvb(iz)+thvb(iz-1)) * (pip(ix,iy,iz,2)-pip(ix,iy,iz-1,2))

                    ! term in w-tendency equation: buoyancy term 
                    w_buoy(ix,iy,iz) = g *((thp(ix,iy,iz,2)+thp(ix,iy,iz-1,2))/(thb(iz)+thb(iz-1))      &   
                                              + (0.61 * 0.5 * (rvp(ix,iy,iz,2)+rvp(ix,iy,iz-1,2)))                                &
                                              - (0.5 * (rcp(ix,iy,iz,2)+rcp(ix,iy,iz-1,2)) + (rrp(ix,iy,iz,2)+rrp(ix,iy,iz-1,2))))

                    ! term in w-tendency equation: x diffusion
                    w_xdiff(ix,iy,iz) = kmx * rdx * rdx * (wp(ix-1,iy,iz,1) - (2*wp(ix,iy,iz,1)) + wp(ix+1,iy,iz,1))

                    ! term in w-tendency equation:yx diffusion
                    w_ydiff(ix,iy,iz) = kmy * rdy * rdy * (wp(ix,iy-1,iz,1) - (2*wp(ix,iy,iz,1)) + wp(ix,iy+1,iz,1))

                    ! term in w-tendency equation: vertical diffusion
                    w_zdiff(ix,iy,iz) = kmz * rdz * rdz * (wp(ix,iy,iz-1,1) - (2*wp(ix,iy,iz,1)) + wp(ix,iy,iz+1,1))
                    
                    w_tend_total(ix,iy,iz) = w_xadv(ix,iy,iz) + w_yadv(ix,iy,iz) + w_zadv(ix,iy,iz) &
                                        + w_pgf(ix,iy,iz) + w_buoy(ix,iy,iz) &
                                        + w_xdiff(ix,iy,iz) +  w_ydiff(ix,iy,iz) + w_zdiff(ix,iy,iz)
                enddo
            enddo ! end x loop
        enddo ! end z loop
    end subroutine calc_tend_w

    subroutine calc_tend_thp
        use run_constants, only: nz,nx,ny,khx,khy,khz,rdx,rdy,rdz
        use model_vars, only:up,vp,wp,thp,rhoub,rhowb,thb                                       &
                            ,thp_xadv,thp_yadv,thp_zadv,thp_meanadv,thp_xdiff,thp_ydiff,thp_zdiff,thp_tend_total                        

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate tendency in theta
        ! this is equation 8 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iy = 2, ny-1
                do iz = 2, nz-1
                    ! term in thp-tendency equation: x advection term = -d(u * thp)/dx
                    thp_xadv(ix,iy,iz) = - rdx                                                                   &
                                        * 0.5 * ((up(ix+1,iy,iz,2) * (thp(ix+1,iy,iz,2) +  thp(ix,iy,iz,2)))     &
                                        -        (up(ix,iy,iz,2)   * (thp(ix,iy,iz,2)   +  thp(ix-1,iy,iz,2))))

                    ! term in thp-tendency equation: y advection term = -d(v * thp)/dy
                    thp_yadv(ix,iy,iz) = - rdy                                                                  &
                                        * 0.5 * ((vp(ix,iy+1,iz,2) * (thp(ix,iy+1,iz,2) +  thp(ix,iy,iz,2)))     &
                                        -        (vp(ix,iy,iz,2)   * (thp(ix,iy,iz,2)   +  thp(ix,iy-1,iz,2))))

                    ! term in thp-tendency equation: vertical advection term = -1/rho * d(rho*w * thp)/dz
                    thp_zadv(ix,iy,iz) = - rdz/(rhoub(iz))                                                 &
                                        * 0.5 *  ((rhowb(iz+1) *    wp(ix,iy,iz+1,2)     * (thp(ix,iy,iz+1,2)    +   thp(ix,iy,iz,2))) &
                                        -         (rhowb(iz)   *    wp(ix,iy,iz,2)     * (thp(ix,iy,iz,2)    +   thp(ix,iy,iz-1,2))))

                    ! term in thp-tendency equation: mean state advection
                    thp_meanadv(ix,iy,iz) = -0.5 * (1/rhoub(iz)) *   rdz                                     &
                                                * ((rhowb(iz) * wp(ix,iy,iz,2) * (thb(iz)-thb(iz-1)))      &
                                                +  (rhowb(iz+1) * wp(ix,iy,iz+1,2) * (thb(iz+1)-thb(iz))))

                    ! term in thp-tendency equation: x diffusion
                    thp_xdiff(ix,iy,iz) = khx * rdx * rdx * (thp(ix-1,iy,iz,1) - (2*thp(ix,iy,iz,1)) + thp(ix+1,iy,iz,1))

                    ! term in thp-tendency equation: y diffusion
                    thp_ydiff(ix,iy,iz) = khy * rdy * rdy * (thp(ix,iy-1,iz,1) - (2*thp(ix,iy,iz,1)) + thp(ix,iy+1,iz,1))

                    ! term in thp-tendency equation: z diffusion
                    thp_zdiff(ix,iy,iz) = khz * rdz * rdz * (thp(ix,iy,iz-1,1) - (2*thp(ix,iy,iz,1)) + thp(ix,iy,iz+1,1))

                    thp_tend_total(ix,iy,iz) = thp_xadv(ix,iy,iz) + thp_yadv(ix,iy,iz) + thp_zadv(ix,iy,iz) &
                                            + thp_meanadv(ix,iy,iz) + thp_xdiff(ix,iy,iz) + thp_ydiff(ix,iy,iz) + thp_zdiff(ix,iy,iz)
                enddo
            enddo ! end x loop
        enddo ! end z loop
    end subroutine calc_tend_thp

    subroutine calc_tend_pip
        use run_constants, only: nz,nx,ny,khx,khy,khz,cs,rdx,rdy,rdz
        use constants, only: cp
        use model_vars, only:up,vp,wp,pip,thp,rhoub,rhowb,thvb                                       &
                            ,pip_xadv,pip_yadv,pip_zadv,pip_xdiff,pip_ydiff,pip_zdiff,pip_tend_total                        

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate tendency in perturbation exner function
        ! this is equation 9 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iy = 2, ny-1
                do iz = 2, nz-1
                    ! term in pip-tendency equation: x advection term
                    pip_xadv(ix,iy,iz) = -((cs**2)*rdx/(cp*thvb(iz))) * (up(ix+1,iy,iz,2)-up(ix,iy,iz,2))

                    ! term in pip-tendency equation: yadvection term
                    pip_yadv(ix,iy,iz) = -((cs**2)*rdy/(cp*thvb(iz))) * (vp(ix,iy+1,iz,2)-vp(ix,iy,iz,2))
                    
                    ! term in pip-tendency equation: z advection term 
                    pip_zadv(ix,iy,iz) = -(((cs**2)*rdz*0.5/(rhoub(iz)*cp*(thvb(iz)**2))))   &
                                    *  ((rhowb(iz+1)*wp(ix,iy,iz+1,2)*(thvb(iz+1)+thvb(iz))) &
                                    -(rhowb(iz)*wp(ix,iy,iz,2)*(thvb(iz)+thvb(iz-1))))

                    ! term in pip-tendency equation: x diffusion
                    pip_xdiff(ix,iy,iz) = khx * rdx * rdx * (pip(ix-1,iy,iz,1) - (2*pip(ix,iy,iz,1)) + pip(ix+1,iy,iz,1))

                    ! term in pip-tendency equation: y diffusion
                    pip_ydiff(ix,iy,iz) = khy * rdy * rdy * (pip(ix,iy-1,iz,1) - (2*pip(ix,iy,iz,1)) + pip(ix,iy+1,iz,1))

                    ! term in pip-tendency equation: x diffusion
                    pip_zdiff(ix,iy,iz) = khz * rdz * rdz * (pip(ix,iy,iz-1,1) - (2*pip(ix,iy,iz,1)) + pip(ix,iy,iz+1,1))

                    pip_tend_total(ix,iy,iz) = pip_xadv(ix,iy,iz)+pip_yadv(ix,iy,iz)+pip_zadv(ix,iy,iz)&
                                            +pip_xdiff(ix,iy,iz)+pip_ydiff(ix,iy,iz)+pip_zdiff(ix,iy,iz)
                enddo ! end x loop
            enddo
        enddo ! end z loop
    end subroutine calc_tend_pip

    subroutine calc_tend_rvp
        use run_constants, only: nz,nx,ny,khx,khy,khz,rdx,rdy,rdz
        use model_vars, only:up,vp,wp,rvp,rvb,rhoub,rhowb,thb                                      &
                            ,rvp_xadv,rvp_yadv,rvp_zadv,rvp_meanadv,rvp_xdiff,rvp_ydiff,rvp_zdiff,rvp_tend_total                        

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate tendency in theta
        ! this is equation 8 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iy = 2, ny-1
                do iz = 2, nz-1
                    ! term in rvp-tendency equation: x advection term = -d(u * rvp)/dx
                    rvp_xadv(ix,iy,iz) = - rdx                                                                   &
                                        * 0.5 * ((up(ix+1,iy,iz,2) * (rvp(ix+1,iy,iz,2) +  rvp(ix,iy,iz,2)))     &
                                        -        (up(ix,iy,iz,2)   * (rvp(ix,iy,iz,2)   +  rvp(ix-1,iy,iz,2))))

                    ! term in rvp-tendency equation: y advection term = -d(v * rvp)/dy
                    rvp_yadv(ix,iy,iz) = - rdy                                                                  &
                                        * 0.5 * ((vp(ix,iy+1,iz,2) * (rvp(ix,iy+1,iz,2) +  rvp(ix,iy,iz,2)))     &
                                        -        (vp(ix,iy,iz,2)   * (rvp(ix,iy,iz,2)   +  rvp(ix,iy-1,iz,2))))

                    ! term in rvp-tendency equation: vertical advection term = -1/rho * d(rho*w * rvp)/dz
                    rvp_zadv(ix,iy,iz) = - rdz/(rhoub(iz))                                                 &
                                        * 0.5 *  ((rhowb(iz+1) *    wp(ix,iy,iz+1,2)     * (rvp(ix,iy,iz+1,2)    +   rvp(ix,iy,iz,2))) &
                                        -         (rhowb(iz)   *    wp(ix,iy,iz,2)     * (rvp(ix,iy,iz,2)    +   rvp(ix,iy,iz-1,2))))

                    ! term in rvp-tendency equation: mean state advection
                    rvp_meanadv(ix,iy,iz) = -0.5 * (1/rhoub(iz)) *   rdz                                     &
                                                * ((rhowb(iz) * wp(ix,iy,iz,2) * (rvb(iz)-rvb(iz-1)))      &
                                                +  (rhowb(iz+1) * wp(ix,iy,iz+1,2) * (rvb(iz+1)-rvb(iz))))

                    ! term in rvp-tendency equation: x diffusion
                    rvp_xdiff(ix,iy,iz) = khx * rdx * rdx * (rvp(ix-1,iy,iz,1) - (2*rvp(ix,iy,iz,1)) + rvp(ix+1,iy,iz,1))

                    ! term in rvp-tendency equation: y diffusion
                    rvp_ydiff(ix,iy,iz) = khy * rdy * rdy * (rvp(ix,iy-1,iz,1) - (2*rvp(ix,iy,iz,1)) + rvp(ix,iy+1,iz,1))

                    ! term in rvp-tendency equation: z diffusion
                    rvp_zdiff(ix,iy,iz) = khz * rdz * rdz * (rvp(ix,iy,iz-1,1) - (2*rvp(ix,iy,iz,1)) + rvp(ix,iy,iz+1,1))

                    rvp_tend_total(ix,iy,iz) = rvp_xadv(ix,iy,iz) + rvp_yadv(ix,iy,iz) + rvp_zadv(ix,iy,iz)                             &
                                            + rvp_meanadv(ix,iy,iz) + rvp_xdiff(ix,iy,iz) + rvp_ydiff(ix,iy,iz) + rvp_zdiff(ix,iy,iz)
                enddo
            enddo ! end x loop
        enddo ! end z loop
    end subroutine calc_tend_rvp

    subroutine calc_tend_rcp
        use run_constants, only: nz,nx,ny,khx,khy,khz,rdx,rdy,rdz
        use model_vars, only:up,vp,wp,rcp,rvb,rhoub,rhowb,thb                                      &
                            ,rcp_xadv,rcp_yadv,rcp_zadv,rcp_xdiff,rcp_ydiff,rcp_zdiff,rcp_tend_total                        

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate tendency in theta
        ! this is equation 8 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iy = 2, ny-1
                do iz = 2, nz-1
                    ! term in rcp-tendency equation: x advection term = -d(u * rcp)/dx
                    rcp_xadv(ix,iy,iz) = - rdx                                                                   &
                                        * 0.5 * ((up(ix+1,iy,iz,2) * (rcp(ix+1,iy,iz,2) +  rcp(ix,iy,iz,2)))     &
                                        -        (up(ix,iy,iz,2)   * (rcp(ix,iy,iz,2)   +  rcp(ix-1,iy,iz,2))))

                    ! term in rcp-tendency equation: y advection term = -d(v * rcp)/dy
                    rcp_yadv(ix,iy,iz) = - rdy                                                                  &
                                        * 0.5 * ((vp(ix,iy+1,iz,2) * (rcp(ix,iy+1,iz,2) +  rcp(ix,iy,iz,2)))     &
                                        -        (vp(ix,iy,iz,2)   * (rcp(ix,iy,iz,2)   +  rcp(ix,iy-1,iz,2))))

                    ! term in rcp-tendency equation: vertical advection term = -1/rho * d(rho*w * rcp)/dz
                    rcp_zadv(ix,iy,iz) = - rdz/(rhoub(iz))                                                 &
                                        * 0.5 *  ((rhowb(iz+1) *    wp(ix,iy,iz+1,2)     * (rcp(ix,iy,iz+1,2)    +   rcp(ix,iy,iz,2))) &
                                        -         (rhowb(iz)   *    wp(ix,iy,iz,2)     * (rcp(ix,iy,iz,2)    +   rcp(ix,iy,iz-1,2))))

                    ! term in rcp-tendency equation: x diffusion
                    rcp_xdiff(ix,iy,iz) = khx * rdx * rdx * (rcp(ix-1,iy,iz,1) - (2*rcp(ix,iy,iz,1)) + rcp(ix+1,iy,iz,1))

                    ! term in rcp-tendency equation: y diffusion
                    rcp_ydiff(ix,iy,iz) = khy * rdy * rdy * (rcp(ix,iy-1,iz,1) - (2*rcp(ix,iy,iz,1)) + rcp(ix,iy+1,iz,1))

                    ! term in rcp-tendency equation: z diffusion
                    rcp_zdiff(ix,iy,iz) = khz * rdz * rdz * (rcp(ix,iy,iz-1,1) - (2*rcp(ix,iy,iz,1)) + rcp(ix,iy,iz+1,1))

                    rcp_tend_total(ix,iy,iz) = rcp_xadv(ix,iy,iz) + rcp_yadv(ix,iy,iz) + rcp_zadv(ix,iy,iz) &
                                            +  rcp_xdiff(ix,iy,iz) + rcp_ydiff(ix,iy,iz) + rcp_zdiff(ix,iy,iz)
                enddo
            enddo ! end x loop
        enddo ! end z loop
    end subroutine calc_tend_rcp

    subroutine calc_tend_rrp
        use run_constants, only: nz,nx,ny,khx,khy,khz,rdx,rdy,rdz
        use model_vars, only:up,vp,wp,rrp,rvb,rhoub,rhowb,thb                                      &
                            ,rrp_xadv,rrp_yadv,rrp_zadv,rrp_xdiff,rrp_ydiff,rrp_zdiff,rrp_tend_total                        
        use microphysics_functions, only: rain_fallspeed

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate tendency in theta
        ! this is equation 8 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iy = 2, ny-1
                do iz = 2, nz-1
                    ! term in rrp-tendency equation: x advection term = -d(u * rrp)/dx
                    rrp_xadv(ix,iy,iz) = - rdx                                                                   &
                                        * 0.5 * ((up(ix+1,iy,iz,2) * (rrp(ix+1,iy,iz,2) +  rrp(ix,iy,iz,2)))     &
                                        -        (up(ix,iy,iz,2)   * (rrp(ix,iy,iz,2)   +  rrp(ix-1,iy,iz,2))))

                    ! term in rrp-tendency equation: y advection term = -d(v * rrp)/dy
                    rrp_yadv(ix,iy,iz) = - rdy                                                                  &
                                        * 0.5 * ((vp(ix,iy+1,iz,2) * (rrp(ix,iy+1,iz,2) +  rrp(ix,iy,iz,2)))     &
                                        -        (vp(ix,iy,iz,2)   * (rrp(ix,iy,iz,2)   +  rrp(ix,iy-1,iz,2))))

                    ! term in rrp-tendency equation: vertical advection term = -1/rho * d(rho*w * rrp)/dz
                    rrp_zadv(ix,iy,iz) = - rdz/(rhoub(iz))                                                 &
                                        * 0.5 *  ((rhowb(iz+1) *    wp(ix,iy,iz+1,2) &!-rain_fallspeed(rhowb(iz+1),0.5*(rrp(ix,iy,iz+1,2)+ rrp(ix,iy,iz,2))))     &
                                                                    * (rrp(ix,iy,iz+1,2)    +   rrp(ix,iy,iz,2)))                                              &
                                        -         (rhowb(iz)   *    wp(ix,iy,iz,2)&!-rain_fallspeed(rhowb(iz),0.5*(rrp(ix,iy,iz,2)+rrp(ix,iy,iz-1,2))))          &
                                                                      * (rrp(ix,iy,iz,2)    +   rrp(ix,iy,iz-1,2))))

                    ! term in rrp-tendency equation: x diffusion
                    rrp_xdiff(ix,iy,iz) = khx * rdx * rdx * (rrp(ix-1,iy,iz,1) - (2*rrp(ix,iy,iz,1)) + rrp(ix+1,iy,iz,1))

                    ! term in rrp-tendency equation: y diffusion
                    rrp_ydiff(ix,iy,iz) = khy * rdy * rdy * (rrp(ix,iy-1,iz,1) - (2*rrp(ix,iy,iz,1)) + rrp(ix,iy+1,iz,1))

                    ! term in rrp-tendency equation: z diffusion
                    rrp_zdiff(ix,iy,iz) = khz * rdz * rdz * (rrp(ix,iy,iz-1,1) - (2*rrp(ix,iy,iz,1)) + rrp(ix,iy,iz+1,1))

                    rrp_tend_total(ix,iy,iz) = rrp_xadv(ix,iy,iz) + rrp_yadv(ix,iy,iz) + rrp_zadv(ix,iy,iz) &
                                            +  rrp_xdiff(ix,iy,iz) + rrp_ydiff(ix,iy,iz) + rrp_zdiff(ix,iy,iz)
                enddo
            enddo ! end x loop
        enddo ! end z loop
    end subroutine calc_tend_rrp

    subroutine apply_tends
        use run_constants, only: nz,nx,ny,dt,d2t
        use model_vars, only:it,thp,pip,up,vp,wp,rvp,rcp,rrp                                   &
                            ,u_tend_total,v_tend_total,w_tend_total,thp_tend_total,pip_tend_total                    &
                            ,rvp_tend_total,rcp_tend_total,rrp_tend_total   

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        ! calculate actual future values
        ! loop over real/unique points
        do ix = 2, nx-1
            do iy = 2, ny-1
                do iz = 2, nz-1
                    if (it==1) then !for very frist timestep, just doing forward upstream
                        up(ix,iy,iz,3) = up(ix,iy,iz,2) + (dt * u_tend_total(ix,iy,iz))
                        vp(ix,iy,iz,3) = vp(ix,iy,iz,2) + (dt * v_tend_total(ix,iy,iz))
                        wp(ix,iy,iz,3) = wp(ix,iy,iz,2) + (dt * w_tend_total(ix,iy,iz))
                        thp(ix,iy,iz,3) = thp(ix,iy,iz,2) + (dt * thp_tend_total(ix,iy,iz))
                        pip(ix,iy,iz,3) = pip(ix,iy,iz,2) + (dt * pip_tend_total(ix,iy,iz))
                    else ! for other timesteps, leapfrog
                        up(ix,iy,iz,3) = up(ix,iy,iz,1) + (d2t * u_tend_total(ix,iy,iz))
                        vp(ix,iy,iz,3) = vp(ix,iy,iz,1) + (d2t * v_tend_total(ix,iy,iz))
                        wp(ix,iy,iz,3) = wp(ix,iy,iz,1) + (d2t * w_tend_total(ix,iy,iz))
                        thp(ix,iy,iz,3) = thp(ix,iy,iz,1) + (d2t * thp_tend_total(ix,iy,iz))
                        pip(ix,iy,iz,3) = pip(ix,iy,iz,1) + (d2t * pip_tend_total(ix,iy,iz))
                    endif
                enddo ! end z loop
            enddo ! end y loop
        enddo ! end x loop

    end subroutine apply_tends

end module solve_prog

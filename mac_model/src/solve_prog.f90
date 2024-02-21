module solve_prog
    !This module contains the full tendency equations for u,w,theta_prime,PI_prime from HW4

    implicit none

    contains
    
    subroutine tendencies
        use run_constants, only: nz,nx,dx,dz0,dt,pbc_x,pbc_z
        use constants, only: cp, cs, g
        use model_vars, only:it,thb,thvb,rhoub,rhowb,thp,pip,pp,up,wp &
                            ,u_tend1,u_tend2,u_tend3,u_tend_total &
                            ,w_tend1,w_tend2,w_tend3,w_tend4,w_tend_total &
                            ,thp_tend1,thp_tend2,thp_tend3,thp_tend_total &
                            ,pip_tend1,pip_tend2,pip_tend_total

        use boundaries, only: enforce_bounds_x, enforce_bounds_z

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate


        real :: rdx,rdz,d2t
        
        rdx      = (1/(dx))  ! reciprocal of dx [1/m]
        rdz      = (1/(dz0))  ! reciprocal of dz [1/m]
        d2t       = (dt+dt)         ! 2*dt

        ! first, reset all tendencies to zero
        do iz = 1, nz
            do ix = 1, nx
                u_tend1(ix,iz)=0.
                u_tend2(ix,iz)=0.
                u_tend3(ix,iz)=0.
                u_tend_total(ix,iz)=0.
                w_tend1(ix,iz)=0.
                w_tend2(ix,iz)=0.
                w_tend3(ix,iz)=0.
                w_tend4(ix,iz)=0.
                w_tend_total(ix,iz)=0.
                thp_tend1(ix,iz)=0.
                thp_tend2(ix,iz)=0.
                thp_tend3(ix,iz)=0.
                thp_tend_total(ix,iz)=0.
                pip_tend1(ix,iz)=0.
                pip_tend2(ix,iz)=0.
                pip_tend_total(ix,iz)=0.
            enddo
        enddo 

        ! calculate tendency in u
        ! this is equation 6 in HW4

        ! loop over real/unique points
        
        do ix = 2, nx-1
            do iz = 2, nz-1
                ! first term in u-tendency equation: horizontal advection term = -d(uu)/dx
                u_tend1(ix,iz) = - rdx * (((0.5*(up(ix+1,iz,2)+up(ix,iz,2)))**2) &
                                        - ((0.5*(up(ix,iz,2)+up(ix-1,iz,2)))**2))

                ! second term in u-tendency equation: vertical advection term = -1/rho * d(rho*u*w)/dz
                u_tend2(ix,iz) = - (rdz/(rhoub(iz)))                                                 &
                                    * 0.25 * ((rhowb(iz+1) *    (wp(ix-1,iz+1,2)    +   wp(ix,iz+1,2))  * (up(ix,iz,2)    +   up(ix,iz+1,2))) &
                                    -         (rhowb(iz)   *    (wp(ix-1,iz,2)      +   wp(ix,iz,2))    * (up(ix,iz-1,2)    +   up(ix,iz,2))))

                ! last term in u-tendency equation: pressure gradient term = -cp * theta_base * d(pi_pert)/ dx
                u_tend3(ix,iz) = -cp*thb(iz)*rdx                                                   & 
                                    * (pip(ix,iz,2)-pip(ix-1,iz,2))

                u_tend_total(ix,iz) = u_tend1(ix,iz) + u_tend2(ix,iz) + u_tend3(ix,iz)
            enddo ! end x loop
        enddo ! end z loop

        ! calculate tendency in w
        ! this is equation 7 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iz = 3, nz-1
                ! first term in w-tendency equation: horizontal advection term 
                w_tend1(ix,iz) = - rdx                                                                              &
                                    * 0.25 * (((up(ix+1,iz,2)+up(ix+1,iz-1,2))*(wp(ix+1,iz,2)+wp(ix,iz,2)))         &                                        
                                             -((up(ix,iz,2)+up(ix,iz-1,2))*(wp(ix,iz,2)+wp(ix-1,iz,2))) )

                ! second term in w-tendency equation: vertical advection term 
                w_tend2(ix,iz) =  - rdz/rhoub(iz) &
                                        * ((rhowb(iz+1)*(0.5*(wp(ix,iz+1,2)+wp(ix,iz,2)))**2) &
                                        - (rhowb(iz)*(0.5*(wp(ix,iz,2)+wp(ix,iz-1,2)))**2))

                ! third term in w-tendency equation: pressure gradient term 
                w_tend3(ix,iz) = -cp * rdz * 0.25 * (thb(iz)+thb(iz-1)) * (pip(ix,iz,2)-pip(ix,iz-1,2))

                ! last term in w-tendency equation: pressure gradient term 
                w_tend4(ix,iz) = g * (thp(ix,iz,2)+thp(ix,iz-1,2))/(thb(iz)+thb(iz-1))

                w_tend_total(ix,iz) = w_tend1(ix,iz) +w_tend2(ix,iz) +w_tend3(ix,iz) + w_tend4(ix,iz)
            enddo ! end x loop
        enddo ! end z loop

    
        ! calculate tendency in theta
        ! this is equation 8 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iz = 2, nz-1
                ! first term in thp-tendency equation: horizontal advection term = -d(u * thp)/dx
                thp_tend1(ix,iz) = - rdx                                                                   &
                                    * 0.5 * ((up(ix+1,iz,2) * (thp(ix+1,iz,2) +  thp(ix,iz,2)))     &
                                    -        (up(ix,iz,2)   * (thp(ix,iz,2)   +  thp(ix-1,iz,2))))

                ! second term in thp-tendency equation: vertical advection term = -1/rho * d(rho*w * thp)/dz
                thp_tend2(ix,iz) = - rdz/(rhoub(iz))                                                 &
                                    * 0.5 *  ((rhowb(iz+1) *    wp(ix,iz+1,2)     * (thp(ix,iz+1,2)    +   thp(ix,iz,2))) &
                                    -         (rhowb(iz)   *    wp(ix,iz,2)     * (thp(ix,iz,2)    +   thp(ix,iz-1,2))))

                ! last term in thp-tendency equation: pressure gradient term = -cp * theta_base * d(pi_pert)/ dx
                thp_tend3(ix,iz) = -0.5 * (1/rhoub(iz)) *   rdz                                     &
                                            * ((rhowb(iz) * wp(ix,iz,2) * (thb(iz)-thb(iz-1)))      &
                                            +  (rhowb(iz+1) * wp(ix,iz+1,2) * (thb(iz+1)-thb(iz))))

                thp_tend_total(ix,iz) = thp_tend1(ix,iz) + thp_tend2(ix,iz) + thp_tend3(ix,iz)
            enddo ! end x loop
        enddo ! end z loop


        ! calculate tendency in perturbation exner function
        ! this is equation 9 in HW4

        ! loop over real/unique points
        do ix = 2, nx-1
            do iz = 2, nz-1
                ! first term in pip-tendency equation: horizontal advection term
                pip_tend1(ix,iz) = -((cs**2)*rdx/(cp*thb(iz))) * (up(ix+1,iz,2)-up(ix,iz,2))
                !pip_tend1(ix,iz) = (rhoub(iz)*thb(iz)*rdx) * (up(ix+1,iz,2)-up(ix,iz,2))

                ! second term in pip-tendency equation: vertical advection term 
                pip_tend2(ix,iz) = -(((cs**2)*rdz*0.5/(rhoub(iz)*cp*(thb(iz)**2))))   &
                                *  ((rhowb(iz+1)*wp(ix,iz+1,2)*(thb(iz+1)+thb(iz))) &
                                  -(rhowb(iz)*wp(ix,iz,2)*(thb(iz)+thb(iz-1))))

                pip_tend_total(ix,iz) = pip_tend1(ix,iz)+pip_tend2(ix,iz)

                !pip_tend_total(ix,iz) = -((cs**2) /(rhoub(iz)*cp*((thb(iz))**2))) * (pip_tend1(ix,iz) + pip_tend2(ix,iz))

            enddo ! end x loop
        enddo ! end z loop

        ! calculate actual future values
        ! loop over real/unique points
        do ix = 2, nx-1
            do iz = 2, nz-1
                if (it==1) then !for very frist timestep
                    up(ix,iz,3) = up(ix,iz,2) + (dt * u_tend_total(ix,iz))
                    wp(ix,iz,3) = wp(ix,iz,2) + (dt * w_tend_total(ix,iz))
                    thp(ix,iz,3) = thp(ix,iz,2) + (dt * thp_tend_total(ix,iz))
                    pip(ix,iz,3) = pip(ix,iz,2) + (dt * pip_tend_total(ix,iz))
                else
                    up(ix,iz,3) = up(ix,iz,1) + (d2t * u_tend_total(ix,iz))
                    wp(ix,iz,3) = wp(ix,iz,1) + (d2t * w_tend_total(ix,iz))
                    thp(ix,iz,3) = thp(ix,iz,1) + (d2t * thp_tend_total(ix,iz))
                    pip(ix,iz,3) = pip(ix,iz,1) + (d2t * pip_tend_total(ix,iz))
                endif
            enddo ! end x loop
        enddo ! end z loop

    end subroutine tendencies

end module solve_prog
module base_state

    ! This module contains the subroutine to calculate the base state.

    implicit none

    contains

    subroutine init_base_state()

        use constants, only: g,cp,rd,p00,cv
        use run_constants, only: nz, nx, ztr, ttr, thtr, psurf, wk_flag, dn_flag
        use model_vars, only: zsn, zwn, tb, thb, rvb, thvb, pib, piwb, pb, rsatb, rhb, satfracb,rhoub,rhowb
        use thermo_functions, only: calc_thv, calc_rsat, calc_satfrac

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter fo x-coordinate

        ! initally, set the base state to be horizontally homogenous
        ! so we just need to loop over all the x values inside the z-loops        

        if (wk_flag==.True.) then
            ! set base state potential temperature to be equal to the WK sounding given
            do iz = 2,nz-1
                if (zsn(iz) <= ztr) then
                    do ix=2,nx-1
                        thb(ix,iz) = 300. + 43.*(zsn(iz)/ztr)**1.25
                    enddo
                else
                    do ix=2,nx-1
                        thb(ix,iz) = thtr * EXP((g*(zsn(iz)-ztr))/(ttr*cp))
                    enddo
                endif
            enddo
            
            ! set base state water vapor mixing ratio to given
            do iz = 2,nz-1
                if (zsn(iz) <= 4000.) then
                    do ix = 2,nx-1
                        rvb(ix,iz) = 1.61E-2 - 3.375E-6*zsn(iz)
                    enddo
                else if (zsn(iz) <= 8000.) then
                    do ix = 2,nx-1
                        rvb(ix,iz) = 2.6E-3 - 6.5E-7*(zsn(iz) - 4000.)
                    enddo
                else 
                    do ix = 2,nx-1
                        rvb(ix,iz) = 0
                    enddo
                endif
            enddo
            
        else if (dn_flag==.True.) then
            !set base state sounding to be dry and neutral
            do iz = 2,nz-1
                do ix = 2,nx-1
                    thb(ix,iz) = 300.
                enddo
            enddo
            
            ! set base state water vapor mixing ratio to 0 as given
            do iz = 2,nz-1
                do ix = 2,nx-1
                    rvb(ix,iz) = 0.
                enddo
            enddo
            
            
        endif 

        ! setting the boundary points same as first real points for zero gradient
        do iz=2,nz-1
            thb(1,iz) = thb(2,iz) 
            thb(nx,iz) = thb(nx-1,iz) 
            rvb(1,iz) = rvb(2,iz) 
            rvb(nx,iz) = rvb(nx-1,iz) 
        enddo
        do ix=1,nx
            thb(ix,1) = thb(ix,2) 
            thb(ix,nz) = thb(ix,nz-1) 
            rvb(ix,1) = rvb(ix,2) 
            rvb(ix,nz) = rvb(ix,nz-1) 
        enddo
        
        do ix = 2,nx-1
            ! base state virtual potential temp using the definition of theta_v
            do iz = 1,nz
                thvb(ix,iz) = calc_thv(thb(ix,iz), rvb(ix,iz)) 
            enddo

            ! base state exner function (non-dimensionalized pressure, pi)
            ! let's calculate both at u grids and w grids (which will be useful later)

            ! we know the surface pressure so can calculate PI for the surface (z=2 on w grid)
            ! this would just be (psurf/p00)**(rd/cp) (see definition of PI in HW1)
            piwb(ix,2) = (psurf/p00)**(rd/cp)

            ! but the grid for PI is offset by half a delta z (which is zun - zwn), 
            ! so we must integrate vertically to get PI(z=2) rather than PI(surf) 
            ! using the hydrostatic equation (see HW1) in terms of PI and thetav
            pib(ix,2) = piwb(ix,2)-((g/(cp*thvb(ix,2)))*(zsn(2)-zwn(2)))

            pib(ix,1) = pib(ix,2) ! the first level is fictitious, so just set it to be same as first real level (constant)
            piwb(ix,1) = piwb(ix,2) ! the first level is fictitious, so just set it to be same as first real level (constant)

            ! can do similar (i.e., integrate PI vertically) for the other levels
            ! for each vertical level, subtract dPI/dz * dz where dz = zun(iz) - zun(iz-1)
            ! and dPI/dz is inversely proportional to the average thetav between the two levels
            ! being integrated over (see hydrostatic equation as in HW1)
            ! note that we are basically assuming thvb scales linearly with height in between the scalar levels
            do iz = 3,nz-1
                pib(ix,iz) = pib(ix,iz-1) - (g/(cp*(thvb(ix,iz-1)+thvb(ix,iz))/2)) * (zsn(iz)-zsn(iz-1))
                piwb(ix,iz) = piwb(ix,iz-1) - (g/(cp*thvb(ix,iz-1))) * (zwn(iz)-zwn(iz-1))

                ! base state pressure
                pb(ix,iz) = p00 * pib(ix,iz)**(cp/rd)

                ! base state temperature
                tb(ix,iz) = pib(ix,iz)*thb(ix,iz)
            
                ! base state saturation vapor mixing ratio 
                ! using Teten's equation as written in HW1
                rsatb(ix,iz) = calc_rsat(tb(ix,iz),pb(ix,iz))

                ! base state saturation fraction
                ! remember RH is just rv/rsat
                satfracb(ix,iz) = calc_satfrac(rvb(ix,iz),rsatb(ix,iz))
                rhb(ix,iz) = satfracb(ix,iz)*100
            
                ! base state air density on u/scalar grid
                ! we can calculate this exactly from the ideal gas law as written in HW1
                rhoub(ix,iz) = (p00*pib(ix,iz)**(cv/rd))/(rd*thvb(ix,iz))

                ! base state air density on w grid

                ! though we have the form for density as above, remember that the w grid is offset
                ! from the u/scalar grid, so we need to use the piwb we calculated earlier

                ! for thv, assuming it scales linearlly between levels and use the average
                rhowb(ix,iz) = (p00*piwb(ix,iz)**(cv/rd))/(rd*(thvb(ix,iz)+thvb(ix,iz-1))/2)
            enddo

        enddo

        do iz = 1,nz
            !set the other fictitious points on model lateral boundaries for zero gradient
            thvb(1,iz) = thvb(2,iz)
            thvb(nx,iz) = thvb(nx-1,iz)
            pib(1,iz) = pib(2,iz)
            pib(nx,iz) = pib(nx-1,iz)
            piwb(1,iz) = piwb(2,iz)
            piwb(nx,iz) = piwb(nx-1,iz)
            pb(1,iz) = pb(2,iz)
            pb(nx,iz) = pb(nx-1,iz)
            tb(1,iz) = tb(2,iz)
            tb(nx,iz) = tb(nx-1,iz)
            rsatb(1,iz) = rsatb(2,iz)
            rsatb(nx,iz) = rsatb(nx-1,iz)
            satfracb(1,iz) = satfracb(2,iz)
            satfracb(nx,iz) = satfracb(nx-1,iz)
            rhb(1,iz) = rhb(2,iz)
            rhb(nx,iz) = rhb(nx-1,iz)
            rhoub(1,iz) = rhoub(2,iz)
            rhoub(nx,iz) = rhoub(nx-1,iz)
            rhowb(1,iz) = rhowb(2,iz)
            rhowb(nx,iz) = rhowb(nx-1,iz)
        enddo

        do ix = 1,nx
            !set the fictitious points at model top and bottom for zero gradient
            thvb(ix,nz) = thvb(ix,nz-1)
            thvb(ix,nz) = thvb(ix,nz-1)
            pib(ix,nz) = pib(ix,nz-1)
            piwb(ix,nz) = piwb(ix,nz-1)
            pb(ix,1) = pb(ix,2)
            pb(ix,nz) = pb(ix,nz-1)
            tb(ix,1) = tb(ix,2)
            tb(ix,nz) = tb(ix,nz-1)
            rsatb(ix,1) = rsatb(ix,2)
            rsatb(ix,nz) = rsatb(ix,nz-1)
            satfracb(ix,1) = satfracb(ix,2)
            satfracb(ix,nz) = satfracb(ix,nz-1)
            rhb(ix,1) = rhb(ix,2)
            rhb(ix,nz) = rhb(ix,nz-1)
            rhoub(ix,1) = rhoub(ix,2)
            rhoub(ix,nz) = rhoub(ix,nz-1)
            rhowb(ix,1) = rhowb(ix,2)
            rhowb(ix,nz) = rhowb(ix,nz-1)
        enddo

    end subroutine init_base_state

end module base_state

module base_state

    ! This module contains the subroutine to calculate the base state.

    implicit none

    contains

    subroutine init_base_state()

        use constants, only: g,cp,rd,p00,cv
        use run_constants, only: nz, nx, ztr, ttr, thtr, psurf, wk_flag, dn_flag, base_out
        use thermo_functions, only: calc_thv, calc_rsat, calc_satfrac
        use model_vars, only:zsn, zwn, thb, rvb, thvb, pib, piwb, rhoub,rhowb,tb,pb,rsatb,rhb,satfracb
    
        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate

        ! initally, set the base state to be horizontally homogenous
        ! so we just need to loop over all the x values inside the z-loops     
        if (wk_flag==.True.) then
            ! use WK sounding
            
            do iz = 2,nz-1
                ! set base state potential temperature to be equal to the WK sounding given
                if (zsn(iz) <= ztr) then
                    thb(iz) = 300. + 43.*(zsn(iz)/ztr)**1.25
                else
                    thb(iz) = thtr * EXP((g*(zsn(iz)-ztr))/(ttr*cp))
                endif
           
                ! set base state water vapor mixing ratio to given
                if (zsn(iz) <= 4000.) then
                    rvb(iz) = 1.61E-2 - 3.375E-6*zsn(iz)
                else if (zsn(iz) <= 8000.) then
                    rvb(iz) = 2.6E-3 - 6.5E-7*(zsn(iz) - 4000.)
                else 
                    rvb(iz) = 0.
                endif 

            enddo ! end z loop
        else if (dn_flag==.True.) then
            !set base state sounding to be dry and neutral
            do iz = 2,nz-1
                ! set base state th to be all 300K
                thb(iz) = 300.

                ! set base state water vapor mixing ratio to 0 as given
                rvb(iz) = 0.
            enddo !end z loop
        endif  ! end WK vs. dry neutral flag

        ! setting the boundary points same as first real points for zero gradient
        thb(1) = thb(2) 
        thb(nz) = thb(nz-1) 
        rvb(1) = rvb(2) 
        rvb(nz) = rvb(nz-1) 
        
        ! base state virtual potential temp using the definition of theta_v
        do iz = 1,nz
            thvb(iz) = calc_thv(thb(iz), rvb(iz)) 
        enddo ! end z loop

        ! base state exner function (non-dimensionalized pressure, pi)
        ! let's calculate both at u grids and w grids (which will be useful later)

        ! we know the surface pressure so can calculate PI for the surface (z=2 on w grid)
        ! this would just be (psurf/p00)**(rd/cp) (see definition of PI in HW1)
        piwb(2) = (psurf/p00)**(rd/cp)

        ! but the grid for PI is offset by half a delta z (which is zun - zwn), 
        ! so we must integrate vertically to get PI(z=2) rather than PI(surf) 
        ! using the hydrostatic equation (see HW1) in terms of PI and thetav
        pib(2) = piwb(2)-((g/(cp*thvb(2)))*(zsn(2)-zwn(2)))

        pib(1) = pib(2) ! the first level is fictitious, so just set it to be same as first real level (constant)
        piwb(1) = piwb(2) ! the first level is fictitious, so just set it to be same as first real level (constant)

        ! can do similar (i.e., integrate PI vertically) for the other levels
        ! for each vertical level, subtract dPI/dz * dz where dz = zun(iz) - zun(iz-1)
        ! and dPI/dz is inversely proportional to the average thetav between the two levels
        ! being integrated over (see hydrostatic equation as in HW1)
        ! note that we are basically assuming thvb scales linearly with height in between the scalar levels
        do iz = 3,nz-1
            pib(iz) = pib(iz-1) - (g/(cp*(thvb(iz-1)+thvb(iz))/2)) * (zsn(iz)-zsn(iz-1))
            piwb(iz) = piwb(iz-1) - (g/(cp*thvb(iz-1))) * (zwn(iz)-zwn(iz-1))
        enddo ! end z loop

        do iz=2,nz-1
            if (base_out) then
                ! base state pressure
                pb(iz) = p00 * pib(iz)**(cp/rd)

                ! base state temperature
                tb(iz) = pib(iz)*thb(iz)

                ! base state saturation vapor mixing ratio 
                ! using Teten's equation as written in HW1
                rsatb(iz) = calc_rsat(tb(iz),pb(iz))

                ! base state saturation fraction
                ! remember RH is just rv/rsat
                satfracb(iz) = calc_satfrac(rvb(iz),rsatb(iz))
                rhb(iz) = satfracb(iz)*100
            endif
            ! base state air density on u/scalar grid
            ! we can calculate this exactly from the ideal gas law as written in HW1
            rhoub(iz) = (p00*pib(iz)**(cv/rd))/(rd*thvb(iz))

            ! base state air density on w grid
            ! though we have the form for density as above, remember that the w grid is offset
            ! from the u/scalar grid, so we need to use the piwb we calculated earlier

            ! for thv, assuming it scales linearlly between levels and use the average
            !rhowb(iz) = (p00*piwb(iz)**(cv/rd))/(rd*(thvb(iz)+thvb(iz-1))/2)
            rhowb(iz) = (p00*((pib(iz)+pib(iz-1))/2)**(cv/rd))/(rd*(thvb(iz)+thvb(iz-1))/2)
        enddo !end z loop



        !set the fictitious points at model top and bottom for zero gradient
        thvb(nz) = thvb(nz-1)
        thvb(nz) = thvb(nz-1)
        pib(nz) = pib(nz-1)
        piwb(nz) = piwb(nz-1)
        rhoub(1) = rhoub(2)
        rhoub(nz) = rhoub(nz-1)
        rhowb(1) = rhowb(2)
        rhowb(nz) = rhowb(nz-1)

        if (base_out) then
            pb(1) = pb(2)
            pb(nz) = pb(nz-1)
            tb(1) = tb(2)
            tb(nz) = tb(nz-1)
            rsatb(1) = rsatb(2)
            rsatb(nz) = rsatb(nz-1)
            satfracb(1) = satfracb(2)
            satfracb(nz) = satfracb(nz-1)
            rhb(1) = rhb(2)
            rhb(nz) = rhb(nz-1)
        endif ! end base_out

    end subroutine init_base_state

end module base_state

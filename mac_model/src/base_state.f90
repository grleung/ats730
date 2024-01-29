module base_state

    ! This module contains the subroutine to calculate the base state.

    implicit none

    contains

    subroutine init_base_state()

        use constants, only: g,cp,rd,p00,cv
        use grid_constants, only: nz, ztr, ttr, thtr, psurf
        use model_vars, only: zun, zwn, tb, thb, rvb, thvb, pib, piwb, pb, rsatb, rhb, satfracb,rhoub,rhowb
        use thermo_functions, only: calc_thv, calc_rsat, calc_satfrac

        implicit none

        integer :: iz ! counter for z-coordinate

        ! set base state potential temperature to be equal to the WK sounding given
        do iz = 2,nz
            if (zun(iz) <= ztr) then
                thb(iz) = 300. + 43.*(zun(iz)/ztr)**1.25
            else
                thb(iz) = thtr * EXP((g*(zun(iz)-ztr))/(ttr*cp))
            endif
        enddo
        thb(1) = thb(2) ! the first level is fictitious, so just set it to be same as first real level (constant)

        ! set base state water vapor mixing ratio to given
        do iz = 2,nz
            if (zun(iz) <= 4000.) then
                rvb(iz) = 1.61E-2 - 3.375E-6*zun(iz)
            else if (zun(iz) <= 8000.) then
                rvb(iz) = 2.6E-3 - 6.5E-7*(zun(iz) - 4000.)
            else 
                rvb(iz) = 0
            endif
        enddo
        rvb(1) = rvb(2) ! the first level is fictitious, so just set it to be same as first real level (constant)

        ! base state virtual potential temp using the definition of theta_v
        do iz = 1,nz
            thvb(iz) = calc_thv(thb(iz), rvb(iz)) 
        enddo

        ! base state exner function (non-dimensionalized pressure, pi)
        ! let's calculate both at u grids and w grids (which will be useful later)

        ! we know the surface pressure so can calculate PI for the surface (z=2 on w grid)
        ! this would just be (psurf/p00)**(rd/cp) (see definition of PI in HW1)
        piwb(2) = (psurf/p00)**(rd/cp)

        ! but the grid for PI is offset by half a delta z (which is zun - zwn), 
        ! so we must integrate vertically to get PI(z=2) rather than PI(surf) 
        ! using the hydrostatic equation (see HW1) in terms of PI and thetav
        pib(2) = piwb(2)-((g/(cp*thvb(2)))*(zun(2)-zwn(2)))

        pib(1) = pib(2) ! the first level is fictitious, so just set it to be same as first real level (constant)
        piwb(1) = piwb(2) ! the first level is fictitious, so just set it to be same as first real level (constant)

        ! can do similar (i.e., integrate PI vertically) for the other levels
        ! for each vertical level, subtract dPI/dz * dz where dz = zun(iz) - zun(iz-1)
        ! and dPI/dz is inversely proportional to the average thetav between the two levels
        ! being integrated over (see hydrostatic equation as in HW1)
        ! note that we are basically assuming thvb scales linearly with height in between the scalar levels
        do iz = 3,nz
            pib(iz) = pib(iz-1) - (g/(cp*(thvb(iz-1)+thvb(iz))/2)) * (zun(iz)-zun(iz-1))
            piwb(iz) = piwb(iz-1) - (g/(cp*thvb(iz-1))) * (zwn(iz)-zwn(iz-1))
        enddo

        ! base state pressure
        do iz = 1,nz
            pb(iz) = p00 * pib(iz)**(cp/rd)
        enddo

        ! base state temperature
        do iz = 1,nz
            tb(iz) = pib(iz)*thb(iz)
        enddo
        
        ! base state saturation vapor mixing ratio 
        ! using Teten's equation as written in HW1
        do iz = 1,nz
            rsatb(iz) = calc_rsat(tb(iz),pb(iz))
        enddo

        ! base state saturation fraction
        ! remember RH is just rv/rsat
        do iz = 1,nz
            satfracb(iz) = calc_satfrac(rvb(iz),rsatb(iz))
            rhb(iz) = satfracb(iz)*100
        enddo

        ! base state air density on u/scalar grid
        ! we can calculate this exactly from the ideal gas law as written in HW1
        do iz = 1,nz
            rhoub(iz) = (p00*pib(iz)**(cv/rd))/(rd*thvb(iz))
        enddo

        ! base state air density on w grid

        ! though we have the form for density as above, remember that the w grid is offset
        ! from the u/scalar grid, so we need to use the piwb we calculated earlier

        ! for thv, we can assume it scales linearlly between levels and use the average
        ! between the value at iz and the next level like earlier integrating the exner function
        ! but pi doesn't scale linearly, so this isn't a good assumption
        ! so we also 
        do iz = 2,nz
            rhowb(iz) = (p00*piwb(iz)**(cv/rd))/(rd*(thvb(iz)+thvb(iz-1))/2)
        enddo

        rhowb(1) = rhowb(2)


    end subroutine init_base_state

end module base_state

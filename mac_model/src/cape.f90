module cape

    ! This module contains the subroutine to calculate CAPE

    implicit none

    contains

    subroutine calculate_parcel_cape()

        use constants, only: cp,lv
        use grid_constants, only: nz, ztr, ttr, thtr, psurf, rvpsurf
        use model_vars, only: zun, thb, rvb, thvb, pib, pb          &  ! base state variables
                            , thp, rvp, thvp, tp, rsatp             &  ! parcel state variables
                            , lclp, elp, capep                         ! parcel sounding parameters
        use thermo_functions, only: calc_rsat, calc_satfrac, calc_thv

        implicit none

        integer :: iz ! counter for z-coordinate
        real    :: phi, c ! saturation adjustment variables

        ! parcel starts at the lowest real scalar point (iz=2 on u-grid)
        ! set the initial parcel vapor mixing ratio (rvp) as given rvpsurf 
        ! and theta (thp) as equal to base state theta at iz=2
        thp(2) = thb(2)
        rvp(2) = rvpsurf

        ! set the fictitious points to just be same as first real level
        thp(1) = thp(2)
        rvp(1) = rvp(2)

        do iz = 2, nz ! loop through vertical levels as parcel is being lifted
            ! for unsaturated parcels, just lift dry adiabatically (conserve theta and rv)
            ! for saturated parcels, lift dry adiabatically then do saturation adjustment 

            ! calculate new temp and saturation pressure at this vertical level
            ! note that by definition of parcel, we assume parcel instantaneously adjusts 
            ! to the environmental pressure, so we can just use base state PI and pressure
            ! to calculate the parcel temp/rsat
            tp(iz) = pib(iz) * thp(iz)
            rsatp(iz) = calc_rsat(tp(iz), pb(iz))

            ! if parcel is saturated, then need to do saturation adjustment
            if (rvp(iz) >= rsatp(iz)) then
                ! if this is first time parcel is saturated,
                ! record the LCL as being the last model level (actually is somewhere
                ! between last model level and this level, because we aren't
                ! raising it continuously)
                if (lclp==0.0) lclp = iz-1

                ! account for change in parcel temp due to latent heat from condensation
                phi = rsatp(iz) * (17.27*237.*lv)/(cp*(tp(iz)-36.)**2) 

                ! actual amount condensed accounting for latent heat
                c = (rvp(iz) - rsatp(iz))/(1+phi) 

                ! raise temp by latent heat from condensing c amount of water
                thp(iz) = thp(iz) + (lv/(cp*pib(iz)))*c
                ! reduce water vapor mixing ratio by amount water condensed
                rvp(iz) = rvp(iz) - c 
                ! make sure temp is consistent with new theta and rv
                tp(iz) = pib(iz) * thp(iz)

            end if

            ! calculate virtual potential temp (thv)
            thvp(iz) = calc_thv(thp(iz),rvp(iz))

            print'(1I2,6F10.3)', iz, zun(iz), thp(iz), thvp(iz), thvp(iz)-thvb(iz), rvp(iz)*1000, rsatp(iz)*1000

            ! lift parcel dry adiabatically to next level
            thp(iz+1) = thp(iz)
            rvp(iz+1) = rvp(iz)
        enddo 

        print*, 'LCL is between', zun(lclp), zun(lclp+1)
        
        ! set the fictitious points to just be same as first real level
        tp(1) = tp(2)
        rsatp(1) = rsatp(2)
        
        

    end subroutine calculate_parcel_cape

end module cape

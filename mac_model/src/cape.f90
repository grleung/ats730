module cape

    ! This module contains the subroutine to calculate CAPE

    implicit none

    contains

    subroutine calculate_parcel_cape()

        use constants, only: g,cp,lv
        use run_constants, only: nz, ztr, ttr, thtr, psurf, rvpsurf
        use model_vars, only: zsn, thb, rvb, thvb, pib, pb          &  ! base state variables
                            , thpcl, rvpcl, thvpcl, tpcl, rsatpcl             &  ! parcel state variables
                            , thvdiff, lclp, elp, capep               ! parcel sounding parameters
        use thermo_functions, only: calc_rsat, calc_satfrac, calc_thv

        implicit none

        integer :: iz ! counter for z-coordinate
        real    :: phi, c ! saturation adjustment variables

        ! parcel starts at the lowest real scalar point (iz=2 on u-grid)
        ! set the initial parcel vapor mixing ratio (rvp) as given rvpsurf 
        ! and theta (thp) as equal to base state theta at iz=2
        thpcl(2) = thb(2)
        rvpcl(2) = rvpsurf

        do iz = 2, nz ! loop through vertical levels as parcel is being lifted
            ! for unsaturated parcels, just lift dry adiabatically (conserve theta and rv)
            ! for saturated parcels, lift dry adiabatically then do saturation adjustment 

            ! calculate new temp and saturation pressure at this vertical level
            ! note that by definition of parcel, we assume parcel instantaneously adjusts 
            ! to the environmental pressure, so we can just use base state PI and pressure
            ! to calculate the parcel temp/rsat
            tpcl(iz) = pib(iz) * thpcl(iz)
            rsatpcl(iz) = calc_rsat(tpcl(iz), pb(iz))

            ! if parcel is saturated, then need to do saturation adjustment
            if (rvpcl(iz) >= rsatpcl(iz)) then
                ! account for change in parcel temp due to latent heat from condensation
                phi = rsatpcl(iz) * (17.27*237.*lv)/(cp*(tpcl(iz)-36.)**2) 

                ! actual amount condensed accounting for latent heat
                c = (rvpcl(iz) - rsatpcl(iz))/(1+phi) 

                ! raise temp by latent heat from condensing c amount of water
                thpcl(iz) = thpcl(iz) + (lv/(cp*pib(iz)))*c
                ! reduce water vapor mixing ratio by amount water condensed
                rvpcl(iz) = rvpcl(iz) - c 
                ! make sure temp is consistent with new theta and rv
                tpcl(iz) = pib(iz) * thpcl(iz)

            end if

            ! calculate virtual potential temp (thv)
            thvpcl(iz) = calc_thv(thpcl(iz),rvpcl(iz))
            thvdiff(iz) = thvpcl(iz) - thvb(iz)

            ! to calculate CAPE, check theta_v difference between parcel and base state
            ! when parcel is positively buoyant (thvdiff>0) then it contributes to CAPE
            if (thvdiff(iz) > 0 ) then
                ! if this is the first time the parcel is positively buoyant, save as LCL
                if (lclp == 0) then 
                    lclp = iz-1
                endif

                !add to CAPE
                capep = capep + g*(thvdiff(iz)/thvb(iz)) * (zsn(iz)-zsn(iz-1))
            end if
            
            ! the first time that parcel is negatively buoyant after LCL is the EL
            if ((lclp>0) .and. (elp==0) .and. thvdiff(iz)<=0) then
                elp = iz-1
            endif

            ! lift parcel dry adiabatically to next level
            thpcl(iz+1) = thpcl(iz)
            rvpcl(iz+1) = rvpcl(iz)
        enddo 

        ! set the fictitious points to just be same as first real level
        thpcl(1) = thpcl(2)
        rvpcl(1) = rvpcl(2)
        tpcl(1) = tpcl(2)
        rsatpcl(1) = rsatpcl(2)
        thvpcl(1) = thvpcl(2)
        thvdiff(1) = thvdiff(2)

    end subroutine calculate_parcel_cape

end module cape

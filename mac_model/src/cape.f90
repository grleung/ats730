module cape

    ! This module contains the subroutine to calculate CAPE

    implicit none

    contains

    subroutine calculate_parcel_cape()

        use constants, only: g,cp,rd,p00,cv
        use grid_constants, only: nz, ztr, ttr, thtr, psurf, rvpsurf
        use model_vars, only: zun, thb, rvb, thvb, pib, thp, rvp, thvp, capep

        implicit none

        integer :: iz ! counter for z-coordinate

        ! parcel starts at the lowest real scalar point (iz=2 on u-grid)
        ! set the initial parcel vapor mixing ratio (rvp) as given rvpsurf 
        ! and theta (thp) as equal to base state theta at iz=2
        thp(2) = thb(2)
        rvp(2) = rvpsurf

        ! set the fictitious points to just be same as first real level
        thp(1) = thp(2)
        rvp(1) = rvp(2)

        ! check that the parcel is unsaturated at iz=2
        

    end subroutine calculate_parcel_cape

end module cape

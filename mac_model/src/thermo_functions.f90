module thermo_functions

    implicit none

    contains


    real function calc_thv(th,rv)
        implicit none
        real :: th, rv

        ! calculate virtual potential temp (thv)
        ! takes potential temperature (th) [K] and vapor mixing ratio (rv) [kg/kg]
        ! to output thv [K]
        calc_thv = th * (1+(0.61*rv))
    end function calc_thv

    real function calc_rsat(t,p)
        implicit none
        real :: t, p

        ! calculate saturation vapor mixing ratio 
        ! using Teten's equation as written in HW1

        ! takes temperature (t) [K] and pressure (p) [Pa] 
        ! to calculate saturation vapor mixing ratio (rsat) [kg/kg]
        calc_rsat = (380./p) * EXP((17.27*(t-273.))/(t-36.))
    end function calc_rsat

    real function calc_satfrac(rv, rsat)
        implicit none
        real :: rv, rsat

        ! calculate saturation fraction

        ! takes vapor mixing ratio (r) [kg/kg] and saturation vapor mixing ratio (rsat) [kg/kg] 
        ! to return satfrac [fraction]
        calc_satfrac = rv/rsat
    end function calc_satfrac

end module thermo_functions
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

    real function rain_fallspeed(rho, rr)
        use run_constants, only: minrain
        use constants, only: g, rhol,trigpi

        implicit none
        real :: rho,rr,lambda, k=1.83, N=8.e6 ! N0 is the slope intercept [/m4]

        ! calculates terminal velocity of raindrops 
        ! given air density (rhob) and slope parameter (lambda)

        ! check that rain is above some minimum value, or else there will be inf values of lambda
        if (rr >= minrain) then
            lambda = ((rhol*N*trigpi)/(rho*rr))**.25
            rain_fallspeed = k * (g*rhol/rho)**.5 * GAMMA(4.5) * lambda**(-0.5)
        else
            rain_fallspeed = 0.
        endif 
    end function rain_fallspeed

    real function rain_vent(rho, rr)
        implicit none
        real :: rho,rr

        rain_vent = 1.6 + 30.39*(rho*rr)**.2046
    end function rain_vent

end module thermo_functions
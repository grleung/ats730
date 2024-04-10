module micro_functions
    !This module contains functions for doing the microphysical calculations

    implicit none

    contains

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
    

end module micro_functions
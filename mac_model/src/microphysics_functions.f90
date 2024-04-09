module microphysics_functions
    !This module contains subroutines for doing the microphysical calculations

    implicit none

    contains

    real function rain_fallspeed(rho, rr)
        use constants, only: g, rhol,trigpi

        implicit none
        real :: rho,rr,lambda, k=1.83, N=8.e6 ! N0 is the slope intercept [/m4]

        ! calculates terminal velocity of raindrops 
        ! given air density (rhob) and slope parameter (lambda)

        lambda = ((rhol*N*trigpi)/(rho*rr))**.25

        rain_fallspeed = k * (g*rhol/rho)**.5 * GAMMA(4.5) * lambda**(-0.5)

    end function rain_fallspeed

end module microphysics_functions
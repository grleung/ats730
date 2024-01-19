module base_state

    ! This module driver contains the driver_subroutine

    implicit none

    contains

    subroutine init_base_state()

        ! Use the 'constants' module in order to get the physical constants
        use constants, only: g,cp
        use grid_constants, only: nz, ztr, ttr, thtr, psurf
        use model_vars, only: zun, zwn, tb, thb

        implicit none

        integer :: iz ! counter for z-coordinate

        ! print statement to check what's happening
        print*,'Gravity'
        print*,g

        ! print statement to check we are getting right model heights
        print*,'zwn'
        print*,zwn

        ! set the potential temperature to be equal to the WK sounding given
        thb(1) = 0. ! this level is fictitious
        do iz = 2,nz
            if (zun(iz) <= ztr) then
                thb(iz) = 300. + 43.*(zun(iz)/ztr)**1.25
            else
                thb(iz) = thtr * EXP((g*(zun(iz)-ztr))/(ttr*cp))
            endif
        enddo

        ! print statement to check we are getting right values
        print*,'thb'
        print*,thb

    end subroutine init_base_state

end module base_state
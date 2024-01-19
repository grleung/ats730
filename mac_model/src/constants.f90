module constants

    ! This module has all of the physical constants. 
    
    implicit none

    real, parameter :: &
            g       =    9.80       & ! gravitational acceleration [m/s2]
        ,   rd      =    287.       & ! dry air gas constant [J/kgK]
        ,   cp      =    1004.      & ! dry air specific heat capacity at constant pressure [J/kgK]
        ,   cv      =    cp - rd    & ! dry air specific heat capacity at constant volume [J/kgK]
        ,   p00     =    1.e5         ! reference pressure [Pa]

end module constants
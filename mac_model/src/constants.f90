module constants

    ! This module has all of the physical constants. 
    
    implicit none

    real, parameter :: &
            g       =    9.80       & ! gravitational acceleration [m/s2]
        ,   R       = 8.3145        & ! ideal gas constant
        ,   rd      =    287.       & ! dry air gas constant [J/kgK]
        ,   cp      =    1004.      & ! dry air specific heat capacity at constant pressure [J/kgK]
        ,   cv      =    cp - rd    & ! dry air specific heat capacity at constant volume [J/kgK]
        ,   p00     =    1.e5       & ! reference pressure [Pa]
        ,   lv      =    2.5e6      & ! latent heat of vaporization [J/kg] 
        ,   trigpi  =    4.*ATAN(1.)&! pi = 3.1415...  
        ,   rhol    =   1.e3        &
        ,   MW_w    =   0.018       &   ! water molecular weight [kg/mol]
        ,   MW_air    =   0.02897       &   ! air molecular weight [kg/mol]
        ,   mu      =   1.8e-5      &   ! air viscosity 
        ,   dg      =     1E-5      & !diffusivity m2/s
        ,   kair    =    2.40E-2    &   !thermal conductivity 
        ,   sigma_w =   0.073          ! water surf tension [J/m2]

end module constants
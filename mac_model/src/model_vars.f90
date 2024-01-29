module model_vars

    use grid_constants, only: nz

    ! This module is where we will declare all the shared arrays 
    
    implicit none

    ! grid coordinate variables
    real, dimension(nz)  :: &
            dzn    & !  deltaZ (level depth) of the w ("momentum") grid [m]
        ,   zun    & !  physical height of u/scalar ("u" for u) grid [m]
        ,   zwn      !  physical height of vertical velocity ("w" for w) grid [m]

    ! base state thermodynamic variables
    real, dimension(nz)  :: &
            tb          & !  base state temperature ("t" for temperature) [K]
        ,   thb         & !  base state potential temperature ("th" for theta) [K]
        ,   rvb         & !  base state water vapor mixing ratio ("r" for ratio, "v" for vapor) [kg/kg]
        ,   thvb        & !  base state virtual potential temperature [K]
        ,   pb          & !  base state pressure [Pa]
        ,   pib         & !  base state non-dimensional pressure on u-grid [no units]
        ,   piwb        & !  base state non-dimensional pressure on w-grid [no units]
        ,   rhoub       & !  base state air density at u/scalar levels [kg/m3]
        ,   rhowb       & !  base state air density at w levels[kg/m3]
        ,   satfracb    & !  base state saturation fraction [frac]
        ,   rhb         & !  base state relative humidity [%]
        ,   rsatb         !  base state saturation mixing ratio [kg/kg]

    ! parcel thermodynamic variables
    real, dimension(nz)  :: &
            tp          & ! parcel temperature ("t" for temperature) [K]
        ,   thp         & ! parcel potential temperature ("th" for theta) [K]
        ,   rvp         & ! parcel water vapor mixing ratio ("r" for ratio, "v" for vapor) [kg/kg]
        ,   thvp        & ! parcel virtual potential temperature [K]
        ,   pp          & !  parcel pressure [Pa]
        ,   pip         & !  parcel non-dimensional pressure on u-grid [no units]
        ,   piwp        & !  parcel non-dimensional pressure on w-grid [no units]
        ,   rhoup       & !  parcel air density at u/scalar levels [kg/m3]
        ,   satfracp    & !  parcel saturation fraction [frac]
        ,   rhp         & !  parcel relative humidity [%]
        ,   rsatp         !  parcel saturation mixing ratio [kg/kg]
    
    real                 :: capep     !  parcel CAPE [J/kg]
    

end module model_vars
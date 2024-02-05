module model_vars

    ! This module is where we will declare all the shared arrays 
    
    implicit none

    ! grid coordinate variables
    real, allocatable, dimension(:)  :: &
            dzn    & !  deltaZ (level depth) of the w  grid [m]
        ,   zsn    & !  physical height of u/scalar ("s" for scalar) grid [m]
        ,   zwn    & !  physical height of vertical velocity ("w" for w) grid [m]
        ,   dxn    & !  deltaX of the u grid [m]
        ,   xsn    & !  physical horizontal position of scalar ("s" for scalar) grid [m]
        ,   xun      !  physical horizontal position of horizontal velocity ("u" for u) grid [m]

    ! base state thermodynamic variables
    real, allocatable, dimension(:,:)  :: &
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
    real, allocatable, dimension(:)  :: &
            tp          & !  parcel temperature ("t" for temperature) [K]
        ,   thp         & !  parcel potential temperature ("th" for theta) [K]
        ,   rvp         & !  parcel water vapor mixing ratio ("r" for ratio, "v" for vapor) [kg/kg]
        ,   thvp        & !  parcel virtual potential temperature [K]
        ,   thvdiff     & !  parcel perturbation virtual potential temperature [K]
        ,   rsatp         !  parcel saturation mixing ratio [kg/kg]
    
    real                 :: &
            capep        !  parcel CAPE [J/kg]

    integer             :: &
           lclp         & !  parcel lifted condensation level [in model levels, u-grid]
        ,   elp           !  parcel equlibirum level [in model levels, u-grid]

    contains

end module model_vars
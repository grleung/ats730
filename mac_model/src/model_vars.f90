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

    ! base state thermodynamic variables (only profile in nz)
    real, allocatable, dimension(:)  :: &
            thb         & !  base state potential temperature ("th" for theta) [K]
        ,   thvb         & !  base state virtual potential temperature ("th" for theta) [K]
        ,   rvb         & !  base state water vapor mixing ratio ("r" for ratio, "v" for vapor) [kg/kg]
        ,   pib         & !  base state non-dimensional pressure on u-grid [no units]
        ,   piwb        & !  base state non-dimensional pressure on w-grid [no units]
        ,   rhoub       & !  base state air density at u/scalar levels [kg/m3]
        ,   rhowb       &  !  base state air density at w levels[kg/m3]
        ,   tb          & ! base state temp [K]
        ,   pb          & ! base state pressure [Pa]
        ,   rsatb       & ! base state saturation vapor pressure [Pa]
        ,   rhb         & ! base state relative humidity [%]
        ,   satfracb      ! base state saturation fraction

    ! parcel thermodynamic variables (only profile in nz)
    real, allocatable, dimension(:)  :: &
            tpcl          & !  parcel temperature ("t" for temperature) [K]
        ,   thpcl         & !  parcel potential temperature ("th" for theta) [K]
        ,   rvpcl         & !  parcel water vapor mixing ratio ("r" for ratio, "v" for vapor) [kg/kg]
        ,   thvpcl        & !  parcel virtual potential temperature [K]
        ,   thvdiff     & !  parcel perturbation virtual potential temperature [K]
        ,   rsatpcl         !  parcel saturation mixing ratio [kg/kg]

    ! parcel sounding derived parameters
    real                 :: &
            capep        !  parcel CAPE [J/kg]
    integer             :: &
           lclp         & !  parcel lifted condensation level [in model levels, u-grid]
        ,   elp           !  parcel equlibirum level [in model levels, u-grid]

    ! prognostic thermodynamic variables (array in nx,nz)
    real, allocatable, dimension(:,:)  :: &
            thp         & !  perturbation potential temperature ("th" for theta) [K]
        ,   rvp         & !  perturbation water vapor mixing ratio ("r" for ratio, "v" for vapor) [kg/kg]
        ,   pip         & !  perturbation non-dimensional pressure on u-grid [no units]
        ,   pp         & !  perturbation  pressure on u-grid [no units]
        ,   up          & !  horizontal velocity [m/s]
        ,   wp            !  vertical velocity [m/s]

    contains

end module model_vars
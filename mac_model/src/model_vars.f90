module model_vars

    ! This module is where we will declare all the shared arrays 
    
    implicit none

    integer :: it !timestep counter

    ! grid coordinate variables
    real, allocatable, dimension(:)  :: &
            dzn    & !  deltaZ (level depth) of the w  grid [m]
        ,   zsn    & !  physical height of u/scalar ("s" for scalar) grid [m]
        ,   zmn    & !  physical height of vertical velocity ("m" for momentum) grid [m]
        ,   dxn    & !  deltaX of the u grid [m]
        ,   xsn    & !  physical x position of scalar ("s" for scalar) grid [m]
        ,   xmn    & !  physical x position of horizontal velocity ("m" for momentum) grid [m]
        ,   dyn    & !  deltaX of the u grid [m]
        ,   ysn    & !  physical y position of scalar ("s" for scalar) grid [m]
        ,   ymn      !  physical y position of horizontal velocity ("m" for momentum) grid [m]

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

    ! prognostic thermodynamic variables (array in nx,ny,nz,3 time dims [past, pres, future])
    real, allocatable, dimension(:,:,:,:)  :: &
            thp         & !  perturbation potential temperature ("th" for theta) [K]
        ,   rvp         & !  perturbation water vapor mixing ratio ("r" for ratio, "v" for vapor) [kg/kg]
        ,   pip         & !  perturbation non-dimensional pressure on u-grid [no units]
        ,   up          & !  x (E-W) velocity [m/s]
        ,   vp          & !  y (N-S) velocity [m/s]
        ,   wp            !  vertical velocity [m/s]

    !(array in nx,ny,nz)    
    real, allocatable, dimension(:,:,:) :: &
           pp         !  perturbation  pressure on u-grid [no units]

    ! tendency  variables (array in nx,ny,nz)
    real, allocatable, dimension(:,:,:)  :: &
             u_xadv,u_yadv,u_zadv,u_pgf,u_xdiff,u_ydiff,u_zdiff,u_tend_total                        &
            ,v_xadv,v_yadv,v_zadv,v_pgf,v_xdiff,v_ydiff,v_zdiff,v_tend_total                        &
            ,w_xadv,w_yadv,w_zadv,w_pgf,w_buoy,w_xdiff,w_ydiff,w_zdiff,w_tend_total                 &
            ,thp_xadv,thp_yadv,thp_zadv,thp_meanadv,thp_xdiff,thp_ydiff,thp_zdiff,thp_tend_total    &
            ,pip_xadv,pip_yadv,pip_zadv,pip_xdiff,pip_ydiff,pip_zdiff,pip_tend_total

    contains

end module model_vars
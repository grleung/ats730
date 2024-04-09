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
        ,   rhowb       & !  base state air density at w levels[kg/m3]
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
        ,   wp          & !  vertical velocity [m/s]
        ,   rcp         & !  cloud water mixing ratio ('r for ratio, 'c' for cloud), tehcnically a perturbation but from base state 0 [kg/kg]
        ,   rrp           !  rain water mixing ratio ('r for ratio, 'r' for rain), tehcnically a perturbation but from base state 0 [kg/kg]

    !(array in nx,ny,nz)    
    real, allocatable, dimension(:,:,:) :: &
           pp               & !  perturbation pressure on u-grid [no units]
        ,  thvp             & !  perturbation virtual potential temperature [K]
        ,  vap2cld          & !  mixing ratio of condensed cloud water this time step (term C in HW9 equations) [kg/kg] -- rvp to rcp (or vice versa)
        ,  rain2vap         & !  evaporation rate of rain water (term E in HW9 equations) [kg/kg/s] -- rrp to rvp
        ,  cld2rain_accr    & !  accretion rate of cloud water to rain (term B in HW9 equations) [kg/kg/s] -- rcp to rrp 
        ,  cld2rain_auto      !  autoconversion rate of cloud water to rain (term A in HW9 equations) [kg/kg/s] -- rcp to rrp

    ! tendency  variables (array in nx,ny,nz)
    real, allocatable, dimension(:,:,:)  :: &
             u_xadv,u_yadv,u_zadv,u_pgf,u_xdiff,u_ydiff,u_zdiff,u_tend_total                        &
            ,v_xadv,v_yadv,v_zadv,v_pgf,v_xdiff,v_ydiff,v_zdiff,v_tend_total                        &
            ,w_xadv,w_yadv,w_zadv,w_pgf,w_buoy,w_xdiff,w_ydiff,w_zdiff,w_tend_total                 &
            ,thp_xadv,thp_yadv,thp_zadv,thp_meanadv,thp_xdiff,thp_ydiff,thp_zdiff,thp_tend_total    &
            ,pip_xadv,pip_yadv,pip_zadv,pip_xdiff,pip_ydiff,pip_zdiff,pip_tend_total                &
            ,rvp_xadv,rvp_yadv,rvp_zadv,rvp_meanadv,rvp_xdiff,rvp_ydiff,rvp_zdiff,rvp_tend_total    &
            ,rcp_xadv,rcp_yadv,rcp_zadv,rcp_xdiff,rcp_ydiff,rcp_zdiff,rcp_tend_total                &
            ,rrp_xadv,rrp_yadv,rrp_zadv,rrp_xdiff,rrp_ydiff,rrp_zdiff,rrp_tend_total                                 

    contains

end module model_vars
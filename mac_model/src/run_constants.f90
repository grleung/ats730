module run_constants

    ! This module is where to set grid constants (number of grid points, spacings, etc.)
    ! as well as any constants for the base state for this run (e.g., height of tropopause, etc.)
    ! Values are from HW1
    
    implicit none

    integer :: &
             nz      =    29       & ! number of vertical (z) levels
        ,    nx      =    29     ! number of horizontal (x) points

    real :: & ! grid parameters
            dz0       =   700.    &  ! depth of first vertical level [m]
        ,   dzrat     =   1.      &  ! ratio of subsequent vertical levels, for now keep all constant
        ,   dx        =   1000.   &  ! width of horizontal grid spacing [m]
        ,   dt        =   3.      &  ! timestep [s]
        ,   endt      =   300.       ! total integration time

    
    logical :: &
            pbc_x     = .True.    & ! use periodic boundaries for x?
        ,   pbc_z     = .True.      ! use periodic boundaries for z?
        
    ! base state settings
    real :: & 
            ztr       =   12000.  & ! height of tropopause [m]
        ,   ttr       =   213.    & ! temperature of tropopause [K]
        ,   thtr      =   343.    & ! theta of tropopause [K]
        ,   psurf     =   96500.   ! surface pressure [Pa]

    logical :: &
            wk_flag   = .False.     & ! flag for Weisman-Klempt sounding
        ,   dn_flag   = .True.        ! flag for dry and neutral environment
    
    ! parcel  settings
    real :: & 
           rvpsurf   =   11.5E-3   ! parcel water vapor mixing ratio at first real level [kg/kg]


    !output settings
    character(len=80) :: & !output paths
            base_outpath = 'hw1_output.txt' &
        ,   parcel_outpath = 'hw2_output.txt' &
        ,   var_outpath = 'hw3_output.txt'   

    logical :: &
            base_out = .False.   &
        ,   parcel_out = .False. &
        ,   var_out = .False.

    integer :: &
            outfreq = 200

    ! perturbation settings
    real :: &
            pert_wind=.True.&   ! which variable to perturb, if true then U, if false then THETA
        ,   radx = 0.       &   ! horizontal radius of perturbation [m]
        ,   radz = 0.       &   ! vertical radius of perturbation [m]
        ,   amp = 0.        &   ! thermal amplitude [K]
        ,   zcnt = 0.       &   ! center of thermal [m above ground]
        ,   xcnt = 0.       &   ! center of termal [m from W side of domain]
        ,   cx = 0.         &   ! horizontal advection speed [m/s]
        ,   cz = 0.         &    ! vertival advection speed [m/s
        ,   cs      =    50.0          ! speed of sounds [m/s] this is too slow but this is given in HW4

    ! diffusion settings
    real :: &
            kmx =   1.0     &   ! momentum exchange coefficient in x-dimension [m2/s]
        ,   kmz =   1.0     &   ! momentum exchange coefficient in z-dimension [m2/s]
        ,   khx =   1.0     &   ! scalar exchange coefficient in x-dimension [m2/s]
        ,   khz =   1.0         ! scalar exchange coefficient in z-dimension [m2/s]
        
end module run_constants
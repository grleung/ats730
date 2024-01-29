module grid_constants

    ! This module is where to set grid constants (number of grid points, spacings, etc.)
    ! as well as any constants for the base state for this run (e.g., height of tropopause, etc.)
    ! Values are from HW1
    
    implicit none

    integer :: &
            nz      =    40     ! number of vertical (z) levels

    real :: & ! grid parameters
            dz0       =   700.    & ! depth of first vertical level [m]
        ,   dzrat     =   1.       ! ratio of subsequent vertical levels, for now keep all constant
     
    real :: & ! base state settings
            ztr       =   12000.  & ! height of tropopause [m]
        ,   ttr       =   213.    & ! temperature of tropopause [K]
        ,   thtr      =   343.    & ! theta of tropopause [K]
        ,   psurf     =   96500.   ! surface pressure [Pa]
    
    real :: & ! parcel  settings
           rvpsurf   =   11.5E-3   ! parcel water vapor mixing ratio at first real level [kg/kg]

    character(len=80) :: & !output paths
            base_outpath = 'hw1_output.txt' &
        ,   parcel_outpath = 'hw2_output.txt'

    logical :: &
            base_out = .False. &
        ,   parcel_out = .False.
        
end module grid_constants
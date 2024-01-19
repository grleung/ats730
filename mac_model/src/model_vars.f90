module model_vars

    use grid_constants, only: nz

    ! This module is where we will declare all the shared arrays 
    
    implicit none

    ! grid coordinate variables
    real, dimension(nz) :: dzn ! deltaZ (level depth) of the w ("momentum") grid [m]
    real, dimension(nz) :: zun ! physical height of u/scalar ("u" for u) grid [m]
    real, dimension(nz) :: zwn ! physical height of vertical velocity ("w" for w) grid [m]

    ! base state thermodynamic variables
    real, dimension(nz) :: tb ! base state temperature ("t" for temperature) [K]
    real, dimension(nz) :: thb ! base state potential temperature ("th" for theta) [K]
    real, dimension(nz) :: rvb ! base state water vapor mixing ratio ("r" for ratio, "v" for vapor) [kg/kg]
    real, dimension(nz) :: thvb ! base state virtual potential temperature [K]
    real, dimension(nz) :: pb ! base state pressure [Pa]
    real, dimension(nz) :: pib ! base state non-dimensional pressure [no units]
    real, dimension(nz) :: rhoub ! base state air density at u/scalar levels [kg/m3]
    real, dimension(nz) :: rhowb ! base state air density at w levels[kg/m3]
    real, dimension(nz) :: rhb ! base state relative humidity [%]
    
    

end module model_vars
module io

    ! This module driver contains subroutines for writing output

    implicit none

    contains

    subroutine write_base_state (filename)
        use grid_constants, only: nz
        use model_vars, only: zun, tb, thb, rvb, thvb, pib, piwb, pb, rhb, rhoub,rhowb

        implicit none
        
        character(len=*) :: filename ! the subroutine takes the filename as an argument
        integer :: iz ! counter for z-coordinate

        ! open the output file we want to write to
        open(unit = 1, file=filename)

        write(1,'(11A10)') 'k','ht','temp','theta','thetav','pres','pi','rhou','rhow','rv','RH'
        write(1,'(11A10)') ' ','[km]','[K]','[K]','[K]','[hPa]','','[kg/m3]','[kg/m3]','[g/kg]','[%]'
        
        do iz=2,nz-1
            write(1,'(1I10,5F10.2,1F10.4,2F10.5,2F10.2)') iz,zun(iz)/1000,tb(iz),thb(iz),thvb(iz),pb(iz)/100,pib(iz),rhoub(iz),rhowb(iz),rvb(iz)*1000,rhb(iz)
        enddo

        close(1)
        

    end subroutine write_base_state

    subroutine write_parcel_traj (filename)
        use grid_constants, only: nz
        use model_vars, only: zun, rvp, thp, thvp, thvdiff, lclp, elp, capep

        implicit none
        
        character(len=*) :: filename ! the subroutine takes the filename as an argument
        integer :: iz ! counter for z-coordinate

        ! open the output file we want to write to
        open(unit = 1, file=filename)

        write(1,'(1A50,1F7.2,1A7,1F7.2,1A5)') "Initial parcel potential temperature at", zun(2)/1000., 'km is: ', thp(2), 'K'
        write(1,'(1A50,1F7.2,1A7,1F7.1,1A5)') "Initial parcel water vapor mixing ratio at", zun(2)/1000., 'km is: ', rvp(2)*1000, 'g/kg'
        write(1,'(1A20,1F7.2,1A7,1F7.2,1A5)') "LFC falls between", zun(lclp)/1000., 'km and', zun(lclp+1)/1000., 'km'
        write(1,'(1A20,1F7.2,1A7,1F7.2,1A5)') "EL falls between", zun(elp)/1000., 'km and', zun(elp+1)/1000., 'km'
        write(1,*) "CAPE is", capep, 'J/kg'

        write(1,'(6A10)') 'k','ht','rvP','thP','thvP','thv_diff'
        write(1,'(6A10)') ' ','[km]','[kg/kg]','[K]','[K]','[K]'
        
        do iz=2,nz-1
            write(1,'(1I10,5F10.2)') iz, zun(iz)/1000, rvp(iz)*1000, thp(iz), thvp(iz), thvdiff(iz)
        enddo

        close(1)

    end subroutine write_parcel_traj

end module io
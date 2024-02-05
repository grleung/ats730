module io

    ! This module driver contains subroutines for reading namelists and writing output

    implicit none

    contains

    subroutine read_namelist 
        ! so far only need these values in namelist, but will probably need more later on
        use run_constants, only: nz,dz0,nx,dx,rvpsurf,base_out,base_outpath,parcel_out,parcel_outpath,wk_flag,dn_flag

        implicit none 

        ! define namelists
        namelist /output/ base_out,base_outpath,parcel_out,parcel_outpath
        namelist /grid/ nz,nx,dz0,dx
        namelist /base/ wk_flag,dn_flag
        namelist /parcel/ rvpsurf
        
        open(unit=1, file='Namelist',action='read')

        !read each namelist group
        read(1,nml=output)
        read(1,nml=parcel)
        read(1,nml=grid)

        close(1)

    end subroutine read_namelist

    subroutine write_base_state
        use run_constants, only: nz, base_out, base_outpath
        use model_vars, only: zsn, tb, thb, rvb, thvb, pib, piwb, pb, rhb, rhoub,rhowb

        implicit none
        
        integer :: iz ! counter for z-coordinate

        if (base_out) then
            ! open the output file we want to write to
            open(unit = 1, file=base_outpath)

            write(1,'(11A10)') 'k','ht','temp','theta','thetav','pres','pi','rhou','rhow','rv','RH'
            write(1,'(11A10)') ' ','[km]','[K]','[K]','[K]','[hPa]','','[kg/m3]','[kg/m3]','[g/kg]','[%]'
            
            do iz=2,nz-1
                write(1,'(1I10,5F10.2,1F10.4,2F10.5,2F10.2)') iz,zsn(iz)/1000,tb(2,iz),thb(2,iz),thvb(2,iz),pb(2,iz)/100,pib(2,iz),rhoub(2,iz),rhowb(2,iz),rvb(2,iz)*1000,rhb(2,iz)
            enddo

            close(1)
        endif 

    end subroutine write_base_state

    subroutine write_parcel_traj 
        use run_constants, only: nz, parcel_out, parcel_outpath
        use model_vars, only: zsn, rvp, thp, thvp, thvdiff, lclp, elp, capep

        implicit none
        
        integer :: iz ! counter for z-coordinate

        if (parcel_out) then 
            ! open the output file we want to write to
            open(unit = 1, file=parcel_outpath)

            write(1,'(1A50,1F7.2,1A7,1F7.2,1A5)') "Initial parcel potential temperature at", zsn(2)/1000., 'km is: ', thp(2), 'K'
            write(1,'(1A50,1F7.2,1A7,1F7.1,1A5)') "Initial parcel water vapor mixing ratio at", zsn(2)/1000., 'km is: ', rvp(2)*1000, 'g/kg'
            write(1,'(1A20,1F7.2,1A7,1F7.2,1A5)') "LFC falls between", zsn(lclp)/1000., 'km and', zsn(lclp+1)/1000., 'km'
            write(1,'(1A20,1F7.2,1A7,1F7.2,1A5)') "EL falls between", zsn(elp)/1000., 'km and', zsn(elp+1)/1000., 'km'
            write(1,*) "CAPE is", capep, 'J/kg'

            write(1,'(6A10)') 'k','ht','rvP','thP','thvP','thv_diff'
            write(1,'(6A10)') ' ','[km]','[kg/kg]','[K]','[K]','[K]'
            
            do iz=2,nz-1
                write(1,'(1I10,5F10.2)') iz, zsn(iz)/1000, rvp(iz)*1000, thp(iz), thvp(iz), thvdiff(iz)
            enddo

            close(1)
        endif

    end subroutine write_parcel_traj

end module io
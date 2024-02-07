module io

    ! This module driver contains subroutines for reading namelists and writing output

    implicit none

    contains

    subroutine read_namelist 
        ! so far only need these values in namelist, but will probably need more later on
        use run_constants, only: nz,dz0,nx,dx,pbc,dt,rvpsurf &
                                ,base_out,base_outpath,parcel_out,parcel_outpath,var_out,var_outpath &
                                ,wk_flag,dn_flag   &
                                ,pert_wind,radx,radz,amp,zcnt,xcnt,cx,cz

        implicit none 

        ! define namelists
        namelist /output/ base_out,base_outpath,parcel_out,parcel_outpath,var_out,var_outpath
        namelist /grid/ nz,nx,dz0,dx,pbc,dt
        namelist /base/ wk_flag,dn_flag
        namelist /parcel/ rvpsurf
        namelist /pert/ pert_wind,radx,radz,amp,zcnt,xcnt,cx,cz
        
        open(unit=1, file='Namelist',action='read')

        !read each namelist group
        read(1,nml=output)
        read(1,nml=grid)
        read(1,nml=base)
        read(1,nml=parcel)
        read(1,nml=pert)
        
        close(1)

    end subroutine read_namelist

    subroutine write_base_state
        use run_constants, only: nz, base_out, base_outpath
        use model_vars, only: zsn, thb, rvb, thvb, pib, piwb, rhoub,rhowb, tb, pb,rhb

        implicit none
        
        integer :: iz ! counter for z-coordinate

        if (base_out) then
            ! open the output file we want to write to
            open(unit = 1, file=base_outpath)

            write(1,'(11A10)') 'k','ht','temp','theta','thetav','pres','pi','rhou','rhow','rv','RH'
            write(1,'(11A10)') ' ','[km]','[K]','[K]','[K]','[hPa]','','[kg/m3]','[kg/m3]','[g/kg]','[%]'
            
            do iz=2,nz-1
                write(1,'(1I10,5F10.2,1F10.4,2F10.5,2F10.2)') iz,zsn(iz)/1000,tb(iz),thb(iz),thvb(iz),pb(iz)/100,pib(iz),rhoub(iz),rhowb(iz),rvb(iz)*1000,rhb(iz)
            enddo

            close(1)
        endif 

    end subroutine write_base_state

    subroutine write_parcel_traj 
        use run_constants, only: nz, parcel_out, parcel_outpath
        use model_vars, only: zsn, rvpcl, thpcl, thvpcl, thvdiff, lclp, elp, capep

        implicit none
        
        integer :: iz ! counter for z-coordinate

        if (parcel_out) then 
            ! open the output file we want to write to
            open(unit = 1, file=parcel_outpath)

            write(1,'(1A50,1F7.2,1A7,1F7.2,1A5)') "Initial parcel potential temperature at", zsn(2)/1000., 'km is: ', thpcl(2), 'K'
            write(1,'(1A50,1F7.2,1A7,1F7.1,1A5)') "Initial parcel water vapor mixing ratio at", zsn(2)/1000., 'km is: ', rvpcl(2)*1000, 'g/kg'
            write(1,'(1A20,1F7.2,1A7,1F7.2,1A5)') "LFC falls between", zsn(lclp)/1000., 'km and', zsn(lclp+1)/1000., 'km'
            write(1,'(1A20,1F7.2,1A7,1F7.2,1A5)') "EL falls between", zsn(elp)/1000., 'km and', zsn(elp+1)/1000., 'km'
            write(1,*) "CAPE is", capep, 'J/kg'

            write(1,'(6A10)') 'k','ht','rvP','thP','thvP','thv_diff'
            write(1,'(6A10)') ' ','[km]','[kg/kg]','[K]','[K]','[K]'
            
            do iz=2,nz-1
                write(1,'(1I10,5F10.2)') iz, zsn(iz)/1000, rvpcl(iz)*1000, thpcl(iz), thvpcl(iz), thvdiff(iz)
            enddo

            close(1)
        endif

    end subroutine write_parcel_traj

    subroutine write_current_state
        use run_constants, only: nz, var_out, var_outpath
        use model_vars, only: it,zsn,xsn,thp,pip,up,wp,pp

        implicit none
        
        integer :: iz ! counter for z-coordinate

        character(len=3) :: timechar
        write(timechar, '(i3)')it

        if (var_out) then
            ! open the output file we want to write to
            open(unit = 1, file=trim(var_outpath)//'timestep_'//trim(adjustl(timechar))//'.txt')

            write(1,*) 'coord zsn'
            write(1, '(1x, *(g0, :, ", "))') zsn(:)

            write(1,*) 'coord xsn'
            write(1, '(1x, *(g0, :, ", "))') xsn(:)

            write(1,*) 'var THP_past'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') thp(iz,:,1)
            enddo

            write(1,*) 'var THP_pres'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') thp(iz,:,2)
            enddo

            write(1,*) 'var THP_fut'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') thp(iz,:,3)
            enddo

            write(1,*) 'var PIP_past'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') pip(iz,:,1)
            enddo

            write(1,*) 'var PIP_pres'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') pip(iz,:,2)
            enddo

            write(1,*) 'var PIP_fut'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') pip(iz,:,3)
            enddo

            write(1,*) 'var UP_past'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') up(iz,:,1)
            enddo

            write(1,*) 'var UP_pres'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') up(iz,:,2)
            enddo

            write(1,*) 'var UP_fut'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') up(iz,:,3)
            enddo

            write(1,*) 'var WP_past'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') wp(iz,:,1)
            enddo

            write(1,*) 'var WP_pres'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') wp(iz,:,2)
            enddo

            write(1,*) 'var WP_fut'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') wp(iz,:,3)
            enddo

            write(1,*) 'var PP_pres'
            do iz=1,nz
                write(1, '(1x, *(g0, :, ", "))') pp(iz,:)
            enddo

            close(1)
        endif 

    end subroutine write_current_state

end module io
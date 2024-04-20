module timestep
    ! This module contains the subroutine to step model vars forward in time

    implicit none

    contains
    
    subroutine step_time
        use model_vars, only: thp,pip,pp,up,wp,rvp,np,mp,nc,mc,mpc,nr,mr,mpr
        use run_constants, only: nz,nx,npb,ndb
        use aerosol, only: check_negs

        implicit none

        integer :: iz ! counter for z
        integer :: ix ! counter for x
        integer :: ipb,idb,it


        print*,'step time start'
        ! step forward in time
        do ix = 1, nx
            do iz = 1, nz
                !do 1 = 1,2
                up(ix,iz,1) = up(ix,iz,2)
                wp(ix,iz,1) = wp(ix,iz,2)
                thp(ix,iz,1) = thp(ix,iz,2)
                pip(ix,iz,1) = pip(ix,iz,2)
                rvp(ix,iz,1) = rvp(ix,iz,2)

                do ipb=1,npb
                    np(ix,iz,ipb,1) = np(ix,iz,ipb,2)
                    mp(ix,iz,ipb,1) = mp(ix,iz,ipb,2)
                    
                    do idb=1,ndb
                        ! smth is happening here idk what but after setting current there's smth wrong
                        mc(ix,iz,ipb,idb,1) = mc(ix,iz,ipb,idb,2)
                        mpc(ix,iz,ipb,idb,1) = mpc(ix,iz,ipb,idb,2)
                        nr(ix,iz,ipb,idb,1) = nr(ix,iz,ipb,idb,2)
                        mr(ix,iz,ipb,idb,1) = mr(ix,iz,ipb,idb,2)
                        mpr(ix,iz,ipb,idb,1) = mpr(ix,iz,ipb,idb,2)
                        nc(ix,iz,ipb,idb,1) = nc(ix,iz,ipb,idb,2)
                    enddo
                enddo

                up(ix,iz,2) = up(ix,iz,3)
                wp(ix,iz,2) = wp(ix,iz,3)
                thp(ix,iz,2) = thp(ix,iz,3)
                pip(ix,iz,2) = pip(ix,iz,3)
                rvp(ix,iz,2) = rvp(ix,iz,3)

                do ipb=1,npb
                    np(ix,iz,ipb,2) = np(ix,iz,ipb,3)
                    mp(ix,iz,ipb,2) = mp(ix,iz,ipb,3)
                    
                    do idb=1,ndb
                        ! smth is happening here idk what but after setting current there's smth wrong
                        mc(ix,iz,ipb,idb,2) = mc(ix,iz,ipb,idb,3)
                        mpc(ix,iz,ipb,idb,2) = mpc(ix,iz,ipb,idb,3)
                        nr(ix,iz,ipb,idb,2) = nr(ix,iz,ipb,idb,3)
                        mr(ix,iz,ipb,idb,2) = mr(ix,iz,ipb,idb,3)
                        mpr(ix,iz,ipb,idb,2) = mpr(ix,iz,ipb,idb,3)

                        if (nc(ix,iz,ipb,idb,3)<0.) then
                            print*,''
                            nc(ix,iz,ipb,idb,2) = 0.
                        else
                            nc(ix,iz,ipb,idb,2) = nc(ix,iz,ipb,idb,3)
                        endif 
                    enddo
                enddo
                !enddo
            enddo ! end z loop
        enddo ! end x loop

        call check_negs


        print*,MAXVAL(wp(:,:,2))
        print*,MAXVAL(rvp(:,:,2))
        !print*,SUM(np(:,:,:,2))+SUM(nc(:,:,:,:,2))
        
    end subroutine step_time

end module timestep
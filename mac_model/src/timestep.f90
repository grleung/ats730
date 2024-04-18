module timestep
    ! This module contains the subroutine to step model vars forward in time

    implicit none

    contains
    
    subroutine step_time
        use model_vars, only: thp,pip,pp,up,wp,rvp,np,mp,nc,mc,mpc,nr,mr,mpr
        use run_constants, only: nz,nx,npb,ndb

        implicit none

        integer :: iz ! counter for z
        integer :: ix ! counter for x
        integer :: iab,idb,it

        ! step forward in time
        do ix = 1, nx
            do iz = 1, nz
                do it = 1,2
                    up(ix,iz,it) = up(ix,iz,it+1)
                    wp(ix,iz,it) = wp(ix,iz,it+1)
                    thp(ix,iz,it) = thp(ix,iz,it+1)
                    pip(ix,iz,it) = pip(ix,iz,it+1)
                    rvp(ix,iz,it) = rvp(ix,iz,it+1)

                    do iab=1,npb
                        np(ix,iz,iab,it) = np(ix,iz,iab,it+1)
                        mp(ix,iz,iab,it) = mp(ix,iz,iab,it+1)
                        
                        do idb=1,ndb
                            nc(ix,iz,iab,idb,it) = nc(ix,iz,iab,idb,it+1)
                            mc(ix,iz,iab,idb,it) = mc(ix,iz,iab,idb,it+1)
                            mpc(ix,iz,iab,idb,it) = mpc(ix,iz,iab,idb,it+1)
                            nr(ix,iz,iab,idb,it) = nr(ix,iz,iab,idb,it+1)
                            mr(ix,iz,iab,idb,it) = mr(ix,iz,iab,idb,it+1)
                            mpr(ix,iz,iab,idb,it) = mpr(ix,iz,iab,idb,it+1)
                        enddo
                    enddo
                enddo
            enddo ! end z loop
        enddo ! end x loop
        
    end subroutine step_time

end module timestep
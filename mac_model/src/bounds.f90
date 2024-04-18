module boundaries
    ! This module contains the subroutines to enforce boundary conditions

    implicit none

    contains
    
    subroutine enforce_bounds_x
        use model_vars, only: thp,pip,pp,up,wp,rvp,np,mp,nc,mc,mpc,nr,mr,mpr
        use run_constants, only: nz,nx,pbc_x,npb,ndb

        implicit none

        integer :: iz ! counter for z
        integer :: ix ! counter for x
        integer :: it ! counter for time
        integer :: iab,idb

        if (pbc_x==.True.) then
            ! if we use periodic lateral boundaries, then we need to set all ix=1 to ix=nx-1 and all ix=nx to ix=2
            do iz=1,nz
                do it=1,3
                    pip(1,iz,it) = pip(nx-1,iz,it)
                    pip(nx,iz,it) = pip(2,iz,it)

                    thp(1,iz,it) = thp(nx-1,iz,it)
                    thp(nx,iz,it) = thp(2,iz,it)

                    up(1,iz,it) = up(nx-1,iz,it)
                    up(nx,iz,it) = up(2,iz,it)

                    wp(1,iz,it) = wp(nx-1,iz,it)
                    wp(nx,iz,it) = wp(2,iz,it)

                    rvp(1,iz,it) = rvp(nx-1,iz,it)
                    rvp(nx,iz,it) = rvp(2,iz,it)

                    do iab=1,npb
                        !np,mp,nc,mc,npc,mpc,nr,mr,mpr
                        np(1,iz,iab,it) = np(nx-1,iz,iab,it)
                        np(nx,iz,iab,it) = np(2,iz,iab,it)

                        mp(1,iz,iab,it) = mp(nx-1,iz,iab,it)
                        mp(nx,iz,iab,it) = mp(2,iz,iab,it)

                        do idb=1,ndb
                            nc(1,iz,iab,idb,it) = nc(nx-1,iz,iab,idb,it)
                            nc(nx,iz,iab,idb,it) = nc(2,iz,iab,idb,it)

                            mc(1,iz,iab,idb,it) = mc(nx-1,iz,iab,idb,it)
                            mc(nx,iz,iab,idb,it) = mc(2,iz,iab,idb,it)

                            mpc(1,iz,iab,idb,it) = mpc(nx-1,iz,iab,idb,it)
                            mpc(nx,iz,iab,idb,it) = mpc(2,iz,iab,idb,it)

                            nr(1,iz,iab,idb,it) = nr(nx-1,iz,iab,idb,it)
                            nr(nx,iz,iab,idb,it) = nr(2,iz,iab,idb,it)

                            mr(1,iz,iab,idb,it) = mr(nx-1,iz,iab,idb,it)
                            mr(nx,iz,iab,idb,it) = mr(2,iz,iab,idb,it)

                            mpr(1,iz,iab,idb,it) = mpr(nx-1,iz,iab,idb,it)
                            mpr(nx,iz,iab,idb,it) = mpr(2,iz,iab,idb,it)
                        enddo
                    enddo 
                enddo
            enddo ! end loop over z
        else
            !if not using periodic lateral boundaries, enforce zero gradient at edges
            do iz=1,nz
                do it=1,3
                    pip(1,iz,it) = pip(2,iz,it)
                    pip(nx,iz,it) = pip(nx-1,iz,it)

                    thp(1,iz,it) = thp(2,iz,it)
                    thp(nx,iz,it) = thp(nx-1,iz,it)

                    up(1,iz,it) = up(2,iz,it)
                    up(nx,iz,it) = up(nx-1,iz,it)

                    wp(1,iz,it) = wp(2,iz,it)
                    wp(nx,iz,it) = wp(nx-1,iz,it)
                    
                    rvp(1,iz,it) = rvp(2,iz,it)
                    rvp(nx,iz,it) = rvp(nx-1,iz,it)

                    do iab=1,npb
                        !np,mp,nc,mc,npc,mpc,nr,mr,mpr
                        np(1,iz,iab,it) = np(2,iz,iab,it)
                        np(nx,iz,iab,it) = np(nx-1,iz,iab,it)

                        mp(1,iz,iab,it) = mp(2,iz,iab,it)
                        mp(nx,iz,iab,it) = mp(nx-1,iz,iab,it)

                        do idb=1,ndb
                            nc(1,iz,iab,idb,it) = nc(2,iz,iab,idb,it)
                            nc(nx,iz,iab,idb,it) = nc(nx-1,iz,iab,idb,it)

                            mc(1,iz,iab,idb,it) = mc(2,iz,iab,idb,it)
                            mc(nx,iz,iab,idb,it) = mc(nx-1,iz,iab,idb,it)

                            mpc(1,iz,iab,idb,it) = mpc(2,iz,iab,idb,it)
                            mpc(nx,iz,iab,idb,it) = mpc(nx-1,iz,iab,idb,it)

                            nr(1,iz,iab,idb,it) = nr(2,iz,iab,idb,it)
                            nr(nx,iz,iab,idb,it) = nr(nx-1,iz,iab,idb,it)

                            mr(1,iz,iab,idb,it) = mr(2,iz,iab,idb,it)
                            mr(nx,iz,iab,idb,it) = mr(nx-1,iz,iab,idb,it)

                            mpr(1,iz,iab,idb,it) = mpr(2,iz,iab,idb,it)
                            mpr(nx,iz,iab,idb,it) = mpr(nx-1,iz,iab,idb,it)
                        enddo
                    enddo 
                enddo
            enddo ! end loop over z
        endif ! end PBC flag

    end subroutine enforce_bounds_x

    subroutine enforce_bounds_z
        use model_vars, only: thp,pip,pp,up,wp,rvp,np,mp,nc,mc,mpc,nr,mr,mpr
        use run_constants, only: nz,nx,pbc_z,npb,ndb

        implicit none

        integer :: iz ! counter for z
        integer :: ix ! counter for x
        integer :: it ! counter for time
        integer :: iab,idb

        if (pbc_z==.True.) then
            ! if we use periodic lateral boundaries, then we need to set all ix=1 to ix=nx-1 and all ix=nx to ix=2
            do ix=2,nx-1
                do it=1,3
                    pip(ix,1,it) = pip(ix,nz-1,it)
                    pip(ix,nz,it) = pip(ix,2,it)

                    thp(ix,1,it) = thp(ix,nz-1,it)
                    thp(ix,nz,it) = thp(ix,2,it)

                    up(ix,1,it) = up(ix,nz-1,it)
                    up(ix,nz,it) = up(ix,2,it)

                    wp(ix,1,it) = wp(ix,nz-1,it)
                    wp(ix,nz,it) = wp(ix,2,it)

                    rvp(ix,1,it) = rvp(ix,nz-1,it)
                    rvp(ix,nz,it) = rvp(ix,2,it)
                    
                    do iab=1,npb
                        !np,mp,nc,mc,npc,mpc,nr,mr,mpr
                        np(ix,1,iab,it) = np(ix,nz-1,iab,it)
                        np(ix,nz,iab,it) = np(ix,2,iab,it)

                        mp(ix,1,iab,it) = mp(ix,nz-1,iab,it)
                        mp(ix,nz,iab,it) = mp(ix,2,iab,it)

                        do idb=1,ndb
                            nc(ix,1,iab,idb,it) = nc(ix,nz-1,iab,idb,it)
                            nc(ix,nz,iab,idb,it) = nc(ix,2,iab,idb,it)

                            mc(ix,1,iab,idb,it) = mc(ix,nz-1,iab,idb,it)
                            mc(ix,nz,iab,idb,it) = mc(ix,2,iab,idb,it)

                            mpc(ix,1,iab,idb,it) = mpc(ix,nz-1,iab,idb,it)
                            mpc(ix,nz,iab,idb,it) = mpc(ix,2,iab,idb,it)

                            nr(ix,1,iab,idb,it) = nr(ix,nz-1,iab,idb,it)
                            nr(ix,nz,iab,idb,it) = nr(ix,2,iab,idb,it)

                            mr(ix,1,iab,idb,it) = mr(ix,nz-1,iab,idb,it)
                            mr(ix,nz,iab,idb,it) = mr(ix,2,iab,idb,it)

                            mpr(ix,1,iab,idb,it) = mpr(ix,nz-1,iab,idb,it)
                            mpr(ix,nz,iab,idb,it) = mpr(ix,2,iab,idb,it)
                        enddo
                    enddo 
                enddo
            enddo ! end loop over x
        else
            !if not using periodic lateral boundaries, enforce zero gradient at top and bottom + zero for w
            do ix=2,nx-1
                do it=1,3
                    pip(ix,1,it) = pip(ix,2,it)
                    pip(ix,nz,it) = pip(ix,nz-1,it)

                    thp(ix,1,it) = thp(ix,2,it)
                    thp(ix,nz,it) = thp(ix,nz-1,it)

                    up(ix,1,it) = up(ix,2,it)
                    up(ix,nz,it) = up(ix,nz-1,it)

                    !enforce zero w through top and bottom
                    wp(ix,1,it) = 0.
                    wp(ix,2,it) = 0.
                    wp(ix,nz,it) = 0.

                    rvp(ix,1,it) = rvp(ix,2,it)
                    rvp(ix,nz,it) = rvp(ix,nz-1,it)

                    do iab=1,npb
                        !np,mp,nc,mc,npc,mpc,nr,mr,mpr
                        np(ix,1,iab,it) = np(ix,2,iab,it)
                        np(ix,nz,iab,it) = np(ix,nz-1,iab,it)

                        mp(ix,1,iab,it) = mp(ix,2,iab,it)
                        mp(ix,nz,iab,it) = mp(ix,nz-1,iab,it)

                        do idb=1,ndb
                            nc(ix,1,iab,idb,it) = nc(ix,2,iab,idb,it)
                            nc(ix,nz,iab,idb,it) = nc(ix,nz-1,iab,idb,it)

                            mc(ix,1,iab,idb,it) = mc(ix,2,iab,idb,it)
                            mc(ix,nz,iab,idb,it) = mc(ix,nz-1,iab,idb,it)

                            mpc(ix,1,iab,idb,it) = mpc(ix,2,iab,idb,it)
                            mpc(ix,nz,iab,idb,it) = mpc(ix,nz-1,iab,idb,it)

                            nr(ix,1,iab,idb,it) = nr(ix,2,iab,idb,it)
                            nr(ix,nz,iab,idb,it) = nr(ix,nz-1,iab,idb,it)

                            mr(ix,1,iab,idb,it) = mr(ix,2,iab,idb,it)
                            mr(ix,nz,iab,idb,it) = mr(ix,nz-1,iab,idb,it)

                            mpr(ix,1,iab,idb,it) = mpr(ix,2,iab,idb,it)
                            mpr(ix,nz,iab,idb,it) = mpr(ix,nz-1,iab,idb,it)
                        enddo
                    enddo 
                enddo
            enddo ! end loop over x
        endif ! end PBC flag

    end subroutine enforce_bounds_z

end module boundaries
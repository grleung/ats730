module boundaries
    ! This module contains the subroutines to enforce boundary conditions

    implicit none

    contains
    
    subroutine enforce_bounds_x
        use model_vars, only: thp,pip,pp,up,wp,rvp,rcp,rrp
        use run_constants, only: nz,nx,pbc_x

        implicit none

        integer :: iz ! counter for z
        integer :: ix ! counter for x
        integer :: it ! counter for time

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

                    rcp(1,iz,it) = rcp(nx-1,iz,it)
                    rcp(nx,iz,it) = rcp(2,iz,it)

                    rrp(1,iz,it) = rrp(nx-1,iz,it)
                    rrp(nx,iz,it) = rrp(2,iz,it)
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

                    rcp(1,iz,it) = rcp(2,iz,it)
                    rcp(nx,iz,it) = rcp(nx-1,iz,it)

                    rrp(1,iz,it) = rrp(2,iz,it)
                    rrp(nx,iz,it) = rrp(nx-1,iz,it)
                enddo
            enddo ! end loop over z
        endif ! end PBC flag

    end subroutine enforce_bounds_x

    subroutine enforce_bounds_z
        use model_vars, only: thp,pip,pp,up,wp,rcp,rrp,rvp
        use run_constants, only: nz,nx,pbc_z

        implicit none

        integer :: iz ! counter for z
        integer :: ix ! counter for x
        integer :: it ! counter for time

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
                    
                    rcp(ix,1,it) = rcp(ix,nz-1,it)
                    rcp(ix,nz,it) = rcp(ix,2,it)
                    
                    rrp(ix,1,it) = rrp(ix,nz-1,it)
                    rrp(ix,nz,it) = rrp(ix,2,it)
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

                    rcp(ix,1,it) = rcp(ix,2,it)
                    rcp(ix,nz,it) = rcp(ix,nz-1,it)

                    rrp(ix,1,it) = rrp(ix,2,it)
                    rrp(ix,nz,it) = rrp(ix,nz-1,it)
                enddo
            enddo ! end loop over x
        endif ! end PBC flag

    end subroutine enforce_bounds_z

end module boundaries
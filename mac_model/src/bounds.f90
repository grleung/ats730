module boundaries
    ! This module contains the subroutines to enforce boundary conditions

    implicit none

    contains
    
    subroutine enforce_bounds_x
        use model_vars, only: thp,pip,up,vp,wp,rvp,rcp,rrp
        use run_constants, only: nz,nx,ny,pbc_x

        implicit none

        integer :: iz ! counter for z
        integer :: iy ! counter for y
        integer :: it ! counter for time

        if (pbc_x==.True.) then
            ! if we use periodic lateral boundaries, then we need to set all ix=1 to ix=nx-1 and all ix=nx to ix=2
            do iz=1,nz
                do iy=1,ny
                    do it=1,3
                        pip(1,iy,iz,it) = pip(nx-1,iy,iz,it)
                        pip(nx,iy,iz,it) = pip(2,iy,iz,it)

                        thp(1,iy,iz,it) = thp(nx-1,iy,iz,it)
                        thp(nx,iy,iz,it) = thp(2,iy,iz,it)

                        up(1,iy,iz,it) = up(nx-1,iy,iz,it)
                        up(nx,iy,iz,it) = up(2,iy,iz,it)

                        vp(1,iy,iz,it) = vp(nx-1,iy,iz,it)
                        vp(nx,iy,iz,it) = vp(2,iy,iz,it)

                        wp(1,iy,iz,it) = wp(nx-1,iy,iz,it)
                        wp(nx,iy,iz,it) = wp(2,iy,iz,it)

                        rvp(1,iy,iz,it) = rvp(nx-1,iy,iz,it)
                        rvp(nx,iy,iz,it) = rvp(2,iy,iz,it)

                        rcp(1,iy,iz,it) = rcp(nx-1,iy,iz,it)
                        rcp(nx,iy,iz,it) = rcp(2,iy,iz,it)

                        rrp(1,iy,iz,it) = rrp(nx-1,iy,iz,it)
                        rrp(nx,iy,iz,it) = rrp(2,iy,iz,it)
                    enddo ! end loop over time
                enddo ! end loop over y
            enddo ! end loop over z
        else
            !if not using periodic lateral boundaries, enforce zero gradient at edges
            do iz=1,nz
                do iy=1,ny
                    do it=1,3
                        pip(1,iy,iz,it) = pip(2,iy,iz,it)
                        pip(nx,iy,iz,it) = pip(nx-1,iy,iz,it)

                        thp(1,iy,iz,it) = thp(2,iy,iz,it)
                        thp(nx,iy,iz,it) = thp(nx-1,iy,iz,it)

                        up(1,iy,iz,it) = up(2,iy,iz,it)
                        up(nx,iy,iz,it) = up(nx-1,iy,iz,it)

                        vp(1,iy,iz,it) = vp(2,iy,iz,it)
                        vp(nx,iy,iz,it) = vp(nx-1,iy,iz,it)

                        wp(1,iy,iz,it) = wp(2,iy,iz,it)
                        wp(nx,iy,iz,it) = wp(nx-1,iy,iz,it)

                        rvp(1,iy,iz,it) = rvp(2,iy,iz,it)
                        rvp(nx,iy,iz,it) = rvp(nx-1,iy,iz,it)

                        rcp(1,iy,iz,it) = rcp(2,iy,iz,it)
                        rcp(nx,iy,iz,it) = rcp(nx-1,iy,iz,it)

                        rrp(1,iy,iz,it) = rrp(2,iy,iz,it)
                        rrp(nx,iy,iz,it) = rrp(nx-1,iy,iz,it)
                    enddo ! end loop over time
                enddo ! end loop over y
            enddo ! end loop over z
        endif ! end PBC flag

    end subroutine enforce_bounds_x

    subroutine enforce_bounds_y
        use model_vars, only: thp,pip,up,vp,wp,rvp,rcp,rrp
        use run_constants, only: nz,nx,ny,pbc_y

        implicit none

        integer :: iz ! counter for z
        integer :: ix ! counter for x
        integer :: it ! counter for time

        if (pbc_y==.True.) then
            ! if we use periodic lateral boundaries, then we need to set all iy=1 to iy=ny-1 and all iy=ny to iy=2
            do iz=1,nz
                do ix=1,nx
                    do it=1,3
                        pip(ix,1,iz,it) = pip(ix,ny-1,iz,it)
                        pip(ix,ny,iz,it) = pip(ix,2,iz,it)

                        thp(ix,1,iz,it) = thp(ix,ny-1,iz,it)
                        thp(ix,ny,iz,it) = thp(ix,2,iz,it)

                        up(ix,1,iz,it) = up(ix,ny-1,iz,it)
                        up(ix,ny,iz,it) = up(ix,2,iz,it)

                        vp(ix,1,iz,it) = vp(ix,ny-1,iz,it)
                        vp(ix,ny,iz,it) = vp(ix,2,iz,it)

                        wp(ix,1,iz,it) = wp(ix,ny-1,iz,it)
                        wp(ix,ny,iz,it) = wp(ix,2,iz,it)
                        
                        rvp(ix,1,iz,it) = rvp(ix,ny-1,iz,it)
                        rvp(ix,ny,iz,it) = rvp(ix,2,iz,it)

                        rcp(ix,1,iz,it) = rcp(ix,ny-1,iz,it)
                        rcp(ix,ny,iz,it) = rcp(ix,2,iz,it)

                        rrp(ix,1,iz,it) = rrp(ix,ny-1,iz,it)
                        rrp(ix,ny,iz,it) = rrp(ix,2,iz,it)
                    enddo ! end loop over time
                enddo ! end loop over x
            enddo ! end loop over z
        else
            !if not using periodic lateral boundaries, enforce zero gradient at edges
            do iz=1,nz
                do ix=1,nx
                    do it=1,3
                        pip(ix,1,iz,it) = pip(ix,2,iz,it)
                        pip(ix,ny,iz,it) = pip(ix,ny-1,iz,it)

                        thp(ix,1,iz,it) = thp(ix,2,iz,it)
                        thp(ix,ny,iz,it) = thp(ix,ny-1,iz,it)

                        up(ix,1,iz,it) = up(ix,2,iz,it)
                        up(ix,ny,iz,it) = up(ix,ny-1,iz,it)

                        vp(ix,1,iz,it) = vp(ix,2,iz,it)
                        vp(ix,ny,iz,it) = vp(ix,ny-1,iz,it)

                        wp(ix,1,iz,it) = wp(ix,2,iz,it)
                        wp(ix,ny,iz,it) = wp(ix,ny-1,iz,it)

                        rvp(ix,1,iz,it) = rvp(ix,2,iz,it)
                        rvp(ix,ny,iz,it) = rvp(ix,ny-1,iz,it)

                        rcp(ix,1,iz,it) = rcp(ix,2,iz,it)
                        rcp(ix,ny,iz,it) = rcp(ix,ny-1,iz,it)

                        rrp(ix,1,iz,it) = rrp(ix,2,iz,it)
                        rrp(ix,ny,iz,it) = rrp(ix,ny-1,iz,it)
                    enddo ! end loop over time
                enddo ! end loop over x
            enddo ! end loop over z
        endif ! end PBC flag

    end subroutine enforce_bounds_y
    

    subroutine enforce_bounds_z
        use model_vars, only: thp,pip,up,vp,wp,rvp,rcp,rrp
        use run_constants, only: nz,nx,ny,pbc_z

        implicit none

        integer :: iy ! counter for y
        integer :: ix ! counter for x
        integer :: it ! counter for time

        if (pbc_z==.True.) then
            ! if we use periodic lateral boundaries, then we need to set all ix=1 to ix=nx-1 and all ix=nx to ix=2
            do iy=2,ny-1
                do ix=2,nx-1
                    do it=1,3
                        pip(ix,iy,1,it) = pip(ix,iy,nz-1,it)
                        pip(ix,iy,nz,it) = pip(ix,iy,2,it)

                        thp(ix,iy,1,it) = thp(ix,iy,nz-1,it)
                        thp(ix,iy,nz,it) = thp(ix,iy,2,it)

                        up(ix,iy,1,it) = up(ix,iy,nz-1,it)
                        up(ix,iy,nz,it) = up(ix,iy,2,it)
                        
                        vp(ix,iy,1,it) = vp(ix,iy,nz-1,it)
                        vp(ix,iy,nz,it) = vp(ix,iy,2,it)

                        wp(ix,iy,1,it) = wp(ix,iy,nz-1,it)
                        wp(ix,iy,nz,it) = wp(ix,iy,2,it)
                        
                        rvp(ix,iy,1,it) = rvp(ix,iy,nz-1,it)
                        rvp(ix,iy,nz,it) = rvp(ix,iy,2,it)
                        
                        rcp(ix,iy,1,it) = rcp(ix,iy,nz-1,it)
                        rcp(ix,iy,nz,it) = rcp(ix,iy,2,it)

                        rrp(ix,iy,1,it) = rrp(ix,iy,nz-1,it)
                        rrp(ix,iy,nz,it) = rrp(ix,iy,2,it)
                    enddo ! end loop over time
                enddo ! end loop over x
            enddo ! end loop over y
        else
            !if not using periodic lateral boundaries, enforce zero gradient at top and bottom + zero for w
            do iy=2,ny-1
                do ix=2,nx-1
                    do it=1,3
                        pip(ix,iy,1,it) = pip(ix,iy,2,it)
                        pip(ix,iy,nz,it) = pip(ix,iy,nz-1,it)

                        thp(ix,iy,1,it) = thp(ix,iy,2,it)
                        thp(ix,iy,nz,it) = thp(ix,iy,nz-1,it)

                        up(ix,iy,1,it) = up(ix,iy,2,it)
                        up(ix,iy,nz,it) = up(ix,iy,nz-1,it)

                        vp(ix,iy,1,it) = vp(ix,iy,2,it)
                        vp(ix,iy,nz,it) = vp(ix,iy,nz-1,it)
                        
                        !enforce zero w through top and bottom
                        wp(ix,iy,1,it) = 0.
                        wp(ix,iy,2,it) = 0.
                        wp(ix,iy,nz,it) = 0.

                        rvp(ix,iy,1,it) = rvp(ix,iy,2,it)
                        rvp(ix,iy,nz,it) = rvp(ix,iy,nz-1,it)
                        
                        rcp(ix,iy,1,it) = rcp(ix,iy,2,it)
                        rcp(ix,iy,nz,it) = rcp(ix,iy,nz-1,it)

                        rrp(ix,iy,1,it) = rrp(ix,iy,2,it)
                        rrp(ix,iy,nz,it) = rrp(ix,iy,nz-1,it)
                    enddo
                enddo
            enddo ! end loop over x
        endif ! end PBC flag

    end subroutine enforce_bounds_z

end module boundaries
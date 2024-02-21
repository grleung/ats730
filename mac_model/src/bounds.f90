module boundaries
    ! This module contains the subroutines to enforce boundary conditions

    implicit none

    contains
    
    subroutine enforce_bounds_x
        use model_vars, only: thp,pip,pp,up,wp
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
                enddo
            enddo ! end loop over z
        endif ! end PBC flag

    end subroutine enforce_bounds_x

    subroutine enforce_bounds_z
        use model_vars, only: thp,pip,pp,up,wp
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
                enddo
            enddo ! end loop over x
        endif ! end PBC flag

    end subroutine enforce_bounds_z

end module boundaries
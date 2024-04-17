module cfl
    !This module contains the subroutine for checking CFL condition

    implicit none

    contains

    subroutine check_cfl
        use model_vars, only: up,vp,wp
        use run_constants, only: dx,dy,dz0,nx,ny,nz,dt

        implicit none

        integer :: ix, iy, iz ! counters

        real :: cflcond = 1. ! set CFL condition to be 1, this might need to be stricter?
        real :: cflx,cfly,cflz

        do iz=1,nz
            do iy=1,ny
                do ix=1,nx
                    cflx = up(ix,iy,iz,2) * dt / dx
                    cfly = vp(ix,iy,iz,2) * dt / dy
                    cflz = wp(ix,iy,iz,2) * dt / dz0
                    
                    if (cflx>cflcond) then
                        print*,'CFL limit exceeded x',ix,iy,iz,up(ix,iy,iz,2),cflx
                        stop
                    endif 

                    if (cfly>cflcond) then
                        print*,'CFL limit exceeded y',ix,iy,iz,vp(ix,iy,iz,2),cfly
                        stop
                    endif 

                    if (cflz>cflcond) then
                        print*,'CFL limit exceeded z',ix,iy,iz,wp(ix,iy,iz,2),cflz
                        stop
                    endif 

                    
                enddo
            enddo
        enddo

    end subroutine check_cfl

end module cfl

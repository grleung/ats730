module cfl
    !This module contains the subroutine for checking CFL condition

    implicit none

    contains

    subroutine check_cfl
        use model_vars, only: up,wp
        use run_constants, only: dx,dz0,nx,nz,dt

        implicit none

        integer :: ix, iz ! counters

        real :: cflcond = 1. ! set CFL condition to be 1, this might need to be stricter?
        real :: cflx,cflz

        do iz=1,nz
                do ix=1,nx
                    cflx = up(ix,iz,2) * dt / dx
                    cflz = wp(ix,iz,2) * dt / dz0
                    
                    if (cflx>cflcond) then
                        print*,'CFL limit exceeded x',ix,iz,up(ix,iz,2),cflx
                        stop
                    endif 

                    if (cflz>cflcond) then
                        print*,'CFL limit exceeded z',ix,iz,wp(ix,iz,2),cflz
                        stop
                    endif 

                    
                enddo
        enddo

    end subroutine check_cfl

end module cfl

module timestep
    ! This module contains the subroutine to step model vars forward in time

    implicit none

    contains
    
    subroutine step_time
        use model_vars, only: thp,pip,pp,up,vp,wp
        use run_constants, only: nz,nx,ny

        implicit none

        integer :: iz ! counter for z
        integer :: iy ! counter for y
        integer :: ix ! counter for x

        ! step forward in time
        do ix = 1, nx
            do iy=1,ny
                do iz = 1, nz
                    up(ix,iy,iz,1) = up(ix,iy,iz,2)
                    up(ix,iy,iz,2) = up(ix,iy,iz,3)
                    vp(ix,iy,iz,1) = vp(ix,iy,iz,2)
                    vp(ix,iy,iz,2) = vp(ix,iy,iz,3)
                    wp(ix,iy,iz,1) = wp(ix,iy,iz,2)
                    wp(ix,iy,iz,2) = wp(ix,iy,iz,3)
                    thp(ix,iy,iz,1) = thp(ix,iy,iz,2)
                    thp(ix,iy,iz,2) = thp(ix,iy,iz,3)
                    pip(ix,iy,iz,1) = pip(ix,iy,iz,2)
                    pip(ix,iy,iz,2) = pip(ix,iy,iz,3)
                enddo ! end z loop
            enddo !end y loop
        enddo ! end x loop
        
    end subroutine step_time

end module timestep
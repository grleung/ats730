module timestep
    ! This module contains the subroutine to step model vars forward in time

    implicit none

    contains
    
    subroutine step_time
        use model_vars, only: thp,pip,pp,up,wp
        use run_constants, only: nz,nx

        implicit none

        integer :: iz ! counter for z
        integer :: ix ! counter for x

        ! step forward in time
        do ix = 1, nx
            do iz = 1, nz
                up(ix,iz,1) = up(ix,iz,2)
                up(ix,iz,2) = up(ix,iz,3)
                wp(ix,iz,1) = wp(ix,iz,2)
                wp(ix,iz,2) = wp(ix,iz,3)
                thp(ix,iz,1) = thp(ix,iz,2)
                thp(ix,iz,2) = thp(ix,iz,3)
                pip(ix,iz,1) = pip(ix,iz,2)
                pip(ix,iz,2) = pip(ix,iz,3)
            enddo ! end x loop
        enddo ! end z loop
        
    end subroutine step_time

end module timestep
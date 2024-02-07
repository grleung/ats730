module advection

    ! This module contains the subroutine to advect a perturbation

    implicit none

    contains

    subroutine advect
        use run_constants, only: nz,nx,dx,dz0,cx,cz,dt,pbc
        use model_vars, only:it,zsn,xsn,thb,thvb,rhoub,thp,pip,pp,up,wp

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: ix ! counter for x-coordinate

        real :: rd2x,rd2z,d2t
        
        rd2x      = (1/(dx+dx))  ! reciprocal of 2dx [1/m]
        rd2z      = (1/(dz0+dz0))  ! reciprocal of 2dz [1/m]
        d2t       = (dt+dt)         ! 2*dt

        

        ! loop over real/unique points
        do iz = 2, nz-1
            do ix = 2, nx-1
                ! leapfrog scheme
                if (it == 1) then
                    up(ix,iz,3) = up(ix,iz,1) - (cx*dt*rd2x*(up(ix+1,iz,2)-up(ix-1,iz,2))) - (cz*dt*rd2z*(up(ix,iz+1,2)-up(ix,iz-1,2))) 
                else
                    up(ix,iz,3) = up(ix,iz,1) - (cx*d2t*rd2x*(up(ix+1,iz,2)-up(ix-1,iz,2))) - (cz*d2t*rd2z*(up(ix,iz+1,2)-up(ix,iz-1,2))) 
                endif
            enddo
        enddo


        ! take care of boundaries
        if (pbc==.True.) then 
            do iz = 2,nz-1
                pip(1,iz,3) = pip(nx-1,iz,3)
                pip(nx,iz,3) = pip(2,iz,3)
                thp(1,iz,3) = thp(nx-1,iz,3)
                thp(nx,iz,3) = thp(2,iz,3)
                up(1,iz,3) = up(nx-1,iz,3)
                up(nx,iz,3) = up(2,iz,3)
                wp(1,iz,3) = wp(nx-1,iz,3)
                wp(nx,iz,3) = wp(2,iz,3)
            enddo

            do ix = 1,nx
                pip(ix,1,3) = pip(ix,nz-1,3)
                pip(ix,nz,3) = pip(ix, 2,3)
                thp(ix,1,3) = thp(ix,nz-1,3)
                thp(ix,nz,3) = thp(ix,2,3)
                up(ix,1,3) = up(ix,nz-1,3)
                up(ix,nz,3) = up(ix,2,3)
                wp(ix,1,3) = wp(ix,nz-1,3)
                wp(ix,nz,3) = wp(ix,2,3)
            enddo
        else 
            ! if not using PBCs, just set 1st point to be same as 2nd, last pt to be same as 2nd to last, etc.
            do iz = 2,nz-1
                pip(1,iz,3) = pip(2,iz,3)
                pip(nx,iz,3) = pip(nx-1,iz,3)
                thp(1,iz,3) = thp(2,iz,3)
                thp(nx,iz,3) = thp(nx-1,iz,3)
                up(1,iz,3) = up(2,iz,3)
                up(nx,iz,3) = up(nx-1,iz,3)
                wp(1,iz,3) = wp(2,iz,3)
                wp(nx,iz,3) = wp(nx-1,iz,3)
            enddo

            do ix = 1,nx
                pip(ix,1,3) = pip(ix,2,3)
                pip(ix,nz,3) = pip(ix, nz-1,3)
                thp(ix,1,3) = thp(ix,2,3)
                thp(ix,nz,3) = thp(ix,nz-1,3)
                up(ix,1,3) = up(ix,2,3)
                up(ix,nz,3) = up(ix,nz-1,3)
                wp(ix,1,3) = wp(ix,2,3)
                wp(ix,nz,3) = wp(ix,nz-1,3)
            enddo
        endif 

        ! step forward in time
        do iz = 1, nz
            do ix = 1, nx
                up(ix,iz,1) = up(ix,iz,2)
                up(ix,iz,2) = up(ix,iz,3)
            enddo
        enddo

    end subroutine advect

end module advection
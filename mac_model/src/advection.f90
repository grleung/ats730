module advection

    ! This module contains the subroutine to advect a perturbation

    implicit none

    contains

    subroutine advect
        use run_constants, only: nz,ny,nx,dx,dy,dz0,cx,cy,cz,dt,pbc_x,pbc_y,pbc_z
        use model_vars, only:it,zsn,xsn,ysn,thb,thvb,rhoub,thp,pip,pp,up,vp,wp

        implicit none

        integer :: iz ! counter for z-coordinate
        integer :: iy ! counter for y-coordinate
        integer :: ix ! counter for x-coordinate

        real :: rd2x,rd2y,rd2z,d2t
        
        rd2x      = (1/(dx+dx))  ! reciprocal of 2dx [1/m]
        rd2y      = (1/(dy+dy))  ! reciprocal of 2dx [1/m]
        rd2z      = (1/(dz0+dz0))  ! reciprocal of 2dz [1/m]
        d2t       = (dt+dt)         ! 2*dt

        ! loop over real/unique points
        do iz = 2, nz-1
            do iy = 2, ny-1
                do ix = 2, nx-1
                    ! leapfrog scheme
                    if (it == 1) then
                        up(ix,iy,iz,3) = up(ix,iy,iz,1) - (cx*dt*rd2x*(up(ix+1,iy,iz,2)-up(ix-1,iy,iz,2))) - (cz*dt*rd2z*(up(ix,iy,iz+1,2)-up(ix,iy,iz-1,2))) 
                    else
                        up(ix,iy,iz,3) = up(ix,iy,iz,1) - (cx*d2t*rd2x*(up(ix+1,iy,iz,2)-up(ix-1,iy,iz,2))) - (cz*d2t*rd2z*(up(ix,iy,iz+1,2)-up(ix,iy,iz-1,2))) 
                    endif
                enddo
            enddo
        enddo

        ! step forward in time
        do iz = 1, nz
            do iy=1,ny
                do ix = 1, nx
                    up(ix,iy,iz,1) = up(ix,iy,iz,2)
                    up(ix,iy,iz,2) = up(ix,iy,iz,3)
                enddo
            enddo
        enddo

    end subroutine advect

end module advection
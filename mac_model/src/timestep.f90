module timestep
    ! This module contains the subroutine to step model vars forward in time

    implicit none

    contains
    
    subroutine step_time
        use model_vars, only: thp,pip,pp,up,wp,rvp,np,mp,nd,mld,mpd
        use run_constants, only: nz,nx,npartbin,ndropbin

        implicit none

        integer :: iz ! counter for z
        integer :: ix ! counter for x
        integer :: ipartbin,idropbin,it


        print*,'step time start'
        ! step forward in time
        do ix = 1, nx
            do iz = 1, nz
                !do 1 = 1,2
                up(ix,iz,1) = up(ix,iz,2)
                wp(ix,iz,1) = wp(ix,iz,2)
                thp(ix,iz,1) = thp(ix,iz,2)
                pip(ix,iz,1) = pip(ix,iz,2)
                rvp(ix,iz,1) = rvp(ix,iz,2)

                do ipartbin=1,npartbin
                    np(ix,iz,ipartbin,1) = np(ix,iz,ipartbin,2)
                    mp(ix,iz,ipartbin,1) = mp(ix,iz,ipartbin,2)
                    
                    do idropbin=1,ndropbin
                        mld(ix,iz,ipartbin,idropbin,1) = mld(ix,iz,ipartbin,idropbin,2)
                        mpd(ix,iz,ipartbin,idropbin,1) = mpd(ix,iz,ipartbin,idropbin,2)
                        nd(ix,iz,ipartbin,idropbin,1) = nd(ix,iz,ipartbin,idropbin,2)
                    enddo
                enddo

                up(ix,iz,2) = up(ix,iz,3)
                wp(ix,iz,2) = wp(ix,iz,3)
                thp(ix,iz,2) = thp(ix,iz,3)
                pip(ix,iz,2) = pip(ix,iz,3)
                rvp(ix,iz,2) = rvp(ix,iz,3)

                do ipartbin=1,npartbin
                    np(ix,iz,ipartbin,2) = np(ix,iz,ipartbin,3)
                    mp(ix,iz,ipartbin,2) = mp(ix,iz,ipartbin,3)
                    
                    do idropbin=1,ndropbin
                        mld(ix,iz,ipartbin,idropbin,2) = mld(ix,iz,ipartbin,idropbin,3)
                        mpd(ix,iz,ipartbin,idropbin,2) = mpd(ix,iz,ipartbin,idropbin,3)
                        
                        if (nd(ix,iz,ipartbin,idropbin,3)<0.) then
                            print*,''
                            nd(ix,iz,ipartbin,idropbin,2) = 0.
                        else
                            nd(ix,iz,ipartbin,idropbin,2) = nd(ix,iz,ipartbin,idropbin,3)
                        endif 
                    enddo
                enddo
                !enddo
            enddo ! end z loop
        enddo ! end x loop

        print*,MAXVAL(wp(:,:,2))
        print*,MAXVAL(rvp(:,:,2))
        !print*,SUM(np(:,:,:,2))+SUM(nc(:,:,:,:,2))
        
    end subroutine step_time

end module timestep
module mem
    use run_constants, only: nz,nx,ny
    use model_vars
    
    implicit none

    contains

    subroutine allocate_mem
        ! grid arrays
        allocate (dzn(nz))
        allocate (zsn(nz))
        allocate (zmn(nz))
        allocate (dxn(nx))
        allocate (xsn(nx))
        allocate (xmn(nx))
        allocate (dyn(ny))
        allocate (ysn(ny))
        allocate (ymn(ny))

        ! base state arrays
        allocate (tb(nz))
        allocate (thb(nz))
        allocate (rvb(nz))
        allocate (thvb(nz))
        allocate (pb(nz))
        allocate (pib(nz))
        allocate (piwb(nz))
        allocate (rhoub(nz))
        allocate (rhowb(nz))
        allocate (satfracb(nz))
        allocate (rhb(nz))
        allocate (rsatb(nz))

        ! parcel arrays
        allocate (tpcl(nz))
        allocate (thpcl(nz))
        allocate (rvpcl(nz))
        allocate (thvpcl(nz))
        allocate (thvdiff(nz))
        allocate (rsatpcl(nz))

        ! prognostic variable arrays
        allocate (thp(nx,ny,nz,3))
        allocate (rvp(nx,ny,nz,3))
        allocate (pip(nx,ny,nz,3))
        allocate (up(nx,ny,nz,3))
        allocate (vp(nx,ny,nz,3))
        allocate (wp(nx,ny,nz,3))
        allocate (pp(nx,ny,nz))

        ! tendency arrays
        allocate (u_xadv(nx,ny,nz))
        allocate (u_yadv(nx,ny,nz))
        allocate (u_zadv(nx,ny,nz))
        allocate (u_pgf(nx,ny,nz))
        allocate (u_xdiff(nx,ny,nz))
        allocate (u_ydiff(nx,ny,nz))
        allocate (u_zdiff(nx,ny,nz))
        allocate (u_tend_total(nx,ny,nz))
        allocate (v_xadv(nx,ny,nz))
        allocate (v_yadv(nx,ny,nz))
        allocate (v_zadv(nx,ny,nz))
        allocate (v_pgf(nx,ny,nz))
        allocate (v_xdiff(nx,ny,nz))
        allocate (v_ydiff(nx,ny,nz))
        allocate (v_zdiff(nx,ny,nz))
        allocate (v_tend_total(nx,ny,nz))
        allocate (w_xadv(nx,ny,nz))
        allocate (w_yadv(nx,ny,nz))
        allocate (w_zadv(nx,ny,nz))
        allocate (w_pgf(nx,ny,nz))
        allocate (w_buoy(nx,ny,nz))
        allocate (w_xdiff(nx,ny,nz))
        allocate (w_ydiff(nx,ny,nz))
        allocate (w_zdiff(nx,ny,nz))
        allocate (w_tend_total(nx,ny,nz))
        allocate (thp_xadv(nx,ny,nz))
        allocate (thp_yadv(nx,ny,nz))
        allocate (thp_zadv(nx,ny,nz))
        allocate (thp_meanadv(nx,ny,nz))
        allocate (thp_xdiff(nx,ny,nz))
        allocate (thp_ydiff(nx,ny,nz))
        allocate (thp_zdiff(nx,ny,nz))
        allocate (thp_tend_total(nx,ny,nz))
        allocate (pip_xadv(nx,ny,nz))
        allocate (pip_yadv(nx,ny,nz))
        allocate (pip_zadv(nx,ny,nz))
        allocate (pip_xdiff(nx,ny,nz))
        allocate (pip_ydiff(nx,ny,nz))
        allocate (pip_zdiff(nx,ny,nz))
        allocate (pip_tend_total(nx,ny,nz))

    end subroutine allocate_mem

    subroutine deallocate_mem
        ! grid arrays
        deallocate (dzn(nz))
        deallocate (zsn(nz))
        deallocate (zmn(nz))
        deallocate (dxn(nx))
        deallocate (xsn(nx))
        deallocate (xmn(nx))
        deallocate (dyn(ny))
        deallocate (ysn(ny))
        deallocate (ymn(ny))
        

        ! base state arrays
        deallocate (tb(nz))
        deallocate (thb(nz))
        deallocate (rvb(nz))
        deallocate (thvb(nz))
        deallocate (pb(nz))
        deallocate (pib(nz))
        deallocate (piwb(nz))
        deallocate (rhoub(nz))
        deallocate (rhowb(nz))
        deallocate (satfracb(nz))
        deallocate (rhb(nz))
        deallocate (rsatb(nz))

        ! parcel arrays
        deallocate (tpcl(nz))
        deallocate (thpcl(nz))
        deallocate (rvpcl(nz))
        deallocate (thvpcl(nz))
        deallocate (thvdiff(nz))
        deallocate (rsatpcl(nz))

        ! prognostic variable arrays
        deallocate (thp(nx,ny,nz,3))
        deallocate (rvp(nx,ny,nz,3))
        deallocate (pip(nx,ny,nz,3))
        deallocate (up(nx,ny,nz,3))
        deallocate (vp(nx,ny,nz,3))
        deallocate (wp(nx,ny,nz,3))
        deallocate (pp(nx,ny,nz))

        ! tendency arrays
        deallocate (u_xadv(nx,ny,nz))
        deallocate (u_yadv(nx,ny,nz))
        deallocate (u_zadv(nx,ny,nz))
        deallocate (u_pgf(nx,ny,nz))
        deallocate (u_xdiff(nx,ny,nz))
        deallocate (u_ydiff(nx,ny,nz))
        deallocate (u_zdiff(nx,ny,nz))
        deallocate (u_tend_total(nx,ny,nz))
        deallocate (v_xadv(nx,ny,nz))
        deallocate (v_yadv(nx,ny,nz))
        deallocate (v_zadv(nx,ny,nz))
        deallocate (v_pgf(nx,ny,nz))
        deallocate (v_xdiff(nx,ny,nz))
        deallocate (v_ydiff(nx,ny,nz))
        deallocate (v_zdiff(nx,ny,nz))
        deallocate (v_tend_total(nx,ny,nz))
        deallocate (w_xadv(nx,ny,nz))
        deallocate (w_yadv(nx,ny,nz))
        deallocate (w_zadv(nx,ny,nz))
        deallocate (w_pgf(nx,ny,nz))
        deallocate (w_buoy(nx,ny,nz))
        deallocate (w_xdiff(nx,ny,nz))
        deallocate (w_ydiff(nx,ny,nz))
        deallocate (w_zdiff(nx,ny,nz))
        deallocate (w_tend_total(nx,ny,nz))
        deallocate (thp_xadv(nx,ny,nz))
        deallocate (thp_yadv(nx,ny,nz))
        deallocate (thp_zadv(nx,ny,nz))
        deallocate (thp_meanadv(nx,ny,nz))
        deallocate (thp_xdiff(nx,ny,nz))
        deallocate (thp_ydiff(nx,ny,nz))
        deallocate (thp_zdiff(nx,ny,nz))
        deallocate (thp_tend_total(nx,ny,nz))
        deallocate (pip_xadv(nx,ny,nz))
        deallocate (pip_yadv(nx,ny,nz))
        deallocate (pip_zadv(nx,ny,nz))
        deallocate (pip_xdiff(nx,ny,nz))
        deallocate (pip_ydiff(nx,ny,nz))
        deallocate (pip_zdiff(nx,ny,nz))
        deallocate (pip_tend_total(nx,ny,nz))
        
    end subroutine deallocate_mem


end module mem
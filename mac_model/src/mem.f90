module mem
    use run_constants, only: nz,nx,npartbin,ndropbin,bin_flag
    use model_vars
    
    implicit none

    contains

    subroutine allocate_mem
        ! grid arrays
        allocate (dzn(nz))
        allocate (zsn(nz))
        allocate (zwn(nz))
        allocate (dxn(nx))
        allocate (xsn(nx))
        allocate (xun(nx))

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
        allocate (thp(nx,nz,3))
        allocate (rvp(nx,nz,3))
        allocate (pip(nx,nz,3))
        allocate (up(nx,nz,3))
        allocate (wp(nx,nz,3))

        allocate (pp(nx,nz))
        allocate (thvp(nx,nz))
        allocate (p(nx,nz))
        allocate (t(nx,nz))
        allocate (samb(nx,nz))

        allocate(mdropbin_lims(ndropbin+1))
        allocate(mpartbin_lims(npartbin+1))
        allocate(np(nx,nz,npartbin,3))
        allocate(mp(nx,nz,npartbin,3))
        allocate(nd(nx,nz,npartbin,ndropbin,3))
        allocate(mld(nx,nz,npartbin,ndropbin,3))
        allocate(mpd(nx,nz,npartbin,ndropbin,3))

        ! tendency arrays
        allocate (u_xadv(nx,nz))
        allocate (u_zadv(nx,nz))
        allocate (u_pgf(nx,nz))
        allocate (u_xdiff(nx,nz))
        allocate (u_zdiff(nx,nz))
        allocate (u_tend_total(nx,nz))
        allocate (w_xadv(nx,nz))
        allocate (w_zadv(nx,nz))
        allocate (w_pgf(nx,nz))
        allocate (w_buoy(nx,nz))
        allocate (w_xdiff(nx,nz))
        allocate (w_zdiff(nx,nz))
        allocate (w_tend_total(nx,nz))
        allocate (thp_xadv(nx,nz))
        allocate (thp_zadv(nx,nz))
        allocate (thp_meanadv(nx,nz))
        allocate (thp_xdiff(nx,nz))
        allocate (thp_zdiff(nx,nz))
        allocate (thp_tend_total(nx,nz))
        allocate (pip_xadv(nx,nz))
        allocate (pip_zadv(nx,nz))
        allocate (pip_xdiff(nx,nz))
        allocate (pip_zdiff(nx,nz))
        allocate (pip_tend_total(nx,nz))

        allocate(np_tend_total(nx,nz,npartbin))
        allocate(mp_tend_total(nx,nz,npartbin))
        allocate(nd_tend_total(nx,nz,npartbin,ndropbin))
        allocate(mld_tend_total(nx,nz,npartbin,ndropbin))
        allocate(mpd_tend_total(nx,nz,npartbin,ndropbin))
        allocate (rvp_xadv(nx,nz))
        allocate (rvp_zadv(nx,nz))
        allocate (rvp_meanadv(nx,nz))
        allocate (rvp_xdiff(nx,nz))
        allocate (rvp_zdiff(nx,nz))
        allocate (rvp_tend_total(nx,nz))
        
    end subroutine allocate_mem

    subroutine deallocate_mem
        ! grid arrays
        deallocate (dzn(nz))
        deallocate (zsn(nz))
        deallocate (zwn(nz))
        deallocate (dxn(nx))
        deallocate (xsn(nx))
        deallocate (xun(nx))

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
        deallocate (thp(nx,nz,3))
        deallocate (rvp(nx,nz,3))
        deallocate (pip(nx,nz,3))
        deallocate (up(nx,nz,3))
        deallocate (wp(nx,nz,3))
        deallocate (pp(nx,nz))
        deallocate (thvp(nx,nz))
        deallocate (p(nx,nz))
        deallocate (t(nx,nz))
        deallocate (samb(nx,nz))

        ! bin microphysics arrays
        deallocate(mdropbin_lims(ndropbin+1))
        deallocate(mpartbin_lims(npartbin+1))
        deallocate(np(nx,nz,npartbin,3))
        deallocate(mp(nx,nz,npartbin,3))
        deallocate(nd(nx,nz,npartbin,ndropbin,3))
        deallocate(mld(nx,nz,npartbin,ndropbin,3))
        deallocate(mpd(nx,nz,npartbin,ndropbin,3))

        ! tendency arrays
        deallocate (u_xadv(nx,nz))
        deallocate (u_zadv(nx,nz))
        deallocate (u_pgf(nx,nz))
        deallocate (u_xdiff(nx,nz))
        deallocate (u_zdiff(nx,nz))
        deallocate (u_tend_total(nx,nz))
        deallocate (w_xadv(nx,nz))
        deallocate (w_zadv(nx,nz))
        deallocate (w_pgf(nx,nz))
        deallocate (w_buoy(nx,nz))
        deallocate (w_xdiff(nx,nz))
        deallocate (w_zdiff(nx,nz))
        deallocate (w_tend_total(nx,nz))
        deallocate (thp_xadv(nx,nz))
        deallocate (thp_zadv(nx,nz))
        deallocate (thp_meanadv(nx,nz))
        deallocate (thp_xdiff(nx,nz))
        deallocate (thp_zdiff(nx,nz))
        deallocate (thp_tend_total(nx,nz))
        deallocate (pip_xadv(nx,nz))
        deallocate (pip_zadv(nx,nz))
        deallocate (pip_xdiff(nx,nz))
        deallocate (pip_zdiff(nx,nz))
        deallocate (pip_tend_total(nx,nz))

        deallocate(np_tend_total(nx,nz,npartbin))
        deallocate(mp_tend_total(nx,nz,npartbin))
        deallocate(nd_tend_total(nx,nz,npartbin,ndropbin))
        deallocate(mld_tend_total(nx,nz,npartbin,ndropbin))
        deallocate(mpd_tend_total(nx,nz,npartbin,ndropbin))
        deallocate (rvp_xadv(nx,nz))
        deallocate (rvp_zadv(nx,nz))
        deallocate (rvp_meanadv(nx,nz))
        deallocate (rvp_xdiff(nx,nz))
        deallocate (rvp_zdiff(nx,nz))
        deallocate (rvp_tend_total(nx,nz))
    
    end subroutine deallocate_mem


end module mem
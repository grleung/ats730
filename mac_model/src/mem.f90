module mem
    use run_constants, only: nz,nx
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
        allocate (rcp(nx,nz,3))
        allocate (rrp(nx,nz,3))

        allocate (pp(nx,nz))
        allocate (thvp(nx,nz))
        allocate (vap2cld(nx,nz))
        allocate (rain2vap(nx,nz))
        allocate (cld2rain_accr(nx,nz))
        allocate (cld2rain_auto(nx,nz))

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

        allocate (rvp_xadv(nx,nz))
        allocate (rvp_zadv(nx,nz))
        allocate (rvp_meanadv(nx,nz))
        allocate (rvp_xdiff(nx,nz))
        allocate (rvp_zdiff(nx,nz))
        allocate (rvp_tend_total(nx,nz))
        allocate (rcp_xadv(nx,nz))
        allocate (rcp_zadv(nx,nz))
        allocate (rcp_xdiff(nx,nz))
        allocate (rcp_zdiff(nx,nz))
        allocate (rcp_tend_total(nx,nz))
        allocate (rrp_xadv(nx,nz))
        allocate (rrp_zadv(nx,nz))
        allocate (rrp_xdiff(nx,nz))
        allocate (rrp_zdiff(nx,nz))
        allocate (rrp_tend_total(nx,nz))

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
        deallocate (rcp(nx,nz,3))
        deallocate (rrp(nx,nz,3))
        deallocate (pp(nx,nz))
        deallocate (thvp(nx,nz))
        deallocate (vap2cld(nx,nz))
        deallocate (rain2vap(nx,nz))
        deallocate (cld2rain_accr(nx,nz))
        deallocate (cld2rain_auto(nx,nz))

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
        
        deallocate (rvp_xadv(nx,nz))
        deallocate (rvp_zadv(nx,nz))
        deallocate (rvp_meanadv(nx,nz))
        deallocate (rvp_xdiff(nx,nz))
        deallocate (rvp_zdiff(nx,nz))
        deallocate (rvp_tend_total(nx,nz))
        deallocate (rcp_xadv(nx,nz))
        deallocate (rcp_zadv(nx,nz))
        deallocate (rcp_xdiff(nx,nz))
        deallocate (rcp_zdiff(nx,nz))
        deallocate (rcp_tend_total(nx,nz))
        deallocate (rrp_xadv(nx,nz))
        deallocate (rrp_zadv(nx,nz))
        deallocate (rrp_xdiff(nx,nz))
        deallocate (rrp_zdiff(nx,nz))
        deallocate (rrp_tend_total(nx,nz))

    end subroutine deallocate_mem


end module mem
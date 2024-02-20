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
        allocate (pp(nx,nz))

        ! tendency arrays
        allocate (u_tend1(nx,nz))
        allocate (u_tend2(nx,nz))
        allocate (u_tend3(nx,nz))
        allocate (u_tend_total(nx,nz))
        allocate (w_tend1(nx,nz))
        allocate (w_tend2(nx,nz))
        allocate (w_tend3(nx,nz))
        allocate (w_tend4(nx,nz))
        allocate (w_tend_total(nx,nz))
        allocate (thp_tend1(nx,nz))
        allocate (thp_tend2(nx,nz))
        allocate (thp_tend3(nx,nz))
        allocate (thp_tend_total(nx,nz))
        allocate (pip_tend1(nx,nz))
        allocate (pip_tend2(nx,nz))
        allocate (pip_tend_total(nx,nz))

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

        ! tendency arrays
        deallocate (u_tend1(nx,nz))
        deallocate (u_tend2(nx,nz))
        deallocate (u_tend3(nx,nz))
        deallocate (u_tend_total(nx,nz))
        deallocate (w_tend1(nx,nz))
        deallocate (w_tend2(nx,nz))
        deallocate (w_tend3(nx,nz))
        deallocate (w_tend4(nx,nz))
        deallocate (w_tend_total(nx,nz))
        deallocate (thp_tend1(nx,nz))
        deallocate (thp_tend2(nx,nz))
        deallocate (thp_tend3(nx,nz))
        deallocate (thp_tend_total(nx,nz))
        deallocate (pip_tend1(nx,nz))
        deallocate (pip_tend2(nx,nz))
        deallocate (pip_tend_total(nx,nz))
        
    end subroutine deallocate_mem


end module mem
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
        allocate (tb(nx,nz))
        allocate (thb(nx,nz))
        allocate (rvb(nx,nz))
        allocate (thvb(nx,nz))
        allocate (pb(nx,nz))
        allocate (pib(nx,nz))
        allocate (piwb(nx,nz))
        allocate (rhoub(nx,nz))
        allocate (rhowb(nx,nz))
        allocate (satfracb(nx,nz))
        allocate (rhb(nx,nz))
        allocate (rsatb(nx,nz))

        ! parcel arrays
        allocate (tp(nz))
        allocate (thp(nz))
        allocate (rvp(nz))
        allocate (thvp(nz))
        allocate (thvdiff(nz))
        allocate (rsatp(nz))
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
        deallocate (tb(nx,nz))
        deallocate (thb(nx,nz))
        deallocate (rvb(nx,nz))
        deallocate (thvb(nx,nz))
        deallocate (pb(nx,nz))
        deallocate (pib(nx,nz))
        deallocate (piwb(nx,nz))
        deallocate (rhoub(nx,nz))
        deallocate (rhowb(nx,nz))
        deallocate (satfracb(nx,nz))
        deallocate (rhb(nx,nz))
        deallocate (rsatb(nx,nz))

        ! parcel arrays
        deallocate (tp(nz))
        deallocate (thp(nz))
        deallocate (rvp(nz))
        deallocate (thvp(nz))
        deallocate (thvdiff(nz))
        deallocate (rsatp(nz))
    end subroutine deallocate_mem


end module mem
module mem
    use grid_constants, only: nz
    use model_vars
    
    implicit none

    contains

    subroutine allocate_mem
        allocate (dzn(nz))
        allocate (zun(nz))
        allocate (zwn(nz))
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
        allocate (tp(nz))
        allocate (thp(nz))
        allocate (rvp(nz))
        allocate (thvp(nz))
        allocate (thvdiff(nz))
        allocate (rsatp(nz))
    end subroutine allocate_mem

    subroutine deallocate_mem
        deallocate (dzn(nz))
        deallocate (zun(nz))
        deallocate (zwn(nz))
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
        deallocate (tp(nz))
        deallocate (thp(nz))
        deallocate (rvp(nz))
        deallocate (thvp(nz))
        deallocate (thvdiff(nz))
        deallocate (rsatp(nz))
    end subroutine deallocate_mem


end module mem
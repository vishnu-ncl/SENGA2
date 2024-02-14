SUBROUTINE ops_data_init()
    use OPS_Fortran_Reference

    use com_senga
    use com_ops_senga

    integer(kind=4) :: d_size(3)
    integer(kind=4) :: d_base(3) = [1,1,1] !array indexing - start from 1
    integer(kind=4) :: d_p(3) !max boundary depths for the dat in the possitive direction
    integer(kind=4) :: d_m(3) !max boundary depths for the dat in the negative direction
    integer(kind=4) :: ispec
    character(len=20) :: buf

!   *-----------------------------------------OPS Declarations-----------------------------------------------*

!   Declare OPS Block
    call ops_decl_block(3, senga_grid, "SENGA_GRID")

!   Declare OPS Dats
    d_size = [nxglbl, nyglbl, nzglbl]
    d_m    = [0,0,0]
    d_p    = [0,0,0]

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store1, d_store1, "real(kind=8)", "STORE1")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store2, d_store2, "real(kind=8)", "STORE2")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store3, d_store3, "real(kind=8)", "STORE3")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store4, d_store4, "real(kind=8)", "STORE4")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store5, d_store5, "real(kind=8)", "STORE5")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store6, d_store6, "real(kind=8)", "STORE6")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, divm, d_divm, "real(kind=8)", "DIVM")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, drun, d_drun, "real(kind=8)", "DRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, urun, d_urun, "real(kind=8)", "URUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vrun, d_vrun, "real(kind=8)", "VRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wrun, d_wrun, "real(kind=8)", "WRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, erun, d_erun, "real(kind=8)", "ERUN")

!---------------------------------------WITH HALOS-----------------------------------------------------------

    d_size = [nxglbl, nyglbl, nzglbl]
    d_m    = [-nhalox,-nhaloy,-nhaloz]
    d_p    = [nhalox,nhaloy,nhaloz]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, drhs, d_drhs, "real(kind=8)", "DRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, urhs, d_urhs, "real(kind=8)", "URHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vrhs, d_vrhs, "real(kind=8)", "VRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wrhs, d_wrhs, "real(kind=8)", "WRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, erhs, d_erhs, "real(kind=8)", "ERHS")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, utmp, d_utmp, "real(kind=8)", "UTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vtmp, d_vtmp, "real(kind=8)", "VTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wtmp, d_wtmp, "real(kind=8)", "WTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, prun, d_prun, "real(kind=8)", "PRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, trun, d_trun, "real(kind=8)", "TRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, transp, d_transp, "real(kind=8)", "TRANSP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store7, d_store7, "real(kind=8)", "STORE7")

!---------------------------------------MULTI-DIM DAT--------------------------------------------------------

    d_size = [nxglbl, nyglbl, nzglbl]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    DO ispec = 1,nspcmx
        WRITE(buf,"(A4,I2.2)") "YRUN",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun(:,:,:,ispec), d_yrun(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A4,I2.2)") "RATE",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, rate(:,:,:,ispec), d_rate(ispec), "real(kind=8)", trim(buf))
        WRITE(buf,"(A4,I2.2)") "RRTE",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, rrte(:,:,:,ispec), d_rrte(ispec), "real(kind=8)", trim(buf))
    END DO

    d_size = [ nxglbl, nyglbl, nzglbl]
    d_m    = [-nhalox,-nhaloy,-nhaloz]
    d_p    = [ nhalox, nhaloy, nhaloz]
    DO ispec = 1,nspcmx
        WRITE(buf,"(A4,I2.2)") "YRHS",ispec
        call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrhs(:,:,:,ispec), d_yrhs(ispec), "real(kind=8)", trim(buf))
    END DO
!------------------------------------------------------------------------------------------------------------
    call ops_partition(" ")
!------------------------------------------------------------------------------------------------------------

END SUBROUTINE ops_data_init

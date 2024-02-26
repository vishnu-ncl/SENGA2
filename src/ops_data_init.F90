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

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ucor, d_ucor, "real(kind=8)", "UCOR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vcor, d_vcor, "real(kind=8)", "VCOR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wcor, d_wcor, "real(kind=8)", "WCOR")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd1x, d_wd1x, "real(kind=8)", "WD1X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd1x, d_pd1x, "real(kind=8)", "PD1X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td1x, d_td1x, "real(kind=8)", "TD1X")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd1y, d_wd1y, "real(kind=8)", "WD1Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd1y, d_pd1y, "real(kind=8)", "PD1Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td1y, d_td1y, "real(kind=8)", "TD1Y")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd1z, d_wd1z, "real(kind=8)", "WD1Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd1z, d_pd1z, "real(kind=8)", "PD1Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td1z, d_td1z, "real(kind=8)", "TD1Z")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd2x, d_wd2x, "real(kind=8)", "WD2X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd2x, d_pd2x, "real(kind=8)", "PD2X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td2x, d_td2x, "real(kind=8)", "TD2X")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd2y, d_wd2y, "real(kind=8)", "WD2Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd2y, d_pd2y, "real(kind=8)", "PD2Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td2y, d_td2y, "real(kind=8)", "TD2Y")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd2z, d_wd2z, "real(kind=8)", "WD2Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd2z, d_pd2z, "real(kind=8)", "PD2Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td2z, d_td2z, "real(kind=8)", "TD2Z")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ufxl, d_ufxl, "real(kind=8)", "UFXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vfxl, d_vfxl, "real(kind=8)", "VFXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wfxl, d_wfxl, "real(kind=8)", "WFXL")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, drun, d_drun, "real(kind=8)", "DRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, urun, d_urun, "real(kind=8)", "URUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vrun, d_vrun, "real(kind=8)", "VRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wrun, d_wrun, "real(kind=8)", "WRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, erun, d_erun, "real(kind=8)", "ERUN")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, derr, d_derr, "real(kind=8)", "DERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, uerr, d_uerr, "real(kind=8)", "UERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, verr, d_verr, "real(kind=8)", "VERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, werr, d_werr, "real(kind=8)", "WERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, eerr, d_eerr, "real(kind=8)", "EERR")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, prun2, d2prun, "real(kind=8)", "PRN2")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, trun2, d2trun, "real(kind=8)", "TRN2")

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

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wmomix, d_wmomix, "real(kind=8)", "WMOMIX")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, difmix, d_difmix, "real(kind=8)", "DIFMIX")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, tdrmix, d_tdrmix, "real(kind=8)", "TDRMIX")

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

!-----------------------------------------Boundary YZ--------------------------------------------------------

    d_size = [1,nyglbl,nzglbl]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1xl, d_bcl1xl, "real(kind=8)", "BCL1XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1xr, d_bcl1xr, "real(kind=8)", "BCL1XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2xl, d_bcl2xl, "real(kind=8)", "BCL2XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2xr, d_bcl2xr, "real(kind=8)", "BCL2XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3xl, d_bcl3xl, "real(kind=8)", "BCL3XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3xr, d_bcl3xr, "real(kind=8)", "BCL3XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4xl, d_bcl4xl, "real(kind=8)", "BCL4XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4xr, d_bcl4xr, "real(kind=8)", "BCL4XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5xl, d_bcl5xl, "real(kind=8)", "BCL5XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5xr, d_bcl5xr, "real(kind=8)", "BCL5XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltxl, d_bcltxl, "real(kind=8)", "BCLTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltxr, d_bcltxr, "real(kind=8)", "BCLTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struxl, d_struxl, "real(kind=8)", "STRUXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struxr, d_struxr, "real(kind=8)", "STRUXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvxl, d_strvxl, "real(kind=8)", "STRVXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvxr, d_strvxr, "real(kind=8)", "STRVXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwxl, d_strwxl, "real(kind=8)", "STRWXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwxr, d_strwxr, "real(kind=8)", "STRWXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpxl, d_strpxl, "real(kind=8)", "STRPXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpxr, d_strpxr, "real(kind=8)", "STRPXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdxl, d_strdxl, "real(kind=8)", "STRDXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdxr, d_strdxr, "real(kind=8)", "STRDXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtxl, d_strtxl, "real(kind=8)", "STRTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtxr, d_strtxr, "real(kind=8)", "STRTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strexl, d_strexl, "real(kind=8)", "STREXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strexr, d_strexr, "real(kind=8)", "STREXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgxl, d_strgxl, "real(kind=8)", "STRGXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgxr, d_strgxr, "real(kind=8)", "STRGXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrxl, d_strrxl, "real(kind=8)", "STRRXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrxr, d_strrxr, "real(kind=8)", "STRRXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtxl, d_dudtxl, "real(kind=8)", "DUDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtxr, d_dudtxr, "real(kind=8)", "DUDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtxl, d_dvdtxl, "real(kind=8)", "DVDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtxr, d_dvdtxr, "real(kind=8)", "DVDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtxl, d_dwdtxl, "real(kind=8)", "DWDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtxr, d_dwdtxr, "real(kind=8)", "DWDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtxl, d_dtdtxl, "real(kind=8)", "DTDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtxr, d_dtdtxr, "real(kind=8)", "DTDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtxl, d_dddtxl, "real(kind=8)", "DDDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtxr, d_dddtxr, "real(kind=8)", "DDDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouxl, d_acouxl, "real(kind=8)", "ACOUXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouxr, d_acouxr, "real(kind=8)", "ACOUXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2xl, d_ova2xl, "real(kind=8)", "OVA2XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2xr, d_ova2xr, "real(kind=8)", "OVA2XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1xl, d_gam1xl, "real(kind=8)", "GAM1XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1xr, d_gam1xr, "real(kind=8)", "GAM1XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmxl, d_ovgmxl, "real(kind=8)", "OVGMXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmxr, d_ovgmxr, "real(kind=8)", "OVGMXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtxl, d_sydtxl, "real(kind=8)", "SYDTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtxr, d_sydtxr, "real(kind=8)", "SYDTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpxl, d_sorpxl, "real(kind=8)", "SORPXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpxr, d_sorpxr, "real(kind=8)", "SORPXR")

!-----------------------------------------Boundary XZ--------------------------------------------------------

    d_size = [nxglbl,1,nzglbl]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1yl, d_bcl1yl, "real(kind=8)", "BCL1YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1yr, d_bcl1yr, "real(kind=8)", "BCL1YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2yl, d_bcl2yl, "real(kind=8)", "BCL2YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2yr, d_bcl2yr, "real(kind=8)", "BCL2YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3yl, d_bcl3yl, "real(kind=8)", "BCL3YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3yr, d_bcl3yr, "real(kind=8)", "BCL3YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4yl, d_bcl4yl, "real(kind=8)", "BCL4YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4yr, d_bcl4yr, "real(kind=8)", "BCL4YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5yl, d_bcl5yl, "real(kind=8)", "BCL5YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5yr, d_bcl5yr, "real(kind=8)", "BCL5YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltyl, d_bcltyl, "real(kind=8)", "BCLTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltyr, d_bcltyr, "real(kind=8)", "BCLTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struyl, d_struyl, "real(kind=8)", "STRUYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struyr, d_struyr, "real(kind=8)", "STRUYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvyl, d_strvyl, "real(kind=8)", "STRVYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvyr, d_strvyr, "real(kind=8)", "STRVYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwyl, d_strwyl, "real(kind=8)", "STRWYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwyr, d_strwyr, "real(kind=8)", "STRWYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpyl, d_strpyl, "real(kind=8)", "STRPYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpyr, d_strpyr, "real(kind=8)", "STRPYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdyl, d_strdyl, "real(kind=8)", "STRDYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdyr, d_strdyr, "real(kind=8)", "STRDYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtyl, d_strtyl, "real(kind=8)", "STRTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtyr, d_strtyr, "real(kind=8)", "STRTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, streyl, d_streyl, "real(kind=8)", "STREYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, streyr, d_streyr, "real(kind=8)", "STREYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgyl, d_strgyl, "real(kind=8)", "STRGYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgyr, d_strgyr, "real(kind=8)", "STRGYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strryl, d_strryl, "real(kind=8)", "STRRYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strryr, d_strryr, "real(kind=8)", "STRRYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtyl, d_dudtyl, "real(kind=8)", "DUDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtyr, d_dudtyr, "real(kind=8)", "DUDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtyl, d_dvdtyl, "real(kind=8)", "DVDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtyr, d_dvdtyr, "real(kind=8)", "DVDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtyl, d_dwdtyl, "real(kind=8)", "DWDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtyr, d_dwdtyr, "real(kind=8)", "DWDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtyl, d_dtdtyl, "real(kind=8)", "DTDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtyr, d_dtdtyr, "real(kind=8)", "DTDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtyl, d_dddtyl, "real(kind=8)", "DDDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtyr, d_dddtyr, "real(kind=8)", "DDDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouyl, d_acouyl, "real(kind=8)", "ACOUYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouyr, d_acouyr, "real(kind=8)", "ACOUYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2yl, d_ova2yl, "real(kind=8)", "OVA2YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2yr, d_ova2yr, "real(kind=8)", "OVA2YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1yl, d_gam1yl, "real(kind=8)", "GAM1YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1yr, d_gam1yr, "real(kind=8)", "GAM1YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmyl, d_ovgmyl, "real(kind=8)", "OVGMYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmyr, d_ovgmyr, "real(kind=8)", "OVGMYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtyl, d_sydtyl, "real(kind=8)", "SYDTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtyr, d_sydtyr, "real(kind=8)", "SYDTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpyl, d_sorpyl, "real(kind=8)", "SORPYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpyr, d_sorpyr, "real(kind=8)", "SORPYR")

!-----------------------------------------Boundary XY--------------------------------------------------------

    d_size = [nxglbl,nyglbl,1]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1zl, d_bcl1zl, "real(kind=8)", "BCL1ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1zr, d_bcl1zr, "real(kind=8)", "BCL1ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2zl, d_bcl2zl, "real(kind=8)", "BCL2ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2zr, d_bcl2zr, "real(kind=8)", "BCL2ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3zl, d_bcl3zl, "real(kind=8)", "BCL3ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3zr, d_bcl3zr, "real(kind=8)", "BCL3ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4zl, d_bcl4zl, "real(kind=8)", "BCL4ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4zr, d_bcl4zr, "real(kind=8)", "BCL4ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5zl, d_bcl5zl, "real(kind=8)", "BCL5ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5zr, d_bcl5zr, "real(kind=8)", "BCL5ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltzl, d_bcltzl, "real(kind=8)", "BCLTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltzr, d_bcltzr, "real(kind=8)", "BCLTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struzl, d_struzl, "real(kind=8)", "STRUZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struzr, d_struzr, "real(kind=8)", "STRUZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvzl, d_strvzl, "real(kind=8)", "STRVZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvzr, d_strvzr, "real(kind=8)", "STRVZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwzl, d_strwzl, "real(kind=8)", "STRWZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwzr, d_strwzr, "real(kind=8)", "STRWZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpzl, d_strpzl, "real(kind=8)", "STRPZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpzr, d_strpzr, "real(kind=8)", "STRPZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdzl, d_strdzl, "real(kind=8)", "STRDZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdzr, d_strdzr, "real(kind=8)", "STRDZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtzl, d_strtzl, "real(kind=8)", "STRTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtzr, d_strtzr, "real(kind=8)", "STRTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strezl, d_strezl, "real(kind=8)", "STREZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strezr, d_strezr, "real(kind=8)", "STREZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgzl, d_strgzl, "real(kind=8)", "STRGZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgzr, d_strgzr, "real(kind=8)", "STRGZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrzl, d_strrzl, "real(kind=8)", "STRRZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrzr, d_strrzr, "real(kind=8)", "STRRZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtzl, d_dudtzl, "real(kind=8)", "DUDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtzr, d_dudtzr, "real(kind=8)", "DUDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtzl, d_dvdtzl, "real(kind=8)", "DVDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtzr, d_dvdtzr, "real(kind=8)", "DVDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtzl, d_dwdtzl, "real(kind=8)", "DWDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtzr, d_dwdtzr, "real(kind=8)", "DWDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtzl, d_dtdtzl, "real(kind=8)", "DTDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtzr, d_dtdtzr, "real(kind=8)", "DTDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtzl, d_dddtzl, "real(kind=8)", "DDDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtzr, d_dddtzr, "real(kind=8)", "DDDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouzl, d_acouzl, "real(kind=8)", "ACOUZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouzr, d_acouzr, "real(kind=8)", "ACOUZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2zl, d_ova2zl, "real(kind=8)", "OVA2ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2zr, d_ova2zr, "real(kind=8)", "OVA2ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1zl, d_gam1zl, "real(kind=8)", "GAM1ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1zr, d_gam1zr, "real(kind=8)", "GAM1ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmzl, d_ovgmzl, "real(kind=8)", "OVGMZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmzr, d_ovgmzr, "real(kind=8)", "OVGMZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtzl, d_sydtzl, "real(kind=8)", "SYDTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtzr, d_sydtzr, "real(kind=8)", "SYDTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpzl, d_sorpzl, "real(kind=8)", "SORPZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpzr, d_sorpzr, "real(kind=8)", "SORPZR")

!------------------------------------Only X-direction--------------------------------------------------------

    d_size = [nxglbl,1,1]
    d_m    = [0,0,0]
    d_p    = [0,0,0]
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, crin, d_crin, "real(kind=8)", "CRIN")

!------------------------------------------------------------------------------------------------------------
    call ops_partition(" ")
!------------------------------------------------------------------------------------------------------------

END SUBROUTINE ops_data_init

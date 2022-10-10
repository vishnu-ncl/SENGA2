SUBROUTINE ops_data_init()
    use OPS_Fortran_Reference

    use data_types
    use com_senga
    use com_ops_senga

    INTEGER :: d_size(3)
    INTEGER d_base(3) /1,1,1/ !array indexing - start from 1
    INTEGER :: d_p(3) !max boundary depths for the dat in the possitive direction
    INTEGER :: d_m(3) !max boundary depths for the dat in the negative direction


    INTEGER a3d_000(3)    /0,0,0/

    integer stride3d_x(3) /1,0,0/
    integer stride3d_y(3) /0,1,0/
    integer stride3d_z(3) /0,0,1/
    
    integer stride3d_xy(3) /1,1,0/
    integer stride3d_xz(3) /1,0,1/
    integer stride3d_yz(3) /0,1,1/

    INTEGER a3d_000_to_p400_x(15) /0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0/
    INTEGER a3d_000_to_m400_x(15) /0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0/

    INTEGER a3d_000_to_p500_x(18) /0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0/
    INTEGER a3d_000_to_m500_x(18) /0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0, -5,0,0/

    INTEGER a3d_p100_to_p400_x(12) /1,0,0, 2,0,0, 3,0,0, 4,0,0/
    INTEGER a3d_m100_to_m400_x(12) /-1,0,0, -2,0,0, -3,0,0, -4,0,0/

    INTEGER a3d_p200_to_m200_x(15) /2,0,0, 1,0,0, 0,0,0, -1,0,0, -2,0,0/
    INTEGER a3d_p300_to_m300_x(21) /3,0,0, 2,0,0, 1,0,0, 0,0,0, -1,0,0, -2,0,0, -3,0,0/
    INTEGER a3d_p400_to_m400_x(27) /4,0,0, 3,0,0, 2,0,0, 1,0,0, 0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0/
    INTEGER a3d_p500_to_m500_x(33) /5,0,0, 4,0,0, 3,0,0, 2,0,0, 1,0,0, 0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0, -5,0,0/

    INTEGER a3d_p300_to_m100_x(15) /3,0,0, 2,0,0, 1,0,0, 0,0,0, -1,0,0/
    INTEGER a3d_p100_to_m300_x(15) /1,0,0, 0,0,0, -1,0,0, -2,0,0, -3,0,0/

    INTEGER a3d_p400_to_m100_x(18) /4,0,0, 3,0,0, 2,0,0, 1,0,0, 0,0,0, -1,0,0/
    INTEGER a3d_p100_to_m400_x(18) /1,0,0, 0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0/
!-------------------------------------------------------------------------------------
    INTEGER a3d_000_to_p040_y(15) /0,0,0, 0,1,0, 0,2,0, 0,3,0, 0,4,0/
    INTEGER a3d_000_to_m040_y(15) /0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0/

    INTEGER a3d_000_to_p050_y(18) /0,0,0, 0,1,0, 0,2,0, 0,3,0, 0,4,0, 0,5,0/
    INTEGER a3d_000_to_m050_y(18) /0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, 0,-5,0/

    INTEGER a3d_p010_to_p040_y(12) /0,1,0, 0,2,0, 0,3,0, 0,4,0/
    INTEGER a3d_m010_to_m040_y(12) /0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0/

    INTEGER a3d_p020_to_m020_y(15) /0,2,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0/
    INTEGER a3d_p030_to_m030_y(21) /0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0/
    INTEGER a3d_p040_to_m040_y(27) /0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0/
    INTEGER a3d_p050_to_m050_y(33) /0,5,0, 0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0, 0,-5,0/

    INTEGER a3d_p030_to_m010_y(15) /0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0/
    INTEGER a3d_p010_to_m030_y(15) /0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0/

    INTEGER a3d_p040_to_m010_y(18) /0,4,0, 0,3,0, 0,2,0, 0,1,0, 0,0,0, 0,-1,0/
    INTEGER a3d_p010_to_m040_y(18) /0,1,0, 0,0,0, 0,-1,0, 0,-2,0, 0,-3,0, 0,-4,0/
!-------------------------------------------------------------------------------------
    INTEGER a3d_000_to_p004_z(15) /0,0,0, 0,0,1, 0,0,2, 0,0,3, 0,0,4/
    INTEGER a3d_000_to_m004_z(15) /0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4/

    INTEGER a3d_000_to_p005_z(18) /0,0,0, 0,0,1, 0,0,2, 0,0,3, 0,0,4, 0,0,5/
    INTEGER a3d_000_to_m005_z(18) /0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, 0,0,-5/

    INTEGER a3d_p001_to_p004_z(12) /0,0,1, 0,0,2, 0,0,3, 0,0,4/
    INTEGER a3d_m001_to_m004_z(12) /0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4/

    INTEGER a3d_p002_to_m002_z(15) /0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2/
    INTEGER a3d_p003_to_m003_z(21) /0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3/
    INTEGER a3d_p004_to_m004_z(27) /0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4/
    INTEGER a3d_p005_to_m005_z(33) /0,0,5, 0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, 0,0,-5/

    INTEGER a3d_p003_to_m001_z(15) /0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1/
    INTEGER a3d_p001_to_m003_z(15) /0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3/

    INTEGER a3d_p004_to_m001_z(18) /0,0,4, 0,0,3, 0,0,2, 0,0,1, 0,0,0, 0,0,-1/
    INTEGER a3d_p001_to_m004_z(18) /0,0,1, 0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4/

!   *----------------------------OPS Declarations----------------------------*
!   Declare OPS Block
    call ops_decl_block(3, senga_grid, "senga grid")

!   Declare OPS Dats
    d_size = (/nxsize, nysize, nzsize/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store1, d_store1, "real(dp)", "STORE1")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store2, d_store2, "real(dp)", "STORE2")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store3, d_store3, "real(dp)", "STORE3")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store4, d_store4, "real(dp)", "STORE4")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store5, d_store5, "real(dp)", "STORE5")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store6, d_store6, "real(dp)", "STORE6")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, divm, d_divm, "real(dp)", "DIVM")
    
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ucor, d_ucor, "real(dp)", "UCOR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vcor, d_vcor, "real(dp)", "VCOR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wcor, d_wcor, "real(dp)", "WCOR")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd1x, d_wd1x, "real(dp)", "WD1X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd1x, d_pd1x, "real(dp)", "PD1X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td1x, d_td1x, "real(dp)", "TD1X")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd1y, d_wd1y, "real(dp)", "WD1Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd1y, d_pd1y, "real(dp)", "PD1Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td1y, d_td1y, "real(dp)", "TD1Y")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd1z, d_wd1z, "real(dp)", "WD1Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd1z, d_pd1z, "real(dp)", "PD1Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td1z, d_td1z, "real(dp)", "TD1Z")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd2x, d_wd2x, "real(dp)", "WD2X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd2x, d_pd2x, "real(dp)", "PD2X")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td2x, d_td2x, "real(dp)", "TD2X")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd2y, d_wd2y, "real(dp)", "WD2Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd2y, d_pd2y, "real(dp)", "PD2Y")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td2y, d_td2y, "real(dp)", "TD2Y")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd2z, d_wd2z, "real(dp)", "WD2Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd2z, d_pd2z, "real(dp)", "PD2Z")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td2z, d_td2z, "real(dp)", "TD2Z")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, drun, d_drun, "real(dp)", "DRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, urun, d_urun, "real(dp)", "URUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vrun, d_vrun, "real(dp)", "VRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wrun, d_wrun, "real(dp)", "WRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, erun, d_erun, "real(dp)", "ERUN")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, derr, d_derr, "real(dp)", "DERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, uerr, d_uerr, "real(dp)", "UERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, verr, d_verr, "real(dp)", "VERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, werr, d_werr, "real(dp)", "WERR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, eerr, d_eerr, "real(dp)", "EERR")

    d_size = (/nxsize, nysize, nzsize/)
    d_m = (/-nhalox,-nhaloy,-nhaloz/)
    d_p = (/nhalox,nhaloy,nhaloz/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, drhs, d_drhs, "real(dp)", "DRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, urhs, d_urhs, "real(dp)", "URHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vrhs, d_vrhs, "real(dp)", "VRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wrhs, d_wrhs, "real(dp)", "WRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, erhs, d_erhs, "real(dp)", "ERHS")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, utmp, d_utmp, "real(dp)", "UTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vtmp, d_vtmp, "real(dp)", "VTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wtmp, d_wtmp, "real(dp)", "WTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, prun, d_prun, "real(dp)", "PRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, trun, d_trun, "real(dp)", "TRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, transp, d_transp, "real(dp)", "TRANSP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store7, d_store7, "real(dp)", "STORE7")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wmomix, d_wmomix, "real(dp)", "WMOMIX")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, difmix, d_difmix, "real(dp)", "DIFMIX")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, tdrmix, d_tdrmix, "real(dp)", "TDRMIX")

!   Boundary
    d_size = (/1, nysize, nzsize/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1xl, d_bcl1xl, "real(dp)", "BCL1XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1xr, d_bcl1xr, "real(dp)", "BCL1XR")    
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2xl, d_bcl2xl, "real(dp)", "BCL2XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2xr, d_bcl2xr, "real(dp)", "BCL2XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3xl, d_bcl3xl, "real(dp)", "BCL3XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3xr, d_bcl3xr, "real(dp)", "BCL3XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4xl, d_bcl4xl, "real(dp)", "BCL4XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4xr, d_bcl4xr, "real(dp)", "BCL4XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5xl, d_bcl5xl, "real(dp)", "BCL5XL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5xr, d_bcl5xr, "real(dp)", "BCL5XR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltxl, d_bcltxl, "real(dp)", "BCLTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltxr, d_bcltxr, "real(dp)", "BCLTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struxl, d_struxl, "real(dp)", "STRUXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struxr, d_struxr, "real(dp)", "STRUXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvxl, d_strvxl, "real(dp)", "STRVXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvxr, d_strvxr, "real(dp)", "STRVXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwxl, d_strwxl, "real(dp)", "STRWXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwxr, d_strwxr, "real(dp)", "STRWXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpxl, d_strpxl, "real(dp)", "STRPXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpxr, d_strpxr, "real(dp)", "STRPXR")    
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdxl, d_strdxl, "real(dp)", "STRDXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdxr, d_strdxr, "real(dp)", "STRDXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtxl, d_strtxl, "real(dp)", "STRTXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtxr, d_strtxr, "real(dp)", "STRTXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strexl, d_strexl, "real(dp)", "STREXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strexr, d_strexr, "real(dp)", "STREXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgxl, d_strgxl, "real(dp)", "STRGXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgxr, d_strgxr, "real(dp)", "STRGXR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrxl, d_strrxl, "real(dp)", "STRRXL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrxr, d_strrxr, "real(dp)", "STRRXR")

    d_size = (/nxsize, 1, nzsize/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1yl, d_bcl1yl, "real(dp)", "BCL1YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1yr, d_bcl1yr, "real(dp)", "BCL1YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2yl, d_bcl2yl, "real(dp)", "BCL2YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2yr, d_bcl2yr, "real(dp)", "BCL2YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3yl, d_bcl3yl, "real(dp)", "BCL3YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3yr, d_bcl3yr, "real(dp)", "BCL3YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4yl, d_bcl4yl, "real(dp)", "BCL4YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4yr, d_bcl4yr, "real(dp)", "BCL4YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5yl, d_bcl5yl, "real(dp)", "BCL5YL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5yr, d_bcl5yr, "real(dp)", "BCL5YR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltyl, d_bcltyl, "real(dp)", "BCLTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltyr, d_bcltyr, "real(dp)", "BCLTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struyl, d_struyl, "real(dp)", "STRUYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struyr, d_struyr, "real(dp)", "STRUYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvyl, d_strvyl, "real(dp)", "STRVYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvyr, d_strvyr, "real(dp)", "STRVYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwyl, d_strwyl, "real(dp)", "STRWYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwyr, d_strwyr, "real(dp)", "STRWYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpyl, d_strpyl, "real(dp)", "STRPYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpyr, d_strpyr, "real(dp)", "STRPYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdyl, d_strdyl, "real(dp)", "STRDYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdyr, d_strdyr, "real(dp)", "STRDYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtyl, d_strtyl, "real(dp)", "STRTYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtyr, d_strtyr, "real(dp)", "STRTYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, streyl, d_streyl, "real(dp)", "STREYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, streyr, d_streyr, "real(dp)", "STREYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgyl, d_strgyl, "real(dp)", "STRGYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgyr, d_strgyr, "real(dp)", "STRGYR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strryl, d_strryl, "real(dp)", "STRRYL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strryr, d_strryr, "real(dp)", "STRRYR")
    
    d_size = (/nxsize, nysize, 1/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1zl, d_bcl1zl, "real(dp)", "BCL1ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1zr, d_bcl1zr, "real(dp)", "BCL1ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2zl, d_bcl2zl, "real(dp)", "BCL2ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2zr, d_bcl2zr, "real(dp)", "BCL2ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3zl, d_bcl3zl, "real(dp)", "BCL3ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3zr, d_bcl3zr, "real(dp)", "BCL3ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4zl, d_bcl4zl, "real(dp)", "BCL4ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4zr, d_bcl4zr, "real(dp)", "BCL4ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5zl, d_bcl5zl, "real(dp)", "BCL5ZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5zr, d_bcl5zr, "real(dp)", "BCL5ZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltzl, d_bcltzl, "real(dp)", "BCLTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltzr, d_bcltzr, "real(dp)", "BCLTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struzl, d_struzl, "real(dp)", "STRUZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struzr, d_struzr, "real(dp)", "STRUZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvzl, d_strvzl, "real(dp)", "STRVZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvzr, d_strvzr, "real(dp)", "STRVZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwzl, d_strwzl, "real(dp)", "STRWZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwzr, d_strwzr, "real(dp)", "STRWZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpzl, d_strpzl, "real(dp)", "STRPZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpzr, d_strpzr, "real(dp)", "STRPZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdzl, d_strdzl, "real(dp)", "STRDZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdzr, d_strdzr, "real(dp)", "STRDZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtzl, d_strtzl, "real(dp)", "STRTZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtzr, d_strtzr, "real(dp)", "STRTZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strezl, d_strezl, "real(dp)", "STREZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strezr, d_strezr, "real(dp)", "STREZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgzl, d_strgzl, "real(dp)", "STRGZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgzr, d_strgzr, "real(dp)", "STRGZR")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrzl, d_strrzl, "real(dp)", "STRRZL")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrzr, d_strrzr, "real(dp)", "STRRZR")

!   Declare OPS Stencil
    call ops_decl_stencil( 3, 1, a3d_000, s3d_000, "0,0,0")

    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_x, s3d_000_strid3d_x, "stride 3D X dir")
    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_y, s3d_000_strid3d_y, "stride 3D Y dir")
    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_z, s3d_000_strid3d_z, "stride 3D Z dir")

    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_xy, s3d_000_strid3d_xy, "stride 3D XY dir")
    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_xz, s3d_000_strid3d_xz, "stride 3D XZ dir")
    call ops_decl_strided_stencil( 3, 1, a3d_000, stride3d_yz, s3d_000_strid3d_yz, "stride 3D YZ dir")
!------------------------------------------------------------------------------------------------------
    call ops_decl_stencil( 3, 5, a3d_000_to_p400_x, s3d_000_to_p400_x, "0,0,0 to 4,0,0")
    call ops_decl_stencil( 3, 5, a3d_000_to_m400_x, s3d_000_to_m400_x, "0,0,0 to -4,0,0")
    
    call ops_decl_stencil( 3, 6, a3d_000_to_p500_x, s3d_000_to_p500_x, "0,0,0 to 5,0,0")
    call ops_decl_stencil( 3, 6, a3d_000_to_m500_x, s3d_000_to_m500_x, "0,0,0 to -5,0,0")

    call ops_decl_stencil( 3, 4, a3d_p100_to_p400_x, s3d_p100_to_p400_x, "1,0,0 to 4,0,0")
    call ops_decl_stencil( 3, 4, a3d_m100_to_m400_x, s3d_m100_to_m400_x, "-1,0,0 to -4,0,0")

    call ops_decl_stencil( 3,  5, a3d_p200_to_m200_x, s3d_p200_to_m200_x, "2,0,0 to -2,0,0")
    call ops_decl_stencil( 3,  7, a3d_p300_to_m300_x, s3d_p300_to_m300_x, "3,0,0 to -3,0,0")
    call ops_decl_stencil( 3,  9, a3d_p400_to_m400_x, s3d_p400_to_m400_x, "4,0,0 to -4,0,0")
    call ops_decl_stencil( 3, 11, a3d_p500_to_m500_x, s3d_p500_to_m500_x, "5,0,0 to -5,0,0")

    call ops_decl_stencil( 3,  5, a3d_p300_to_m100_x, s3d_p300_to_m100_x, "3,0,0 to -1,0,0")
    call ops_decl_stencil( 3,  5, a3d_p100_to_m300_x, s3d_p100_to_m300_x, "1,0,0 to -3,0,0")

    call ops_decl_stencil( 3,  6, a3d_p400_to_m100_x, s3d_p400_to_m100_x, "4,0,0 to -1,0,0")
    call ops_decl_stencil( 3,  6, a3d_p100_to_m400_x, s3d_p100_to_m400_x, "1,0,0 to -4,0,0")
!-----------------------------------------------------------------------------------------------------
    call ops_decl_stencil( 3, 5, a3d_000_to_p040_y, s3d_000_to_p040_y, "0,0,0 to 0,4,0")
    call ops_decl_stencil( 3, 5, a3d_000_to_m040_y, s3d_000_to_m040_y, "0,0,0 to  0,-4,0")

    call ops_decl_stencil( 3, 6, a3d_000_to_p050_y, s3d_000_to_p050_y, "0,0,0 to 0,5,0")
    call ops_decl_stencil( 3, 6, a3d_000_to_m050_y, s3d_000_to_m050_y, "0,0,0 to  0,-5,0")

    call ops_decl_stencil( 3, 4, a3d_p010_to_p040_y, s3d_p010_to_p040_y, "0,1,0 to 0,4,0")
    call ops_decl_stencil( 3, 4, a3d_m010_to_m040_y, s3d_m010_to_m040_y, "0,-1,0 to 0,-4,0")

    call ops_decl_stencil( 3,  5, a3d_p020_to_m020_y, s3d_p020_to_m020_y, "0,2,0 to  0,-2,0")
    call ops_decl_stencil( 3,  7, a3d_p030_to_m030_y, s3d_p030_to_m030_y, "0,3,0 to  0,-3,0")
    call ops_decl_stencil( 3,  9, a3d_p040_to_m040_y, s3d_p040_to_m040_y, "0,4,0 to  0,-4,0")
    call ops_decl_stencil( 3, 11, a3d_p050_to_m050_y, s3d_p050_to_m050_y, "0,5,0 to  0,-5,0")

    call ops_decl_stencil( 3,  5, a3d_p030_to_m010_y, s3d_p030_to_m010_y, "0,3,0 to  0,-1,0")
    call ops_decl_stencil( 3,  5, a3d_p010_to_m030_y, s3d_p010_to_m030_y, "0,1,0 to  0,-3,0")

    call ops_decl_stencil( 3,  6, a3d_p040_to_m010_y, s3d_p040_to_m010_y, "0,4,0 to  0,-1,0")
    call ops_decl_stencil( 3,  6, a3d_p010_to_m040_y, s3d_p010_to_m040_y, "0,1,0 to  0,-4,0")
!-------------------------------------------------------------------------------------
    call ops_decl_stencil( 3, 5, a3d_000_to_p004_z, s3d_000_to_p004_z, "0,0,0 to 0,0,4")
    call ops_decl_stencil( 3, 5, a3d_000_to_m004_z, s3d_000_to_m004_z, "0,0,0 to  0,0,-4")

    call ops_decl_stencil( 3, 6, a3d_000_to_p005_z, s3d_000_to_p005_z, "0,0,0 to 0,0,5")
    call ops_decl_stencil( 3, 6, a3d_000_to_m005_z, s3d_000_to_m005_z, "0,0,0 to  0,0,-5")

    call ops_decl_stencil( 3, 4, a3d_p001_to_p004_z, s3d_p001_to_p004_z, "0,0,1 to 0,0,4")
    call ops_decl_stencil( 3, 4, a3d_m001_to_m004_z, s3d_m001_to_m004_z, "0,0,-1 to 0,0,-4")

    call ops_decl_stencil( 3,  5, a3d_p002_to_m002_z, s3d_p002_to_m002_z, "0,0,2 to  0,0,-2")
    call ops_decl_stencil( 3,  7, a3d_p003_to_m003_z, s3d_p003_to_m003_z, "0,0,3 to  0,0,-3")
    call ops_decl_stencil( 3,  9, a3d_p004_to_m004_z, s3d_p004_to_m004_z, "0,0,4 to  0,0,-4")
    call ops_decl_stencil( 3, 11, a3d_p005_to_m005_z, s3d_p005_to_m005_z, "0,0,5 to  0,0,-5")

    call ops_decl_stencil( 3,  5, a3d_p003_to_m001_z, s3d_p003_to_m001_z, "0,0,3 to  0,0,-1")
    call ops_decl_stencil( 3,  5, a3d_p001_to_m003_z, s3d_p001_to_m003_z, "0,0,1 to  0,0,-3")

    call ops_decl_stencil( 3,  6, a3d_p004_to_m001_z, s3d_p004_to_m001_z, "0,0,4 to  0,0,-1")
    call ops_decl_stencil( 3,  6, a3d_p001_to_m004_z, s3d_p001_to_m004_z, "0,0,1 to  0,0,-4")

    call ops_partition(" ")

END SUBROUTINE ops_data_init

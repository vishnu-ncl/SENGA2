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

!-------------------------------------------------------------------------------------
    INTEGER a3d_p320_m120_mixed_xy(60) /3,2,0, 3,1,0, 3,-1,0, 3,-2,0, 2,2,0, 2,1,0, 2,-1,0, 2,-2,0, 1,2,0, 1,1,0, 1,-1,0, 1,-2,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0, -1,2,0, -1,1,0, -1,-1,0, -1,-2,0/
    INTEGER a3d_p120_m320_mixed_xy(60) /1,2,0, 1,1,0, 1,-1,0, 1,-2,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0, -1,2,0, -1,1,0, -1,-1,0, -1,-2,0, -2,2,0, -2,1,0, -2,-1,0, -2,-2,0, -3,2,0, -3,1,0, -3,-1,0, -3,-2,0/

    INTEGER a3d_p420_m020_mixed_xy(60) /4,2,0, 4,1,0, 4,-1,0, 4,-2,0, 3,2,0, 3,1,0, 3,-1,0, 3,-2,0, 2,2,0, 2,1,0, 2,-1,0, 2,-2,0, 1,2,0, 1,1,0, 1,-1,0, 1,-2,0, 0,2,0, 0,1,0, 0,-1,0, 0,-2,0/
    INTEGER a3d_p020_m420_mixed_xy(60) /0,2,0, 0,1,0, 0,-1,0, 0,-2,0, -1,2,0, -1,1,0, -1,-1,0, -1,-2,0, -2,2,0, -2,1,0, -2,-1,0, -2,-2,0, -3,2,0, -3,1,0, -3,-1,0, -3,-2,0, -4,2,0, -4,1,0, -4,-1,0, -4,-2,0/

    INTEGER a3d_p220_m220_mixed_xy(24) /2,2,0, 2,-2,0, 1,1,0, 1,-1,0, -1,1,0, -1,-1,0, -2,2,0, -2,-2,0/
    INTEGER a3d_p330_m330_mixed_xy(36) /3,3,0, 3,-3,0, 2,2,0, 2,-2,0, 1,1,0, 1,-1,0, -1,1,0, -1,-1,0, -2,2,0, -2,-2,0, -3,3,0, -3,-3,0/
    INTEGER a3d_p440_m440_mixed_xy(48) /4,4,0, 4,-4,0, 3,3,0, 3,-3,0, 2,2,0, 2,-2,0, 1,1,0, 1,-1,0, -1,1,0, -1,-1,0, -2,2,0, -2,-2,0, -3,3,0, -3,-3,0, -4,4,0, -4,-4,0/

    INTEGER a3d_p550_to_m550_xy(33) /5,5,0, 4,4,0, 3,3,0, 2,2,0, 1,1,0, 0,0,0, -1,-1,0, -2,-2,0, -3,-3,0, -4,-4,0, -5,-5,0/


!   *----------------------------OPS Declarations----------------------------*
!   Declare OPS Block
    call ops_decl_block(3, senga_grid, "senga grid")

!   Declare OPS Dats
    d_size = (/nxsize, nysize, nzsize/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)

    allocate (store1(nxsize,nysize,nzsize))    
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store1, d_store1, "real(8)", "STORE1")
    allocate (store2(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store2, d_store2, "real(8)", "STORE2")
    allocate (store3(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store3, d_store3, "real(8)", "STORE3")
    allocate (store4(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store4, d_store4, "real(8)", "STORE4")
    allocate (store5(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store5, d_store5, "real(8)", "STORE5")
    allocate (store6(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store6, d_store6, "real(8)", "STORE6")
    allocate (divm(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, divm, d_divm, "real(8)", "DIVM")
    
    allocate (ucor(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ucor, d_ucor, "real(8)", "UCOR")
    allocate (vcor(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vcor, d_vcor, "real(8)", "VCOR")
    allocate (wcor(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wcor, d_wcor, "real(8)", "WCOR")

    allocate (wd1x(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd1x, d_wd1x, "real(8)", "WD1X")
    allocate (pd1x(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd1x, d_pd1x, "real(8)", "PD1X")
    allocate (td1x(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td1x, d_td1x, "real(8)", "TD1X")

    allocate (wd1y(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd1y, d_wd1y, "real(8)", "WD1Y")
    allocate (pd1y(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd1y, d_pd1y, "real(8)", "PD1Y")
    allocate (td1y(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td1y, d_td1y, "real(8)", "TD1Y")

    allocate (wd1z(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd1z, d_wd1z, "real(8)", "WD1Z")
    allocate (pd1z(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd1z, d_pd1z, "real(8)", "PD1Z")
    allocate (td1z(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td1z, d_td1z, "real(8)", "TD1Z")

    allocate (wd2x(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd2x, d_wd2x, "real(8)", "WD2X")
    allocate (pd2x(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd2x, d_pd2x, "real(8)", "PD2X")
    allocate (td2x(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td2x, d_td2x, "real(8)", "TD2X")

    allocate (wd2y(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd2y, d_wd2y, "real(8)", "WD2Y")
    allocate (pd2y(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd2y, d_pd2y, "real(8)", "PD2Y")
    allocate (td2y(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td2y, d_td2y, "real(8)", "TD2Y")

    allocate (wd2z(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wd2z, d_wd2z, "real(8)", "WD2Z")
    allocate (pd2z(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, pd2z, d_pd2z, "real(8)", "PD2Z")
    allocate (td2z(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, td2z, d_td2z, "real(8)", "TD2Z")

    allocate (ufxl(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ufxl, d_ufxl, "real(8)", "UFXL")
    allocate (vfxl(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vfxl, d_vfxl, "real(8)", "VFXL")
    allocate (wfxl(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wfxl, d_wfxl, "real(8)", "WFXL")

    allocate (drun(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, drun, d_drun, "real(8)", "DRUN")
    allocate (urun(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, urun, d_urun, "real(8)", "URUN")
    allocate (vrun(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vrun, d_vrun, "real(8)", "VRUN")
    allocate (wrun(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wrun, d_wrun, "real(8)", "WRUN")
    allocate (erun(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, erun, d_erun, "real(8)", "ERUN")

    allocate (derr(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, derr, d_derr, "real(8)", "DERR")
    allocate (uerr(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, uerr, d_uerr, "real(8)", "UERR")
    allocate (verr(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, verr, d_verr, "real(8)", "VERR")
    allocate (werr(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, werr, d_werr, "real(8)", "WERR")
    allocate (eerr(nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, eerr, d_eerr, "real(8)", "EERR")

!---------------------------------------MULTI-DIM DAT-------------------------------------------------------
    d_size = (/nxsize, nysize, nzsize/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    allocate (yrun(nspcmx,nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, yrun, d_yrun, "real(8)", "YRUN")
    allocate (yerr(nspcmx,nxsize,nysize,nzsize)) 
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, yerr, d_yerr, "real(8)", "YERR")
    allocate (rate(nspcmx,nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, rate, d_rate, "real(8)", "RATE")
    allocate (rrte(nspcmx,nxsize,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, rrte, d_rrte, "real(8)", "RRTE")

    d_size = (/nxsize, nysize, nzsize/)
    d_m = (/-nhalox,-nhaloy,-nhaloz/)
    d_p = (/nhalox,nhaloy,nhaloz/)
    allocate (itndex(nintmx,1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, nintmx, d_size, d_base, d_m, d_p, itndex, d_itndex, "integer", "ITNDEX")
    allocate (yrhs(nspcmx,1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, yrhs, d_yrhs, "real(8)", "YRHS")

    d_size = (/1, nysize, nzsize/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    allocate (bclyxl(nspcmx,1,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, bclyxl, d_bclyxl, "real(8)", "BCLYXL")
    allocate (bclyxr(nspcmx,1,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, bclyxr, d_bclyxr, "real(8)", "BCLYXR")
    allocate (stryxl(nspcmx,1,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, stryxl, d_stryxl, "real(8)", "STRYXL")
    allocate (stryxr(nspcmx,1,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, stryxr, d_stryxr, "real(8)", "STRYXR")
    allocate (dydtxl(nspcmx,1,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, dydtxl, d_dydtxl, "real(8)", "DYDTXL")
    allocate (dydtxr(nspcmx,1,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, dydtxr, d_dydtxr, "real(8)", "DYDTXR")
    allocate (ratexl(nspcmx,1,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, ratexl, d_ratexl, "real(8)", "RATEXL")
    allocate (ratexr(nspcmx,1,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, ratexr, d_ratexr, "real(8)", "RATEXR")
    allocate (strhxl(nspcmx,1,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, strhxl, d_strhxl, "real(8)", "STRHXL")
    allocate (strhxr(nspcmx,1,nysize,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, strhxr, d_strhxr, "real(8)", "STRHXR")

    d_size = (/nxsize, 1, nzsize/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    allocate (bclyyl(nspcmx,nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, bclyyl, d_bclyyl, "real(8)", "BCLYYL")
    allocate (bclyyr(nspcmx,nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, bclyyr, d_bclyyr, "real(8)", "BCLYYR")
    allocate (stryyl(nspcmx,nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, stryyl, d_stryyl, "real(8)", "STRYYL")
    allocate (stryyr(nspcmx,nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, stryyr, d_stryyr, "real(8)", "STRYYR")
    allocate (dydtyl(nspcmx,nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, dydtyl, d_dydtyl, "real(8)", "DYDTYL")
    allocate (dydtyr(nspcmx,nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, dydtyr, d_dydtyr, "real(8)", "DYDTYR")
    allocate (rateyl(nspcmx,nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, rateyl, d_rateyl, "real(8)", "RATEYL")
    allocate (rateyr(nspcmx,nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, rateyr, d_rateyr, "real(8)", "RATEYR")
    allocate (strhyl(nspcmx,nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, strhyl, d_strhyl, "real(8)", "STRHYL")
    allocate (strhyr(nspcmx,nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, strhyr, d_strhyr, "real(8)", "STRHYR")

    d_size = (/nxsize, nysize, 1/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    allocate (bclyzl(nspcmx,nxsize,nysize,1))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, bclyzl, d_bclyzl, "real(8)", "BCLYZL")
    allocate (bclyzr(nspcmx,nxsize,nysize,1))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, bclyzr, d_bclyzr, "real(8)", "BCLYZR")
    allocate (stryzl(nspcmx,nxsize,nysize,1))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, stryzl, d_stryzl, "real(8)", "STRYZL")
    allocate (stryzr(nspcmx,nxsize,nysize,1))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, stryzr, d_stryzr, "real(8)", "STRYZR")
    allocate (dydtzl(nspcmx,nxsize,nysize,1))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, dydtzl, d_dydtzl, "real(8)", "DYDTZL")
    allocate (dydtzr(nspcmx,nxsize,nysize,1))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, dydtzr, d_dydtzr, "real(8)", "DYDTZR")
    allocate (ratezl(nspcmx,nxsize,nysize,1))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, ratezl, d_ratezl, "real(8)", "RATEZL")
    allocate (ratezr(nspcmx,nxsize,nysize,1))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, ratezr, d_ratezr, "real(8)", "RATEZR")
    allocate (strhzl(nspcmx,nxsize,nysize,1))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, strhzl, d_strhzl, "real(8)", "STRHZL")
    allocate (strhzr(nspcmx,nxsize,nysize,1))
    call ops_decl_dat(senga_grid, nspcmx, d_size, d_base, d_m, d_p, strhzr, d_strhzr, "real(8)", "STRHZR")

!---------------------------------------WITH HALOS---------------------------------------------------------
    d_size = (/nxsize, nysize, nzsize/)
    d_m = (/-nhalox,-nhaloy,-nhaloz/)
    d_p = (/nhalox,nhaloy,nhaloz/)
    allocate (drhs(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, drhs, d_drhs, "real(8)", "DRHS")
    allocate (urhs(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, urhs, d_urhs, "real(8)", "URHS")
    allocate (vrhs(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vrhs, d_vrhs, "real(8)", "VRHS")
    allocate (wrhs(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wrhs, d_wrhs, "real(8)", "WRHS")
    allocate (erhs(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, erhs, d_erhs, "real(8)", "ERHS")

    allocate (utmp(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, utmp, d_utmp, "real(8)", "UTMP")
    allocate (vtmp(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vtmp, d_vtmp, "real(8)", "VTMP")
    allocate (wtmp(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wtmp, d_wtmp, "real(8)", "WTMP")
    allocate (prun(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, prun, d_prun, "real(8)", "PRUN")
    allocate (trun(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, trun, d_trun, "real(8)", "TRUN")
    allocate (transp(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, transp, d_transp, "real(8)", "TRANSP")
    allocate (store7(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store7, d_store7, "real(8)", "STORE7")

    allocate (wmomix(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wmomix, d_wmomix, "real(8)", "WMOMIX")
    allocate (difmix(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, difmix, d_difmix, "real(8)", "DIFMIX")
    allocate (tdrmix(1-nhalox:nxsize+nhalox,1-nhaloy:nysize+nhaloy,1-nhaloz:nzsize+nhaloz))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, tdrmix, d_tdrmix, "real(8)", "TDRMIX")

!-----------------------------------------Boundary YZ--------------------------------------------------
    d_size = (/1, nysize, nzsize/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    allocate (bcl1xl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1xl, d_bcl1xl, "real(8)", "BCL1XL")
    allocate (bcl1xr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1xr, d_bcl1xr, "real(8)", "BCL1XR")    
    allocate (bcl2xl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2xl, d_bcl2xl, "real(8)", "BCL2XL")
    allocate (bcl2xr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2xr, d_bcl2xr, "real(8)", "BCL2XR")
    allocate (bcl3xl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3xl, d_bcl3xl, "real(8)", "BCL3XL")
    allocate (bcl3xr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3xr, d_bcl3xr, "real(8)", "BCL3XR")
    allocate (bcl4xl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4xl, d_bcl4xl, "real(8)", "BCL4XL")
    allocate (bcl4xr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4xr, d_bcl4xr, "real(8)", "BCL4XR")
    allocate (bcl5xl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5xl, d_bcl5xl, "real(8)", "BCL5XL")
    allocate (bcl5xr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5xr, d_bcl5xr, "real(8)", "BCL5XR")
    allocate (bcltxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltxl, d_bcltxl, "real(8)", "BCLTXL")
    allocate (bcltxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltxr, d_bcltxr, "real(8)", "BCLTXR")
    allocate (struxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struxl, d_struxl, "real(8)", "STRUXL")
    allocate (struxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struxr, d_struxr, "real(8)", "STRUXR")
    allocate (strvxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvxl, d_strvxl, "real(8)", "STRVXL")
    allocate (strvxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvxr, d_strvxr, "real(8)", "STRVXR")
    allocate (strwxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwxl, d_strwxl, "real(8)", "STRWXL")
    allocate (strwxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwxr, d_strwxr, "real(8)", "STRWXR")
    allocate (strpxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpxl, d_strpxl, "real(8)", "STRPXL")
    allocate (strpxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpxr, d_strpxr, "real(8)", "STRPXR")    
    allocate (strdxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdxl, d_strdxl, "real(8)", "STRDXL")
    allocate (strdxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdxr, d_strdxr, "real(8)", "STRDXR")
    allocate (strtxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtxl, d_strtxl, "real(8)", "STRTXL")
    allocate (strtxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtxr, d_strtxr, "real(8)", "STRTXR")
    allocate (strexl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strexl, d_strexl, "real(8)", "STREXL")
    allocate (strexr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strexr, d_strexr, "real(8)", "STREXR")
    allocate (strgxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgxl, d_strgxl, "real(8)", "STRGXL")
    allocate (strgxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgxr, d_strgxr, "real(8)", "STRGXR")
    allocate (strrxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrxl, d_strrxl, "real(8)", "STRRXL")
    allocate (strrxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrxr, d_strrxr, "real(8)", "STRRXR")
    allocate (dudtxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtxl, d_dudtxl, "real(8)", "DUDTXL")
    allocate (dudtxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtxr, d_dudtxr, "real(8)", "DUDTXR")
    allocate (dvdtxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtxl, d_dvdtxl, "real(8)", "DVDTXL")
    allocate (dvdtxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtxr, d_dvdtxr, "real(8)", "DVDTXR")
    allocate (dwdtxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtxl, d_dwdtxl, "real(8)", "DWDTXL")
    allocate (dwdtxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtxr, d_dwdtxr, "real(8)", "DWDTXR")
    allocate (dtdtxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtxl, d_dtdtxl, "real(8)", "DTDTXL")
    allocate (dtdtxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtxr, d_dtdtxr, "real(8)", "DTDTXR")
    allocate (dddtxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtxl, d_dddtxl, "real(8)", "DDDTXL")
    allocate (dddtxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtxr, d_dddtxr, "real(8)", "DDDTXR")
    allocate (acouxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouxl, d_acouxl, "real(8)", "ACOUXL")
    allocate (acouxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouxr, d_acouxr, "real(8)", "ACOUXR")
    allocate (ova2xl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2xl, d_ova2xl, "real(8)", "OVA2XL")
    allocate (ova2xr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2xr, d_ova2xr, "real(8)", "OVA2XR")
    allocate (gam1xl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1xl, d_gam1xl, "real(8)", "GAM1XL")
    allocate (gam1xr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1xr, d_gam1xr, "real(8)", "GAM1XR")
    allocate (ovgmxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmxl, d_ovgmxl, "real(8)", "OVGMXL")
    allocate (ovgmxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmxr, d_ovgmxr, "real(8)", "OVGMXR")
    allocate (sydtxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtxl, d_sydtxl, "real(8)", "SYDTXL")
    allocate (sydtxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtxr, d_sydtxr, "real(8)", "SYDTXR")
    allocate (sorpxl(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpxl, d_sorpxl, "real(8)", "SORPXL")
    allocate (sorpxr(1,nysize,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpxr, d_sorpxr, "real(8)", "SORPXR")


!-----------------------------------------Boundary XZ--------------------------------------------------
    d_size = (/nxsize, 1, nzsize/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    allocate (bcl1yl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1yl, d_bcl1yl, "real(8)", "BCL1YL")
    allocate (bcl1yr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1yr, d_bcl1yr, "real(8)", "BCL1YR")
    allocate (bcl2yl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2yl, d_bcl2yl, "real(8)", "BCL2YL")
    allocate (bcl2yr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2yr, d_bcl2yr, "real(8)", "BCL2YR")
    allocate (bcl3yl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3yl, d_bcl3yl, "real(8)", "BCL3YL")
    allocate (bcl3yr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3yr, d_bcl3yr, "real(8)", "BCL3YR")
    allocate (bcl4yl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4yl, d_bcl4yl, "real(8)", "BCL4YL")
    allocate (bcl4yr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4yr, d_bcl4yr, "real(8)", "BCL4YR")
    allocate (bcl5yl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5yl, d_bcl5yl, "real(8)", "BCL5YL")
    allocate (bcl5yr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5yr, d_bcl5yr, "real(8)", "BCL5YR")
    allocate (bcltyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltyl, d_bcltyl, "real(8)", "BCLTYL")
    allocate (bcltyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltyr, d_bcltyr, "real(8)", "BCLTYR")
    allocate (struyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struyl, d_struyl, "real(8)", "STRUYL")
    allocate (struyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struyr, d_struyr, "real(8)", "STRUYR")
    allocate (strvyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvyl, d_strvyl, "real(8)", "STRVYL")
    allocate (strvyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvyr, d_strvyr, "real(8)", "STRVYR")
    allocate (strwyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwyl, d_strwyl, "real(8)", "STRWYL")
    allocate (strwyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwyr, d_strwyr, "real(8)", "STRWYR")
    allocate (strpyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpyl, d_strpyl, "real(8)", "STRPYL")
    allocate (strpyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpyr, d_strpyr, "real(8)", "STRPYR")
    allocate (strdyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdyl, d_strdyl, "real(8)", "STRDYL")
    allocate (strdyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdyr, d_strdyr, "real(8)", "STRDYR")
    allocate (strtyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtyl, d_strtyl, "real(8)", "STRTYL")
    allocate (strtyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtyr, d_strtyr, "real(8)", "STRTYR")
    allocate (streyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, streyl, d_streyl, "real(8)", "STREYL")
    allocate (streyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, streyr, d_streyr, "real(8)", "STREYR")
    allocate (strgyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgyl, d_strgyl, "real(8)", "STRGYL")
    allocate (strgyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgyr, d_strgyr, "real(8)", "STRGYR")
    allocate (strryl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strryl, d_strryl, "real(8)", "STRRYL")
    allocate (strryr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strryr, d_strryr, "real(8)", "STRRYR")
    allocate (dudtyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtyl, d_dudtyl, "real(8)", "DUDTYL")
    allocate (dudtyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtyr, d_dudtyr, "real(8)", "DUDTYR")
    allocate (dvdtyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtyl, d_dvdtyl, "real(8)", "DVDTYL")
    allocate (dvdtyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtyr, d_dvdtyr, "real(8)", "DVDTYR")
    allocate (dwdtyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtyl, d_dwdtyl, "real(8)", "DWDTYL")
    allocate (dwdtyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtyr, d_dwdtyr, "real(8)", "DWDTYR")
    allocate (dtdtyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtyl, d_dtdtyl, "real(8)", "DTDTYL")
    allocate (dtdtyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtyr, d_dtdtyr, "real(8)", "DTDTYR")
    allocate (dddtyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtyl, d_dddtyl, "real(8)", "DDDTYL")
    allocate (dddtyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtyr, d_dddtyr, "real(8)", "DDDTYR")
    allocate (acouyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouyl, d_acouyl, "real(8)", "ACOUYL")
    allocate (acouyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouyr, d_acouyr, "real(8)", "ACOUYR")
    allocate (ova2yl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2yl, d_ova2yl, "real(8)", "OVA2YL")
    allocate (ova2yr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2yr, d_ova2yr, "real(8)", "OVA2YR")
    allocate (gam1yl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1yl, d_gam1yl, "real(8)", "GAM1YL")
    allocate (gam1yr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1yr, d_gam1yr, "real(8)", "GAM1YR")
    allocate (ovgmyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmyl, d_ovgmyl, "real(8)", "OVGMYL")
    allocate (ovgmyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmyr, d_ovgmyr, "real(8)", "OVGMYR")
    allocate (sydtyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtyl, d_sydtyl, "real(8)", "SYDTYL")
    allocate (sydtyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtyr, d_sydtyr, "real(8)", "SYDTYR")
    allocate (sorpyl(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpyl, d_sorpyl, "real(8)", "SORPYL")
    allocate (sorpyr(nxsize,1,nzsize))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpyr, d_sorpyr, "real(8)", "SORPYR")


!-----------------------------------------Boundary XY--------------------------------------------------
    d_size = (/nxsize, nysize, 1/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    allocate (bcl1zl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1zl, d_bcl1zl, "real(8)", "BCL1ZL")
    allocate (bcl1zr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl1zr, d_bcl1zr, "real(8)", "BCL1ZR")
    allocate (bcl2zl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2zl, d_bcl2zl, "real(8)", "BCL2ZL")
    allocate (bcl2zr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl2zr, d_bcl2zr, "real(8)", "BCL2ZR")
    allocate (bcl3zl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3zl, d_bcl3zl, "real(8)", "BCL3ZL")
    allocate (bcl3zr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl3zr, d_bcl3zr, "real(8)", "BCL3ZR")
    allocate (bcl4zl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4zl, d_bcl4zl, "real(8)", "BCL4ZL")
    allocate (bcl4zr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl4zr, d_bcl4zr, "real(8)", "BCL4ZR")
    allocate (bcl5zl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5zl, d_bcl5zl, "real(8)", "BCL5ZL")
    allocate (bcl5zr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcl5zr, d_bcl5zr, "real(8)", "BCL5ZR")
    allocate (bcltzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltzl, d_bcltzl, "real(8)", "BCLTZL")
    allocate (bcltzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, bcltzr, d_bcltzr, "real(8)", "BCLTZR")
    allocate (struzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struzl, d_struzl, "real(8)", "STRUZL")
    allocate (struzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, struzr, d_struzr, "real(8)", "STRUZR")
    allocate (strvzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvzl, d_strvzl, "real(8)", "STRVZL")
    allocate (strvzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strvzr, d_strvzr, "real(8)", "STRVZR")
    allocate (strwzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwzl, d_strwzl, "real(8)", "STRWZL")
    allocate (strwzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strwzr, d_strwzr, "real(8)", "STRWZR")
    allocate (strpzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpzl, d_strpzl, "real(8)", "STRPZL")
    allocate (strpzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strpzr, d_strpzr, "real(8)", "STRPZR")
    allocate (strdzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdzl, d_strdzl, "real(8)", "STRDZL")
    allocate (strdzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strdzr, d_strdzr, "real(8)", "STRDZR")
    allocate (strtzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtzl, d_strtzl, "real(8)", "STRTZL")
    allocate (strtzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strtzr, d_strtzr, "real(8)", "STRTZR")
    allocate (strezl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strezl, d_strezl, "real(8)", "STREZL")
    allocate (strezr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strezr, d_strezr, "real(8)", "STREZR")
    allocate (strgzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgzl, d_strgzl, "real(8)", "STRGZL")
    allocate (strgzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strgzr, d_strgzr, "real(8)", "STRGZR")
    allocate (strrzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrzl, d_strrzl, "real(8)", "STRRZL")
    allocate (strrzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, strrzr, d_strrzr, "real(8)", "STRRZR")
    allocate (dudtzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtzl, d_dudtzl, "real(8)", "DUDTZL")
    allocate (dudtzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dudtzr, d_dudtzr, "real(8)", "DUDTZR")
    allocate (dvdtzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtzl, d_dvdtzl, "real(8)", "DVDTZL")
    allocate (dvdtzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dvdtzr, d_dvdtzr, "real(8)", "DVDTZR")
    allocate (dwdtzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtzl, d_dwdtzl, "real(8)", "DWDTZL")
    allocate (dwdtzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dwdtzr, d_dwdtzr, "real(8)", "DWDTZR")
    allocate (dtdtzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtzl, d_dtdtzl, "real(8)", "DTDTZL")
    allocate (dtdtzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dtdtzr, d_dtdtzr, "real(8)", "DTDTZR")
    allocate (dddtzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtzl, d_dddtzl, "real(8)", "DDDTZL")
    allocate (dddtzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, dddtzr, d_dddtzr, "real(8)", "DDDTZR")
    allocate (acouzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouzl, d_acouzl, "real(8)", "ACOUZL")
    allocate (acouzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, acouzr, d_acouzr, "real(8)", "ACOUZR")
    allocate (ova2zl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2zl, d_ova2zl, "real(8)", "OVA2ZL")
    allocate (ova2zr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ova2zr, d_ova2zr, "real(8)", "OVA2ZR")
    allocate (gam1zl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1zl, d_gam1zl, "real(8)", "GAM1ZL")
    allocate (gam1zr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, gam1zr, d_gam1zr, "real(8)", "GAM1ZR")
    allocate (ovgmzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmzl, d_ovgmzl, "real(8)", "OVGMZL")
    allocate (ovgmzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, ovgmzr, d_ovgmzr, "real(8)", "OVGMZR")
    allocate (sydtzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtzl, d_sydtzl, "real(8)", "SYDTZL")
    allocate (sydtzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sydtzr, d_sydtzr, "real(8)", "SYDTZR")
    allocate (sorpzl(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpzl, d_sorpzl, "real(8)", "SORPZL")
    allocate (sorpzr(nxsize,nysize,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, sorpzr, d_sorpzr, "real(8)", "SORPZR")

!------------------------------------Only X-direction-------------------------------------------------
    d_size = (/nxsize, 1, 1/)
    d_m = (/0,0,0/)
    d_p = (/0,0,0/)
    allocate (crin(nxsize,1,1))
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, crin, d_crin, "real(8)", "CRIN")

!------------------------------------OPS Reduction Handles--------------------------------------------
    call ops_decl_reduction_handle(8, h_erdtot, "real(8)", "erdtot")
    call ops_decl_reduction_handle(8, h_erutot, "real(8)", "erutot")
    call ops_decl_reduction_handle(8, h_ervtot, "real(8)", "ervtot")
    call ops_decl_reduction_handle(8, h_erwtot, "real(8)", "erwtot")
    call ops_decl_reduction_handle(8, h_eretot, "real(8)", "eretot")
    call ops_decl_reduction_handle(8, h_erytot, "real(8)", "erytot")
    call ops_decl_reduction_handle(8, h_prefer, "real(8)", "prefer")
    call ops_decl_reduction_handle(8, h_tket, "real(8)", "tket")
    call ops_decl_reduction_handle(8, h_ubart, "real(8)", "ubart")
    call ops_decl_reduction_handle(8, h_vbart, "real(8)", "vbart")
    call ops_decl_reduction_handle(8, h_wbart, "real(8)", "wbart")
    call ops_decl_reduction_handle(8, h_uvart, "real(8)", "uvart")
    call ops_decl_reduction_handle(8, h_vvart, "real(8)", "vvart")
    call ops_decl_reduction_handle(8, h_wvart, "real(8)", "wvart")

!------------------------------------OPS Stencil---------------------------------------------------
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
!-------------------------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------------------------
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
!--------------------------------------------------------------------------------------------------
    call ops_decl_stencil( 3, 20, a3d_p320_m120_mixed_xy, s3d_p320_m120_mixed_xy, "3,2,0 to -1,-2,0")
    call ops_decl_stencil( 3, 20, a3d_p120_m320_mixed_xy, s3d_p120_m320_mixed_xy, "1,2,0 to -3,-2,0")

    call ops_decl_stencil( 3, 20, a3d_p420_m020_mixed_xy, s3d_p420_m020_mixed_xy, "4,2,0 to 0,-2,0")
    call ops_decl_stencil( 3, 20, a3d_p020_m420_mixed_xy, s3d_p020_m420_mixed_xy, "0,2,0 to -4,-2,0")

    call ops_decl_stencil( 3, 8, a3d_p220_m220_mixed_xy, s3d_p220_m220_mixed_xy, "2,2,0 to -2,-2,0")
    call ops_decl_stencil( 3, 12, a3d_p330_m330_mixed_xy, s3d_p330_m330_mixed_xy, "3,3,0 to -3,-3,0")
    call ops_decl_stencil( 3, 16, a3d_p440_m440_mixed_xy, s3d_p440_m440_mixed_xy, "4,4,0 to -4,-4,0")

    call ops_decl_stencil( 3, 11, a3d_p550_to_m550_xy, s3d_p550_to_m550_xy, "5,5,0 to -5,-5,0")
!-------------------------------------------------------------------------------------------------
    call ops_partition(" ")

END SUBROUTINE ops_data_init

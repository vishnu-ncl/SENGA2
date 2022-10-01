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
    INTEGER a3d_000_to_p400_x(15) /0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0/
    INTEGER a3d_000_to_m400_x(15) /0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0/

    INTEGER a3d_000_to_p500_x(18) /0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0/
    INTEGER a3d_000_to_m500_x(18) /0,0,0, -1,0,0, -2,0,0, -3,0,0, -4,0,0, -5,0,0/

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

    d_size = (/nxsize, nysize, nzsize/)
    d_m = (/-nhalox,-nhaloy,-nhaloz/)
    d_p = (/nhalox,nhaloy,nhaloz/)
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, urhs, d_urhs, "real(dp)", "URHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, drhs, d_drhs, "real(dp)", "DRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, erhs, d_erhs, "real(dp)", "ERHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wrhs, d_wrhs, "real(dp)", "WRHS")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vrhs, d_vrhs, "real(dp)", "VRHS")

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, utmp, d_utmp, "real(dp)", "UTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vtmp, d_vtmp, "real(dp)", "VTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wtmp, d_wtmp, "real(dp)", "WTMP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, trun, d_trun, "real(dp)", "TRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, transp, d_transp, "real(dp)", "TRANSP")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, store7, d_store7, "real(dp)", "STORE7")

!   Declare OPS Stencil
    call ops_decl_stencil( 3, 1, a3d_000, s3d_000, "0,0,0")

    call ops_decl_stencil( 3, 5, a3d_000_to_p400_x, s3d_000_to_p400_x, "0,0,0 to 4,0,0")
    call ops_decl_stencil( 3, 5, a3d_000_to_m400_x, s3d_000_to_m400_x, "0,0,0 to -4,0,0")
    
    call ops_decl_stencil( 3, 6, a3d_000_to_p500_x, s3d_000_to_p500_x, "0,0,0 to 5,0,0")
    call ops_decl_stencil( 3, 6, a3d_000_to_m500_x, s3d_000_to_m500_x, "0,0,0 to -5,0,0")

    call ops_decl_stencil( 3,  5, a3d_p200_to_m200_x, s3d_p200_to_m200_x, "2,0,0 to -2,0,0")
    call ops_decl_stencil( 3,  7, a3d_p300_to_m300_x, s3d_p300_to_m300_x, "3,0,0 to -3,0,0")
    call ops_decl_stencil( 3,  9, a3d_p400_to_m400_x, s3d_p400_to_m400_x, "4,0,0 to -4,0,0")
    call ops_decl_stencil( 3, 11, a3d_p500_to_m500_x, s3d_p500_to_m500_x, "5,0,0 to -5,0,0")

    call ops_decl_stencil( 3,  5, a3d_p300_to_m100_x, s3d_p300_to_m100_x, "3,0,0 to -1,0,0")
    call ops_decl_stencil( 3,  5, a3d_p100_to_m300_x, s3d_p100_to_m300_x, "1,0,0 to -3,0,0")

    call ops_decl_stencil( 3,  6, a3d_p400_to_m100_x, s3d_p400_to_m100_x, "4,0,0 to -1,0,0")
    call ops_decl_stencil( 3,  6, a3d_p100_to_m400_x, s3d_p100_to_m400_x, "1,0,0 to -4,0,0")
!-------------------------------------------------------------------------------------
    call ops_decl_stencil( 3, 5, a3d_000_to_p040_y, s3d_000_to_p040_y, "0,0,0 to 0,4,0")
    call ops_decl_stencil( 3, 5, a3d_000_to_m040_y, s3d_000_to_m040_y, "0,0,0 to  0,-4,0")

    call ops_decl_stencil( 3, 6, a3d_000_to_p050_y, s3d_000_to_p050_y, "0,0,0 to 0,5,0")
    call ops_decl_stencil( 3, 6, a3d_000_to_m050_y, s3d_000_to_m050_y, "0,0,0 to  0,-5,0")

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

MODULE com_ops_senga

    use OPS_Fortran_Reference
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    implicit none

!------------------------------------------------------------------------------------------------------------
!   Senga2 OPS vars

!------------------------------------------------------------------------------------------------------------
!   OPS block
    TYPE(ops_block) :: senga_grid

!------------------------------------------------------------------------------------------------------------
!   OPS dats
    TYPE(ops_dat) :: d_store1, d_store2, d_store3, d_store4, d_store5, d_store6, d_divm

    TYPE(ops_dat) :: d_ucor, d_vcor, d_wcor

    TYPE(ops_dat) :: d_wd1x, d_pd1x, d_td1x
    TYPE(ops_dat) :: d_wd1y, d_pd1y, d_td1y
    TYPE(ops_dat) :: d_wd1z, d_pd1z, d_td1z

    TYPE(ops_dat) :: d_wd2x, d_pd2x, d_td2x
    TYPE(ops_dat) :: d_wd2y, d_pd2y, d_td2y
    TYPE(ops_dat) :: d_wd2z, d_pd2z, d_td2z

    TYPE(ops_dat) :: d_ufxl, d_vfxl, d_wfxl

    TYPE(ops_dat) :: d_drun, d_urun, d_vrun, d_wrun, d_erun
    TYPE(ops_dat) :: d_drun_dump, d_urun_dump, d_vrun_dump, d_wrun_dump, d_erun_dump
    TYPE(ops_dat) :: d_trun_dump
    TYPE(ops_dat) :: d_drhs, d_urhs, d_vrhs, d_wrhs, d_erhs
    TYPE(ops_dat) :: d_derr, d_uerr, d_verr, d_werr, d_eerr
    TYPE(ops_dat) :: d2prun, d2trun

    TYPE(ops_dat) :: d_wmomix, d_difmix, d_tdrmix
    TYPE(ops_dat) :: d_combo1, d_combo2, d_combo3

    TYPE(ops_dat) :: d_utmp, d_vtmp, d_wtmp, d_prun, d_trun, d_transp, d_store7

    TYPE(ops_dat) :: d_itndex(nintmx)
    TYPE(ops_dat) :: d_yrhs_mdim
    TYPE(ops_dat) :: d_yrhs(nspcmx), d_yrun(nspcmx), d_yerr(nspcmx), d_rate(nspcmx), d_rrte(nspcmx)
    TYPE(ops_dat) :: d_yrun_dump(nspcmx)
    TYPE(ops_dat) :: d_ctrans(nspcmx)
    TYPE(ops_dat) :: d_tcoeff, d_tderiv;

    TYPE(ops_dat) :: d_bclyxl(nspcmx), d_stryxl(nspcmx), d_dydtxl(nspcmx), d_ratexl(nspcmx), d_strhxl(nspcmx)
    TYPE(ops_dat) :: d_t6bxl(nspcmx), d_tt6xl(nspcmx)
    TYPE(ops_dat) :: d_yinf1(nspcmx), d_yinf2(nspcmx), d_stlyxl(nspcmx)
    TYPE(ops_dat) :: d_bclyxr(nspcmx), d_stryxr(nspcmx), d_dydtxr(nspcmx), d_ratexr(nspcmx), d_strhxr(nspcmx)
    TYPE(ops_dat) :: d_t6bxr(nspcmx), d_tt6xr(nspcmx)
    TYPE(ops_dat) :: d_bclyyl(nspcmx), d_stryyl(nspcmx), d_dydtyl(nspcmx), d_rateyl(nspcmx), d_strhyl(nspcmx)
    TYPE(ops_dat) :: d_t6byl(nspcmx), d_tt6yl(nspcmx)
    TYPE(ops_dat) :: d_bclyyr(nspcmx), d_stryyr(nspcmx), d_dydtyr(nspcmx), d_rateyr(nspcmx), d_strhyr(nspcmx)
    TYPE(ops_dat) :: d_t6byr(nspcmx), d_tt6yr(nspcmx)
    TYPE(ops_dat) :: d_bclyzl(nspcmx), d_stryzl(nspcmx), d_dydtzl(nspcmx), d_ratezl(nspcmx), d_strhzl(nspcmx)
    TYPE(ops_dat) :: d_t6bzl(nspcmx), d_tt6zl(nspcmx)
    TYPE(ops_dat) :: d_bclyzr(nspcmx), d_stryzr(nspcmx), d_dydtzr(nspcmx), d_ratezr(nspcmx), d_strhzr(nspcmx)
    TYPE(ops_dat) :: d_t6bzr(nspcmx), d_tt6zr(nspcmx)

    TYPE(ops_dat) :: d_totyxl

    TYPE(ops_dat) :: d_bcl1xl, d_bcl2xl, d_bcl3xl, d_bcl4xl, d_bcl5xl, d_bcltxl
    TYPE(ops_dat) :: d_bcl1xr, d_bcl2xr, d_bcl3xr, d_bcl4xr, d_bcl5xr, d_bcltxr

    TYPE(ops_dat) :: d_struxl, d_strvxl, d_strwxl, d_strpxl, d_strdxl, d_strtxl
    TYPE(ops_dat) :: d_strexl, d_strgxl, d_strrxl
    TYPE(ops_dat) :: d_struxr, d_strvxr, d_strwxr, d_strpxr, d_strdxr, d_strtxr
    TYPE(ops_dat) :: d_strexr, d_strgxr, d_strrxr

    TYPE(ops_dat) :: d_dudtxl, d_dvdtxl, d_dwdtxl, d_dtdtxl, d_dddtxl
    TYPE(ops_dat) :: d_dudtxr, d_dvdtxr, d_dwdtxr, d_dtdtxr, d_dddtxr

    TYPE(ops_dat) :: d_acouxl, d_ova2xl, d_gam1xl, d_ovgmxl, d_sydtxl, d_sorpxl
    TYPE(ops_dat) :: d_acouxr, d_ova2xr, d_gam1xr, d_ovgmxr, d_sydtxr, d_sorpxr

    TYPE(ops_dat) :: d_t1bxl, d_t2bxl, d_t3bxl, d_t4bxl, d_t51bxl, d_t52bxl
    TYPE(ops_dat) :: d_t1bxr, d_t2bxr, d_t3bxr, d_t4bxr, d_t51bxr, d_t52bxr

    TYPE(ops_dat) :: d_uinf1, d_uinf2, d_vinf1, d_vinf2, d_winf1, d_winf2, d_ustead

    TYPE(ops_dat) :: d_tt1xl, d_tt2xl, d_tt3xl, d_tt4xl, d_tt5xl
    TYPE(ops_dat) :: d_tt1xr, d_tt2xr, d_tt3xr, d_tt4xr, d_tt5xr

    TYPE(ops_dat) :: d_stluxl, d_stlvxl, d_stlwxl, d_stltxl

    TYPE(ops_dat) :: d_bcl1yl, d_bcl2yl, d_bcl3yl, d_bcl4yl, d_bcl5yl, d_bcltyl
    TYPE(ops_dat) :: d_bcl1yr, d_bcl2yr, d_bcl3yr, d_bcl4yr, d_bcl5yr, d_bcltyr

    TYPE(ops_dat) :: d_struyl, d_strvyl, d_strwyl, d_strpyl, d_strdyl, d_strtyl
    TYPE(ops_dat) :: d_streyl, d_strgyl, d_strryl
    TYPE(ops_dat) :: d_struyr, d_strvyr, d_strwyr, d_strpyr, d_strdyr, d_strtyr
    TYPE(ops_dat) :: d_streyr, d_strgyr, d_strryr

    TYPE(ops_dat) :: d_dudtyl, d_dvdtyl, d_dwdtyl, d_dtdtyl, d_dddtyl
    TYPE(ops_dat) :: d_dudtyr, d_dvdtyr, d_dwdtyr, d_dtdtyr, d_dddtyr

    TYPE(ops_dat) :: d_acouyl, d_ova2yl, d_gam1yl, d_ovgmyl, d_sydtyl, d_sorpyl
    TYPE(ops_dat) :: d_acouyr, d_ova2yr, d_gam1yr, d_ovgmyr, d_sydtyr, d_sorpyr

    TYPE(ops_dat) :: d_t1byl, d_t2byl, d_t3byl, d_t4byl, d_t51byl, d_t52byl
    TYPE(ops_dat) :: d_t1byr, d_t2byr, d_t3byr, d_t4byr, d_t51byr, d_t52byr

    TYPE(ops_dat) :: d_tt1yl, d_tt2yl, d_tt3yl, d_tt4yl, d_tt5yl
    TYPE(ops_dat) :: d_tt1yr, d_tt2yr, d_tt3yr, d_tt4yr, d_tt5yr

    TYPE(ops_dat) :: d_bcl1zl, d_bcl2zl, d_bcl3zl, d_bcl4zl, d_bcl5zl, d_bcltzl
    TYPE(ops_dat) :: d_bcl1zr, d_bcl2zr, d_bcl3zr, d_bcl4zr, d_bcl5zr, d_bcltzr

    TYPE(ops_dat) :: d_struzl, d_strvzl, d_strwzl, d_strpzl, d_strdzl, d_strtzl
    TYPE(ops_dat) :: d_strezl, d_strgzl, d_strrzl
    TYPE(ops_dat) :: d_struzr, d_strvzr, d_strwzr, d_strpzr, d_strdzr, d_strtzr
    TYPE(ops_dat) :: d_strezr, d_strgzr, d_strrzr

    TYPE(ops_dat) :: d_dudtzl, d_dvdtzl, d_dwdtzl, d_dtdtzl, d_dddtzl
    TYPE(ops_dat) :: d_dudtzr, d_dvdtzr, d_dwdtzr, d_dtdtzr, d_dddtzr

    TYPE(ops_dat) :: d_acouzl, d_ova2zl, d_gam1zl, d_ovgmzl, d_sydtzl, d_sorpzl
    TYPE(ops_dat) :: d_acouzr, d_ova2zr, d_gam1zr, d_ovgmzr, d_sydtzr, d_sorpzr

    TYPE(ops_dat) :: d_t1bzl, d_t2bzl, d_t3bzl, d_t4bzl, d_t51bzl, d_t52bzl
    TYPE(ops_dat) :: d_t1bzr, d_t2bzr, d_t3bzr, d_t4bzr, d_t51bzr, d_t52bzr

    TYPE(ops_dat) :: d_tt1zl, d_tt2zl, d_tt3zl, d_tt4zl, d_tt5zl
    TYPE(ops_dat) :: d_tt1zr, d_tt2zr, d_tt3zr, d_tt4zr, d_tt5zr

    TYPE(ops_dat) :: d_crin

!------------------------------------------------------------------------------------------------------------

    TYPE(ops_dat) :: d_ufilt, d_vfilt, d_wfilt, d_yfilt(nspcmx)
    TYPE(ops_dat) :: d_ufold, d_vfold, d_wfold, d_yfold(nspcmx)
    TYPE(ops_dat) :: d_urand, d_vrand, d_wrand, d_yrand(nspcmx)

!------------------------------------------------------------------------------------------------------------
!   OPS reduction handles
    TYPE(ops_reduction) :: h_erdtot, h_erutot, h_ervtot, h_erwtot, h_eretot, h_erytot
    TYPE(ops_reduction) :: h_tket, h_tkes, h_enstro, h_ubart, h_vbart, h_wbart, h_uvart, h_vvart, h_wvart
    TYPE(ops_reduction) :: h_umean, h_denom

!------------------------------------------------------------------------------------------------------------
!   OPS stencils
    TYPE(ops_stencil) :: s3d_000

    TYPE(ops_stencil) :: s3d_000_strid3d_x, s3d_000_strid3d_y, s3d_000_strid3d_z
    TYPE(ops_stencil) :: s3d_000_strid3d_xy, s3d_000_strid3d_xz, s3d_000_strid3d_yz

    TYPE(ops_stencil) :: s3d_000_to_p400_x, s3d_000_to_m400_x
    TYPE(ops_stencil) :: s3d_p100_to_p400_x, s3d_m100_to_m400_x
    TYPE(ops_stencil) :: s3d_p500_to_m500_x

    TYPE(ops_stencil) :: s3d_000_to_p040_y, s3d_000_to_m040_y
    TYPE(ops_stencil) :: s3d_p010_to_p040_y, s3d_m010_to_m040_y
    TYPE(ops_stencil) :: s3d_p050_to_m050_y

    TYPE(ops_stencil) :: s3d_000_to_p004_z, s3d_000_to_m004_z
    TYPE(ops_stencil) :: s3d_p001_to_p004_z, s3d_m001_to_m004_z
    TYPE(ops_stencil) :: s3d_p005_to_m005_z

!------------------------------------------------------------------------------------------------------------

    TYPE(ops_stencil) :: s3d_p000_m440_mixed_xy, s3d_p010_m430_mixed_xy, s3d_p020_m420_mixed_xy, s3d_p030_m410_mixed_xy, s3d_p040_m400_mixed_xy
    TYPE(ops_stencil) :: s3d_p100_m300_mixed_xy, s3d_p100_m340_mixed_xy, s3d_p110_m330_mixed_xy, s3d_p120_m320_mixed_xy, s3d_p140_m300_mixed_xy
    TYPE(ops_stencil) :: s3d_p200_m240_mixed_xy, s3d_p210_m230_mixed_xy, s3d_p220_m220_mixed_xy, s3d_p230_m210_mixed_xy, s3d_p240_m200_mixed_xy
    TYPE(ops_stencil) :: s3d_p300_m100_mixed_xy, s3d_p300_m140_mixed_xy, s3d_p320_m120_mixed_xy, s3d_p330_m110_mixed_xy, s3d_p330_m330_mixed_small_xy, s3d_p330_m330_mixed_xy, s3d_p340_m100_mixed_xy
    TYPE(ops_stencil) :: s3d_p400_p040_mixed_xy, s3d_p410_p030_mixed_xy, s3d_p420_m020_mixed_xy, s3d_p430_m010_mixed_xy, s3d_p440_p000_mixed_xy, s3d_p440_m440_mixed_small_xy, s3d_p440_m440_mixed_xy
    TYPE(ops_stencil) :: s3d_p550_m550_mixed_xy

!------------------------------------------------------------------------------------------------------------

    TYPE(ops_stencil) :: s3d_p000_m404_mixed_xz, s3d_p001_m403_mixed_xz, s3d_p002_m402_mixed_xz, s3d_p003_m401_mixed_xz, s3d_p004_m400_mixed_xz
    TYPE(ops_stencil) :: s3d_p100_m300_mixed_xz, s3d_p100_m304_mixed_xz, s3d_p101_m303_mixed_xz, s3d_p102_m302_mixed_xz, s3d_p104_m300_mixed_xz
    TYPE(ops_stencil) :: s3d_p200_m204_mixed_xz, s3d_p201_m203_mixed_xz, s3d_p202_m202_mixed_xz, s3d_p203_m201_mixed_xz, s3d_p204_m200_mixed_xz
    TYPE(ops_stencil) :: s3d_p300_m100_mixed_xz, s3d_p300_m104_mixed_xz, s3d_p302_m102_mixed_xz, s3d_p303_m101_mixed_xz, s3d_p303_m303_mixed_small_xz, s3d_p303_m303_mixed_xz, s3d_p304_m100_mixed_xz
    TYPE(ops_stencil) :: s3d_p400_p004_mixed_xz, s3d_p401_p003_mixed_xz, s3d_p402_m002_mixed_xz, s3d_p403_m001_mixed_xz, s3d_p404_p000_mixed_xz, s3d_p404_m404_mixed_small_xz, s3d_p404_m404_mixed_xz
    TYPE(ops_stencil) :: s3d_p505_m505_mixed_xz

!------------------------------------------------------------------------------------------------------------

    TYPE(ops_stencil) :: s3d_p000_m044_mixed_yz, s3d_p001_m043_mixed_yz, s3d_p002_m042_mixed_yz, s3d_p003_m041_mixed_yz, s3d_p004_m040_mixed_yz
    TYPE(ops_stencil) :: s3d_p010_m030_mixed_yz, s3d_p010_m034_mixed_yz, s3d_p011_m033_mixed_yz, s3d_p012_m032_mixed_yz, s3d_p014_m030_mixed_yz
    TYPE(ops_stencil) :: s3d_p020_m024_mixed_yz, s3d_p021_m023_mixed_yz, s3d_p022_m022_mixed_yz, s3d_p023_m021_mixed_yz, s3d_p024_m020_mixed_yz
    TYPE(ops_stencil) :: s3d_p030_m010_mixed_yz, s3d_p030_m014_mixed_yz, s3d_p032_m012_mixed_yz, s3d_p033_m011_mixed_yz, s3d_p033_m033_mixed_small_yz, s3d_p033_m033_mixed_yz, s3d_p034_m010_mixed_yz
    TYPE(ops_stencil) :: s3d_p040_p004_mixed_yz, s3d_p041_p003_mixed_yz, s3d_p042_m002_mixed_yz, s3d_p043_m001_mixed_yz, s3d_p044_p000_mixed_yz, s3d_p044_m044_mixed_small_yz, s3d_p044_m044_mixed_yz
    TYPE(ops_stencil) :: s3d_p055_m055_mixed_yz

!------------------------------------------------------------------------------------------------------------
!   ops_halos
    TYPE(ops_halo) :: halos_x(10+2*nspcmx), halos_y(10+2*nspcmx), halos_z(10+2*nspcmx)

!------------------------------------------------------------------------------------------------------------
!   ops_halo group
    TYPE(ops_halo_group) :: halos_grp_x, halos_grp_y, halos_grp_z

END MODULE com_ops_senga

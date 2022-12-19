MODULE com_ops_senga

    use OPS_Fortran_Reference
    
    implicit none

    !Senga2 OPS vars

    TYPE(ops_block) :: senga_grid

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
    TYPE(ops_dat) :: d_drhs, d_urhs, d_vrhs, d_wrhs, d_erhs
    TYPE(ops_dat) :: d_derr, d_uerr, d_verr, d_werr, d_eerr 
    
    TYPE(ops_dat) :: d_wmomix, d_difmix, d_tdrmix
    
    TYPE(ops_dat) :: d_utmp, d_vtmp, d_wtmp, d_prun, d_trun, d_transp, d_store7

    TYPE(ops_dat) :: d_itndex
    TYPE(ops_dat) :: d_yrhs, d_yrun, d_yerr, d_rate, d_rrte

    TYPE(ops_dat) :: d_bclyxl, d_stryxl, d_dydtxl, d_ratexl, d_strhxl
    TYPE(ops_dat) :: d_bclyxr, d_stryxr, d_dydtxr, d_ratexr, d_strhxr
    TYPE(ops_dat) :: d_bclyyl, d_stryyl, d_dydtyl, d_rateyl, d_strhyl
    TYPE(ops_dat) :: d_bclyyr, d_stryyr, d_dydtyr, d_rateyr, d_strhyr
    TYPE(ops_dat) :: d_bclyzl, d_stryzl, d_dydtzl, d_ratezl, d_strhzl
    TYPE(ops_dat) :: d_bclyzr, d_stryzr, d_dydtzr, d_ratezr, d_strhzr

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

    TYPE(ops_dat) :: d_crin

    TYPE(ops_reduction) :: h_erdtot, h_erutot, h_ervtot, h_erwtot, h_eretot, h_erytot, h_prefer
    TYPE(ops_reduction) :: h_tket, h_ubart, h_vbart, h_wbart, h_uvart, h_vvart, h_wvart

    TYPE(ops_stencil) :: s3d_000

    TYPE(ops_stencil) :: s3d_000_strid3d_x, s3d_000_strid3d_y, s3d_000_strid3d_z
    TYPE(ops_stencil) :: s3d_000_strid3d_xy, s3d_000_strid3d_xz, s3d_000_strid3d_yz

    TYPE(ops_stencil) :: s3d_000_to_p400_x, s3d_000_to_m400_x
    TYPE(ops_stencil) :: s3d_000_to_p500_x, s3d_000_to_m500_x
    TYPE(ops_stencil) :: s3d_p100_to_p400_x, s3d_m100_to_m400_x
    TYPE(ops_stencil) :: s3d_p200_to_m200_x, s3d_p300_to_m300_x, s3d_p400_to_m400_x, s3d_p500_to_m500_x
    TYPE(ops_stencil) :: s3d_p300_to_m100_x, s3d_p100_to_m300_x
    TYPE(ops_stencil) :: s3d_p400_to_m100_x, s3d_p100_to_m400_x

    TYPE(ops_stencil) :: s3d_000_to_p040_y, s3d_000_to_m040_y
    TYPE(ops_stencil) :: s3d_000_to_p050_y, s3d_000_to_m050_y
    TYPE(ops_stencil) :: s3d_p010_to_p040_y, s3d_m010_to_m040_y
    TYPE(ops_stencil) :: s3d_p020_to_m020_y, s3d_p030_to_m030_y, s3d_p040_to_m040_y, s3d_p050_to_m050_y
    TYPE(ops_stencil) :: s3d_p030_to_m010_y, s3d_p010_to_m030_y
    TYPE(ops_stencil) :: s3d_p040_to_m010_y, s3d_p010_to_m040_y

    TYPE(ops_stencil) :: s3d_000_to_p004_z, s3d_000_to_m004_z
    TYPE(ops_stencil) :: s3d_000_to_p005_z, s3d_000_to_m005_z
    TYPE(ops_stencil) :: s3d_p001_to_p004_z, s3d_m001_to_m004_z
    TYPE(ops_stencil) :: s3d_p002_to_m002_z, s3d_p003_to_m003_z, s3d_p004_to_m004_z, s3d_p005_to_m005_z
    TYPE(ops_stencil) :: s3d_p003_to_m001_z, s3d_p001_to_m003_z
    TYPE(ops_stencil) :: s3d_p004_to_m001_z, s3d_p001_to_m004_z

    TYPE(ops_stencil) :: s3d_p320_m120_mixed_xy, s3d_p120_m320_mixed_xy
    TYPE(ops_stencil) :: s3d_p420_m020_mixed_xy, s3d_p020_m420_mixed_xy
    TYPE(ops_stencil) :: s3d_p220_m220_mixed_xy, s3d_p330_m330_mixed_xy, s3d_p440_m440_mixed_xy
    TYPE(ops_stencil) :: s3d_p550_to_m550_xy

END MODULE com_ops_senga

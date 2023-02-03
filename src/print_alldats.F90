SUBROUTINE print_alldats()
    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

    CHARACTER (LEN=60) :: fname
    CHARACTER (LEN=3) :: pnxres
    PARAMETER(pnxres = '.h5')
    CHARACTER(len=4) :: citime

    WRITE(citime,'(I4.4)') itime

    fname = 'test_dir/store1_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store1, trim(fname))

    fname = 'test_dir/store2_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store2, trim(fname))

    fname = 'test_dir/store3_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store3, trim(fname))

    fname = 'test_dir/store4_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store4, trim(fname))

    fname = 'test_dir/store5_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store5, trim(fname))

    fname = 'test_dir/store6_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store6, trim(fname))
    
    fname = 'test_dir/divm_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_divm, trim(fname))

    fname = 'test_dir/ucor_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ucor, trim(fname))

    fname = 'test_dir/vcor_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vcor, trim(fname))

    fname = 'test_dir/wcor_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wcor, trim(fname))

    fname = 'test_dir/wd1x_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd1x, trim(fname))
    
    fname = 'test_dir/pd1x_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd1x, trim(fname))

    fname = 'test_dir/td1x_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td1x, trim(fname))

    fname = 'test_dir/wd1y_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd1y, trim(fname))

    fname = 'test_dir/pd1y_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd1y, trim(fname))

    fname = 'test_dir/td1y_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td1y, trim(fname))

    fname = 'test_dir/wd1z_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd1z, trim(fname))

    fname = 'test_dir/pd1z_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd1z, trim(fname))

    fname = 'test_dir/td1z_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td1z, trim(fname))

    fname = 'test_dir/wd2x_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd2x, trim(fname))

    fname = 'test_dir/pd2x_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd2x, trim(fname))

    fname = 'test_dir/td2x_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td2x, trim(fname))

    fname = 'test_dir/wd2y_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd2y, trim(fname))

    fname = 'test_dir/pd2y_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd2y, trim(fname))

    fname = 'test_dir/td2y_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td2y, trim(fname))

    fname = 'test_dir/wd2z_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd2z, trim(fname))

    fname = 'test_dir/pd2z_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd2z, trim(fname))

    fname = 'test_dir/td2z_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td2z, trim(fname))

    fname = 'test_dir/ufxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ufxl, trim(fname))

    fname = 'test_dir/vfxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vfxl, trim(fname))

    fname = 'test_dir/wfxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wfxl, trim(fname))

    fname = 'test_dir/drun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_drun, trim(fname))

    fname = 'test_dir/urun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urun, trim(fname))

    fname = 'test_dir/vrun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vrun, trim(fname))

    fname = 'test_dir/wrun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wrun, trim(fname))

    fname = 'test_dir/erun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_erun, trim(fname))

    fname = 'test_dir/derr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_derr, trim(fname))

    fname = 'test_dir/uerr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_uerr, trim(fname))

    fname = 'test_dir/verr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_verr, trim(fname))

    fname = 'test_dir/werr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_werr, trim(fname))

    fname = 'test_dir/eerr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_eerr, trim(fname))

    fname = 'test_dir/yrun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_yrun, trim(fname))

    fname = 'test_dir/yerr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_yerr, trim(fname))

    fname = 'test_dir/rate_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_rate, trim(fname))

    fname = 'test_dir/rrte_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_rrte, trim(fname))

!    fname = 'test_dir/itndex_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_itndex, trim(fname))

    fname = 'test_dir/yrhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_yrhs, trim(fname))

    fname = 'test_dir/bclyxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bclyxl, trim(fname))

    fname = 'test_dir/bclyxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bclyxr, trim(fname))

    fname = 'test_dir/stryxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_stryxl, trim(fname))

    fname = 'test_dir/stryxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_stryxr, trim(fname))

    fname = 'test_dir/dydtxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dydtxl, trim(fname))

    fname = 'test_dir/dydtxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dydtxr, trim(fname))

    fname = 'test_dir/ratexl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ratexl, trim(fname))

    fname = 'test_dir/ratexr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ratexr, trim(fname))

    fname = 'test_dir/strhxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strhxl, trim(fname))

    fname = 'test_dir/strhxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strhxr, trim(fname))

    fname = 'test_dir/bclyyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bclyyl, trim(fname))

    fname = 'test_dir/bclyyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bclyyr, trim(fname))

    fname = 'test_dir/stryyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_stryyl, trim(fname))

    fname = 'test_dir/stryyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_stryyr, trim(fname))

    fname = 'test_dir/dydtyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dydtyl, trim(fname))

    fname = 'test_dir/dydtyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dydtyr, trim(fname))

    fname = 'test_dir/rateyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_rateyl, trim(fname))

    fname = 'test_dir/rateyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_rateyr, trim(fname))

    fname = 'test_dir/strhyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strhyl, trim(fname))

    fname = 'test_dir/strhyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strhyr, trim(fname))

    fname = 'test_dir/bclyzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bclyzl, trim(fname))

    fname = 'test_dir/bclyzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bclyzr, trim(fname))

    fname = 'test_dir/stryzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_stryzl, trim(fname))

    fname = 'test_dir/stryzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_stryzr, trim(fname))

    fname = 'test_dir/dydtzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dydtzl, trim(fname))

    fname = 'test_dir/dydtzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dydtzr, trim(fname))

    fname = 'test_dir/ratezl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ratezl, trim(fname))

    fname = 'test_dir/ratezr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ratezr, trim(fname))

    fname = 'test_dir/strhzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strhzl, trim(fname))

    fname = 'test_dir/strhzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strhzr, trim(fname))

    fname = 'test_dir/drhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_drhs, trim(fname))

    fname = 'test_dir/urhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urhs, trim(fname))
    
    fname = 'test_dir/vrhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vrhs, trim(fname))

    fname = 'test_dir/wrhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wrhs, trim(fname))

    fname = 'test_dir/erhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_erhs, trim(fname))

    fname = 'test_dir/utmp_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_utmp, trim(fname))

    fname = 'test_dir/vtmp_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vtmp, trim(fname))

    fname = 'test_dir/wtmp_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wtmp, trim(fname))

    fname = 'test_dir/trun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_trun, trim(fname))

    fname = 'test_dir/prun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_prun, trim(fname))

    fname = 'test_dir/transp_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_transp, trim(fname))

    fname = 'test_dir/store7_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store7, trim(fname))    

    fname = 'test_dir/wmomix_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wmomix, trim(fname))

    fname = 'test_dir/difmix_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_difmix, trim(fname))

    fname = 'test_dir/tdrmix_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_tdrmix, trim(fname))

    fname = 'test_dir/bcl1xl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl1xl, trim(fname))

    fname = 'test_dir/bcl1xr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl1xr, trim(fname))

    fname = 'test_dir/bcl2xl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl2xl, trim(fname))

    fname = 'test_dir/bcl2xr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl2xr, trim(fname))

    fname = 'test_dir/bcl3xl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl3xl, trim(fname))

    fname = 'test_dir/bcl3xr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl3xr, trim(fname))

    fname = 'test_dir/bcl4xl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl4xl, trim(fname))

    fname = 'test_dir/bcl4xr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl4xr, trim(fname))

    fname = 'test_dir/bcl5xl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl5xl, trim(fname))

    fname = 'test_dir/bcl5xr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl5xr, trim(fname))

    fname = 'test_dir/bcltxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcltxl, trim(fname))

    fname = 'test_dir/bcltxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcltxr, trim(fname))

    fname = 'test_dir/struxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_struxl, trim(fname))

    fname = 'test_dir/struxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_struxr, trim(fname))

    fname = 'test_dir/strvxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strvxl, trim(fname))

    fname = 'test_dir/strvxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strvxr, trim(fname))

    fname = 'test_dir/strwxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strwxl, trim(fname))

    fname = 'test_dir/strwxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strwxr, trim(fname))

    fname = 'test_dir/strpxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strpxl, trim(fname))

    fname = 'test_dir/strpxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strpxr, trim(fname))

    fname = 'test_dir/strdxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strdxl, trim(fname))

    fname = 'test_dir/strdxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strdxr, trim(fname))

    fname = 'test_dir/strtxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strtxl, trim(fname))

    fname = 'test_dir/strtxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strtxr, trim(fname))

    fname = 'test_dir/strexl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strexl, trim(fname))

    fname = 'test_dir/strexr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strexr, trim(fname))

    fname = 'test_dir/strgxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strgxl, trim(fname))

    fname = 'test_dir/strgxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strgxr, trim(fname))

    fname = 'test_dir/strrxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strrxl, trim(fname))

    fname = 'test_dir/strrxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strrxr, trim(fname))

    fname = 'test_dir/dudtxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dudtxl, trim(fname))

    fname = 'test_dir/dudtxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dudtxr, trim(fname))

    fname = 'test_dir/dvdtxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dvdtxl, trim(fname))

    fname = 'test_dir/dvdtxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dvdtxr, trim(fname))

    fname = 'test_dir/dwdtxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dwdtxl, trim(fname))

    fname = 'test_dir/dwdtxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dwdtxr, trim(fname))

    fname = 'test_dir/dtdtxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dtdtxl, trim(fname))

    fname = 'test_dir/dtdtxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dtdtxr, trim(fname))

    fname = 'test_dir/dddtxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dddtxl, trim(fname))

    fname = 'test_dir/dddtxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dddtxr, trim(fname))

    fname = 'test_dir/acouxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_acouxl, trim(fname))

    fname = 'test_dir/acouxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_acouxr, trim(fname))

    fname = 'test_dir/ova2xl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ova2xl, trim(fname))

    fname = 'test_dir/ova2xr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ova2xr, trim(fname))

    fname = 'test_dir/gam1xl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_gam1xl, trim(fname))

    fname = 'test_dir/gam1xr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_gam1xr, trim(fname))

    fname = 'test_dir/ovgmxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ovgmxl, trim(fname))

    fname = 'test_dir/ovgmxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ovgmxr, trim(fname))

    fname = 'test_dir/sydtxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sydtxl, trim(fname))

    fname = 'test_dir/sydtxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sydtxr, trim(fname))    

    fname = 'test_dir/sorpxl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sorpxl, trim(fname))

    fname = 'test_dir/sorpxr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sorpxr, trim(fname))
    
    fname = 'test_dir/bcl1yl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl1yl, trim(fname))

    fname = 'test_dir/bcl1yr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl1yr, trim(fname))

    fname = 'test_dir/bcl2yl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl2yl, trim(fname))

    fname = 'test_dir/bcl2yr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl2yr, trim(fname))

    fname = 'test_dir/bcl3yl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl3yl, trim(fname))

    fname = 'test_dir/bcl3yr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl3yr, trim(fname))

    fname = 'test_dir/bcl4yl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl4yl, trim(fname))

    fname = 'test_dir/bcl4yr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl4yr, trim(fname))

    fname = 'test_dir/bcl5yl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl5yl, trim(fname))

    fname = 'test_dir/bcl5yr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl5yr, trim(fname))

    fname = 'test_dir/bcltyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcltyl, trim(fname))

    fname = 'test_dir/bcltyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcltyr, trim(fname))

    fname = 'test_dir/struyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_struyl, trim(fname))

    fname = 'test_dir/struyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_struyr, trim(fname))

    fname = 'test_dir/strvyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strvyl, trim(fname))

    fname = 'test_dir/strvyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strvyr, trim(fname))

    fname = 'test_dir/strwyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strwyl, trim(fname))

    fname = 'test_dir/strwyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strwyr, trim(fname))

    fname = 'test_dir/strpyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strpyl, trim(fname))

    fname = 'test_dir/strpyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strpyr, trim(fname))

    fname = 'test_dir/strdyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strdyl, trim(fname))

    fname = 'test_dir/strdyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strdyr, trim(fname))

    fname = 'test_dir/strtyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strtyl, trim(fname))

    fname = 'test_dir/strtyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strtyr, trim(fname))

    fname = 'test_dir/streyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_streyl, trim(fname))

    fname = 'test_dir/streyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_streyr, trim(fname))

    fname = 'test_dir/strgyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strgyl, trim(fname))

    fname = 'test_dir/strgyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strgyr, trim(fname))

    fname = 'test_dir/strryl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strryl, trim(fname))

    fname = 'test_dir/strryr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strryr, trim(fname))

    fname = 'test_dir/dudtyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dudtyl, trim(fname))

    fname = 'test_dir/dudtyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dudtyr, trim(fname))

    fname = 'test_dir/dvdtyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dvdtyl, trim(fname))

    fname = 'test_dir/dvdtyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dvdtyr, trim(fname))

    fname = 'test_dir/dwdtyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dwdtyl, trim(fname))

    fname = 'test_dir/dwdtyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dwdtyr, trim(fname))

    fname = 'test_dir/dtdtyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dtdtyl, trim(fname))

    fname = 'test_dir/dtdtyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dtdtyr, trim(fname))

    fname = 'test_dir/dddtyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dddtyl, trim(fname))

    fname = 'test_dir/dddtyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dddtyr, trim(fname))

    fname = 'test_dir/acouyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_acouyl, trim(fname))

    fname = 'test_dir/acouyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_acouyr, trim(fname))

    fname = 'test_dir/ova2yl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ova2yl, trim(fname))

    fname = 'test_dir/ova2yr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ova2yr, trim(fname))

    fname = 'test_dir/gam1yl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_gam1yl, trim(fname))

    fname = 'test_dir/gam1yr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_gam1yr, trim(fname))

    fname = 'test_dir/ovgmyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ovgmyl, trim(fname))

    fname = 'test_dir/ovgmyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ovgmyr, trim(fname))

    fname = 'test_dir/sydtyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sydtyl, trim(fname))

    fname = 'test_dir/sydtyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sydtyr, trim(fname))

    fname = 'test_dir/sorpyl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sorpyl, trim(fname))

    fname = 'test_dir/sorpyr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sorpyr, trim(fname))

    fname = 'test_dir/bcl1zl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl1zl, trim(fname))

    fname = 'test_dir/bcl1zr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl1zr, trim(fname))

    fname = 'test_dir/bcl2zl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl2zl, trim(fname))

    fname = 'test_dir/bcl2zr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl2zr, trim(fname))

    fname = 'test_dir/bcl3zl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl3zl, trim(fname))

    fname = 'test_dir/bcl3zr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl3zr, trim(fname))

    fname = 'test_dir/bcl4zl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl4zl, trim(fname))

    fname = 'test_dir/bcl4zr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl4zr, trim(fname))

    fname = 'test_dir/bcl5zl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl5zl, trim(fname))

    fname = 'test_dir/bcl5zr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcl5zr, trim(fname))

    fname = 'test_dir/bcltzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcltzl, trim(fname))

    fname = 'test_dir/bcltzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_bcltzr, trim(fname))

    fname = 'test_dir/struzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_struzl, trim(fname))

    fname = 'test_dir/struzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_struzr, trim(fname))

    fname = 'test_dir/strvzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strvzl, trim(fname))

    fname = 'test_dir/strvzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strvzr, trim(fname))

    fname = 'test_dir/strwzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strwzl, trim(fname))

    fname = 'test_dir/strwzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strwzr, trim(fname))

    fname = 'test_dir/strpzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strpzl, trim(fname))

    fname = 'test_dir/strpzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strpzr, trim(fname))

    fname = 'test_dir/strdzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strdzl, trim(fname))

    fname = 'test_dir/strdzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strdzr, trim(fname))

    fname = 'test_dir/strtzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strtzl, trim(fname))

    fname = 'test_dir/strtzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strtzr, trim(fname))

    fname = 'test_dir/strezl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strezl, trim(fname))

    fname = 'test_dir/strezr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strezr, trim(fname))

    fname = 'test_dir/strgzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strgzl, trim(fname))

    fname = 'test_dir/strgzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strgzr, trim(fname))

    fname = 'test_dir/strrzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strrzl, trim(fname))

    fname = 'test_dir/strrzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_strrzr, trim(fname))

    fname = 'test_dir/dudtzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dudtzl, trim(fname))

    fname = 'test_dir/dudtzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dudtzr, trim(fname))

    fname = 'test_dir/dvdtzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dvdtzl, trim(fname))

    fname = 'test_dir/dvdtzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dvdtzr, trim(fname))

    fname = 'test_dir/dwdtzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dwdtzl, trim(fname))

    fname = 'test_dir/dwdtzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dwdtzr, trim(fname))

    fname = 'test_dir/dtdtzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dtdtzl, trim(fname))

    fname = 'test_dir/dtdtzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dtdtzr, trim(fname))

    fname = 'test_dir/dddtzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dddtzl, trim(fname))

    fname = 'test_dir/dddtzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_dddtzr, trim(fname))

    fname = 'test_dir/acouzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_acouzl, trim(fname))

    fname = 'test_dir/acouzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_acouzr, trim(fname))

    fname = 'test_dir/ova2zl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ova2zl, trim(fname))

    fname = 'test_dir/ova2zr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ova2zr, trim(fname))

    fname = 'test_dir/gam1zl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_gam1zl, trim(fname))

    fname = 'test_dir/gam1zr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_gam1zr, trim(fname))

    fname = 'test_dir/ovgmzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ovgmzl, trim(fname))

    fname = 'test_dir/ovgmzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ovgmzr, trim(fname))

    fname = 'test_dir/sydtzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sydtzl, trim(fname))

    fname = 'test_dir/sydtzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sydtzr, trim(fname))

    fname = 'test_dir/sorpzl_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sorpzl, trim(fname))

    fname = 'test_dir/sorpzr_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_sorpzr, trim(fname))

    fname = 'test_dir/crin_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_crin, trim(fname))

END SUBROUTINE print_alldats

SUBROUTINE print_alldats()
    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

    integer(kind=4) :: ispec,iindex,dtime
    character(len=60) :: fname
    character(len=3) :: pnxhdf
    parameter(pnxhdf = '.h5')
    character(len=8) :: citime

    !dtime=itime
    dtime=INT(itime/ntdump)
    WRITE(citime,'(I8.8)') dtime

!   store1, store2, store3, store4, store5, store6, divm
    fname = 'test_dir/dats1_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store1, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store2, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store3, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store4, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store5, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store6, trim(fname))
    call ops_fetch_dat_hdf5_file(d_divm, trim(fname))

!   ucor, vcor, wcor
    fname = 'test_dir/dats2_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ucor, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vcor, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wcor, trim(fname))

!   wd1x, wd1y, wd1z, wd2x, wd2y, wd2z, pd1x, pd1y, pd1z
!   pd2x, pd2y, pd2z, td1x, td1y, td1z td2x, td2y, td2z
    fname = 'test_dir/dats3_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd1x, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd1y, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd1z, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd2x, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd2y, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wd2z, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd1x, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd1y, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd1z, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd2x, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd2y, trim(fname))
    call ops_fetch_dat_hdf5_file(d_pd2z, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td1x, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td1y, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td1z, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td2x, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td2y, trim(fname))
    call ops_fetch_dat_hdf5_file(d_td2z, trim(fname))

!   ufxl, vfxl, wfxl
    fname = 'test_dir/dats4_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_ufxl, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vfxl, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wfxl, trim(fname))

!   drun, urun, vrun, wrun, erun, yrun
    fname = 'test_dir/dats5_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_drun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vrun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wrun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_erun, trim(fname))
    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_yrun(ispec), trim(fname))
    END DO

!   derr, uerr, verr, werr, eerr, yerr, rate, rrte
    fname = 'test_dir/dats6_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_derr, trim(fname))
    call ops_fetch_dat_hdf5_file(d_uerr, trim(fname))
    call ops_fetch_dat_hdf5_file(d_verr, trim(fname))
    call ops_fetch_dat_hdf5_file(d_werr, trim(fname))
    call ops_fetch_dat_hdf5_file(d_eerr, trim(fname))
    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_yerr(ispec), trim(fname))
    END DO
    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_rate(ispec), trim(fname))
    END DO
    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_rrte(ispec), trim(fname))
    END DO

!   drhs, urhs, vrhs
    fname = 'test_dir/dats7_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_drhs, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urhs, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vrhs, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wrhs, trim(fname))
    call ops_fetch_dat_hdf5_file(d_erhs, trim(fname))
    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_yrhs(ispec), trim(fname))
    END DO
    DO iindex = 1,nintmx
        call ops_fetch_dat_hdf5_file(d_itndex(iindex), trim(fname))
    END DO

!   utmp, vtmp, wtmp, trun, prun, transp, store7
    fname = 'test_dir/dats8_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_utmp, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vtmp, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wtmp, trim(fname))
    call ops_fetch_dat_hdf5_file(d_trun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_prun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_transp, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store7, trim(fname))

!   wmomix, difmix, tdrmix
    fname = 'test_dir/dats9_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wmomix, trim(fname))
    call ops_fetch_dat_hdf5_file(d_difmix, trim(fname))
    call ops_fetch_dat_hdf5_file(d_tdrmix, trim(fname))

!   x-boundary array
!    fname = 'test_dir/dats10_timestep'//citime//pnxhdf
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl1xl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl1xr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl2xl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl2xr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl3xl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl3xr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl4xl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl4xr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl5xl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl5xr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcltxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcltxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_struxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_struxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strvxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strvxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strwxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strwxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strpxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strpxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strdxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strdxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strtxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strtxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strexl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strexr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strgxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strgxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strrxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strrxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dudtxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dudtxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dvdtxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dvdtxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dwdtxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dwdtxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dtdtxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dtdtxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dddtxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_acouxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_acouxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ova2xl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ova2xr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_gam1xl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_gam1xr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ovgmxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ovgmxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sydtxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sydtxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sorpxl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sorpxr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_tt1xl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_tt1xr, trim(fname))
!    DO ispec = 1,nspec
!        call ops_fetch_dat_hdf5_file(d_strhxl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_strhxr(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_ratexl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_ratexr(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_stryxl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_stryxr(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_bclyxl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_bclyxr(ispec), trim(fname))
!    END DO
!
!!   y-boundary array
!    fname = 'test_dir/dats11_timestep'//citime//pnxhdf
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl1yl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl1yr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl2yl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl2yr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl3yl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl3yr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl4yl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl4yr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl5yl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl5yr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcltyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcltyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_struyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_struyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strvyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strvyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strwyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strwyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strpyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strpyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strdyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strdyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strtyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strtyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_streyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_streyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strgyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strgyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strryl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strryr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dudtyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dudtyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dvdtyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dvdtyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dwdtyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dwdtyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dtdtyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dtdtyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dddtyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dddtyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_acouyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_acouyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ova2yl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ova2yr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_gam1yl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_gam1yr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ovgmyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ovgmyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sydtyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sydtyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sorpyl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sorpyr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_tt1yl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_tt1yr, trim(fname))
!    DO ispec = 1,nspec
!        call ops_fetch_dat_hdf5_file(d_strhyl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_strhyr(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_rateyl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_rateyr(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_stryyl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_stryyr(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_bclyyl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_bclyyr(ispec), trim(fname))
!    END DO
!
!!   z-boundary array
!    fname = 'test_dir/dats12_timestep'//citime//pnxhdf
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl1zl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl1zr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl2zl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl2zr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl3zl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl3zr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl4zl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl4zr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl5zl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcl5zr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcltzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_bcltzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_struzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_struzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strvzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strvzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strwzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strwzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strpzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strpzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strdzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strdzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strtzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strtzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strezl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strezr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strgzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strgzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strrzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_strrzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dudtzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dudtzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dvdtzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dvdtzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dwdtzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dwdtzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dtdtzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dtdtzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dddtzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_dddtzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_acouzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_acouzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ova2zl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ova2zr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_gam1zl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_gam1zr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ovgmzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_ovgmzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sydtzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sydtzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sorpzl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_sorpzr, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_tt1zl, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_tt1zr, trim(fname))
!    DO ispec = 1,nspec
!        call ops_fetch_dat_hdf5_file(d_strhzl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_strhzr(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_ratezl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_ratezr(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_stryzl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_stryzr(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_bclyzl(ispec), trim(fname))
!        call ops_fetch_dat_hdf5_file(d_bclyzr(ispec), trim(fname))
!    END DO
!
!    fname = 'test_dir/crin_timestep'//citime//pnxhdf
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_crin, trim(fname))

END SUBROUTINE print_alldats

SUBROUTINE print_dats()
    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

    integer(kind=4) :: ispec
    character(len=60) :: fname
    character(len=3) :: pnxhdf
    parameter(pnxhdf = '.h5')
    character(len=8) :: citime
    character(len=4) :: pnxtxt
    parameter(pnxtxt = '.txt')


    WRITE(citime,'(I8.8)') itime

    fname = 'test_dir/drhs_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_drhs, trim(fname))

    fname = 'test_dir/urhs_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urhs, trim(fname))

    fname = 'test_dir/vrhs_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vrhs, trim(fname))

    fname = 'test_dir/wrhs_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wrhs, trim(fname))

    fname = 'test_dir/erhs_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_erhs, trim(fname))

    fname = 'test_dir/yrhs_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_yrhs(ispec), trim(fname))
    END DO

    fname = 'test_dir/drun_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_drun, trim(fname))

!    fname = 'test_dir/drun_timestep'//citime//pnxtxt
!    call ops_print_dat_to_txtfile(d_drun, trim(fname))

    fname = 'test_dir/urun_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urun, trim(fname))

!    fname = 'test_dir/urun_timestep'//citime//pnxtxt
!    call ops_print_dat_to_txtfile(d_urun, trim(fname))

    fname = 'test_dir/vrun_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vrun, trim(fname))

!    fname = 'test_dir/vrun_timestep'//citime//pnxtxt
!    call ops_print_dat_to_txtfile(d_vrun, trim(fname))

    fname = 'test_dir/wrun_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wrun, trim(fname))

!    fname = 'test_dir/wrun_timestep'//citime//pnxtxt
!    call ops_print_dat_to_txtfile(d_wrun, trim(fname))

    fname = 'test_dir/erun_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_erun, trim(fname))

!    fname = 'test_dir/erun_timestep'//citime//pnxtxt
!    call ops_print_dat_to_txtfile(d_erun, trim(fname))

    fname = 'test_dir/yrun_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_yrun(ispec), trim(fname))
    END DO

    fname = 'test_dir/rate_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_rate(ispec), trim(fname))
    END DO

    fname = 'test_dir/rrte_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_rrte(ispec), trim(fname))
    END DO

    fname = 'test_dir/trun_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_trun, trim(fname))

    fname = 'test_dir/prun_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_prun, trim(fname))

    fname = 'test_dir/store1_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store1, trim(fname))

    fname = 'test_dir/store2_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store2, trim(fname))

    fname = 'test_dir/store3_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store3, trim(fname))

    fname = 'test_dir/store4_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store4, trim(fname))

    fname = 'test_dir/store5_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store5, trim(fname))

    fname = 'test_dir/store6_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store6, trim(fname))

    fname = 'test_dir/store7_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_store7, trim(fname))

    fname = 'test_dir/utmp_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_utmp, trim(fname))

    fname = 'test_dir/vtmp_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vtmp, trim(fname))

    fname = 'test_dir/wtmp_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wtmp, trim(fname))

    fname = 'test_dir/transp_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_transp, trim(fname))

    fname = 'test_dir/divm_timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_divm, trim(fname))

END SUBROUTINE print_dats

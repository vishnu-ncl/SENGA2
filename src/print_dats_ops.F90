
! Auto-generated at 2026-04-28 18:43:09.657702 by ops-translator
SUBROUTINE print_dats
  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE OPS_Fortran_hdf5_Declarations
  USE OPS_CONSTANTS

  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  INTEGER(KIND = 4) :: ispec
  CHARACTER(LEN = 60) :: fname
  CHARACTER(LEN = 3) :: pnxhdf
  PARAMETER(pnxhdf = '.h5')
  CHARACTER(LEN = 8) :: citime
  CHARACTER(LEN = 4) :: pnxtxt
  PARAMETER(pnxtxt = '.txt')


  WRITE(citime, '(I8.8)') itime

  fname = 'test_dir/drhs_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_drhs, TRIM(fname))

  fname = 'test_dir/urhs_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_urhs, TRIM(fname))

  fname = 'test_dir/vrhs_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_vrhs, TRIM(fname))

  fname = 'test_dir/wrhs_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_wrhs, TRIM(fname))

  fname = 'test_dir/erhs_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_erhs, TRIM(fname))

  fname = 'test_dir/yrhs_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  DO ispec = 1, nspcmx
    CALL ops_fetch_dat_hdf5_file(d_yrhs(ispec), TRIM(fname))
  END DO

  fname = 'test_dir/drun_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_drun, TRIM(fname))

  !    fname = 'test_dir/drun_timestep'//citime//pnxtxt
  !    call ops_print_dat_to_txtfile(d_drun, trim(fname))

  fname = 'test_dir/urun_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_urun, TRIM(fname))

  !    fname = 'test_dir/urun_timestep'//citime//pnxtxt
  !    call ops_print_dat_to_txtfile(d_urun, trim(fname))

  fname = 'test_dir/vrun_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_vrun, TRIM(fname))

  !    fname = 'test_dir/vrun_timestep'//citime//pnxtxt
  !    call ops_print_dat_to_txtfile(d_vrun, trim(fname))

  fname = 'test_dir/wrun_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_wrun, TRIM(fname))

  !    fname = 'test_dir/wrun_timestep'//citime//pnxtxt
  !    call ops_print_dat_to_txtfile(d_wrun, trim(fname))

  fname = 'test_dir/erun_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_erun, TRIM(fname))

  !    fname = 'test_dir/erun_timestep'//citime//pnxtxt
  !    call ops_print_dat_to_txtfile(d_erun, trim(fname))

  fname = 'test_dir/yrun_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  DO ispec = 1, nspcmx
    CALL ops_fetch_dat_hdf5_file(d_yrun(ispec), TRIM(fname))
  END DO

  fname = 'test_dir/rate_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  DO ispec = 1, nspcmx
    CALL ops_fetch_dat_hdf5_file(d_rate(ispec), TRIM(fname))
  END DO

  fname = 'test_dir/rrte_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  DO ispec = 1, nspcmx
    CALL ops_fetch_dat_hdf5_file(d_rrte(ispec), TRIM(fname))
  END DO

  fname = 'test_dir/trun_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_trun, TRIM(fname))

  fname = 'test_dir/prun_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_prun, TRIM(fname))

  fname = 'test_dir/store1_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_store1, TRIM(fname))

  fname = 'test_dir/store2_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_store2, TRIM(fname))

  fname = 'test_dir/store3_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_store3, TRIM(fname))

  fname = 'test_dir/store4_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_store4, TRIM(fname))

  fname = 'test_dir/store5_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_store5, TRIM(fname))

  fname = 'test_dir/store6_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_store6, TRIM(fname))

  fname = 'test_dir/store7_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_store7, TRIM(fname))

  fname = 'test_dir/utmp_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_utmp, TRIM(fname))

  fname = 'test_dir/vtmp_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_vtmp, TRIM(fname))

  fname = 'test_dir/wtmp_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_wtmp, TRIM(fname))

  fname = 'test_dir/transp_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_transp, TRIM(fname))

  fname = 'test_dir/divm_timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_divm, TRIM(fname))

END SUBROUTINE print_dats
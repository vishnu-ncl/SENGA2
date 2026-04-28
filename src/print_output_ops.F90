
! Auto-generated at 2026-04-28 18:43:09.669668 by ops-translator
SUBROUTINE print_output
  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE copy_kernel_module
  USE OPS_Fortran_hdf5_Declarations
  USE OPS_CONSTANTS

  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  INTEGER(KIND = 4) :: ispec, dtime
  CHARACTER(LEN = 60) :: fname
  CHARACTER(LEN = 3) :: pnxhdf
  PARAMETER(pnxhdf = '.h5')
  CHARACTER(LEN = 8) :: citime
  INTEGER(KIND = 4) :: rangexyz(6)

  dtime = INT(itime / ntdump)
  WRITE(citime, '(I8.8)') dtime

  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d2prun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ))

  CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d2trun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ))

  fname = 'output/timestep' // citime // pnxhdf
  CALL ops_fetch_block_hdf5_file(senga_grid, TRIM(fname))

  CALL ops_fetch_dat_hdf5_file(d_drun, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_urun, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_vrun, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_wrun, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d_erun, TRIM(fname))

    DO ispec = 1, nspcmx
    CALL ops_fetch_dat_hdf5_file(d_yrun(ispec), TRIM(fname))
    CALL ops_fetch_dat_hdf5_file(d_rrte(ispec), TRIM(fname))
  END DO

  CALL ops_fetch_dat_hdf5_file(d2trun, TRIM(fname))
  CALL ops_fetch_dat_hdf5_file(d2prun, TRIM(fname))

END SUBROUTINE print_output
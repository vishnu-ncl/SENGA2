SUBROUTINE print_dats()
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

    fname = 'drhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_drhs, trim(fname))

    fname = 'urhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urhs, trim(fname))
    
    fname = 'urhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urhs, trim(fname))

    fname = 'wrhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wrhs, trim(fname))

    fname = 'erhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_erhs, trim(fname))

    fname = 'yrhs_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_yrhs, trim(fname))

    fname = 'drun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_drun, trim(fname))

    fname = 'urun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urun, trim(fname))

    fname = 'vrun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vrun, trim(fname))

    fname = 'wrun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wrun, trim(fname))

    fname = 'erun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_erun, trim(fname))

    fname = 'yrun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_yrun, trim(fname))

    fname = 'trun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_trun, trim(fname))

    fname = 'prun_timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_prun, trim(fname))

END SUBROUTINE print_dats

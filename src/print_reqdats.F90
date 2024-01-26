SUBROUTINE print_reqdats()
    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

    integer(kind=4) :: ispec,dtime
    CHARACTER (LEN=60) :: fname
    CHARACTER (LEN=3) :: pnxres
    PARAMETER(pnxres = '.h5')
    CHARACTER(len=4) :: citime
    INTEGER :: rangexyz(6)
    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)

    dtime=int(itime/ntdump)

    WRITE(citime,'(I4.4)') dtime

!    fname = 'test_dir/drhs_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_drhs, trim(fname))
!
!    fname = 'test_dir/urhs_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_urhs, trim(fname))
!
!    fname = 'test_dir/vrhs_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_vrhs, trim(fname))
!
!    fname = 'test_dir/wrhs_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_wrhs, trim(fname))
!
!    fname = 'test_dir/erhs_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_erhs, trim(fname))
!
!    fname = 'test_dir/yrhs_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    DO ispec = 1,nspcmx
!        call ops_fetch_dat_hdf5_file(d_yrhs(ispec), trim(fname))
!    END DO


    fname = 'test_dir/timestep'//citime//pnxres
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
    call ops_fetch_dat_hdf5_file(d_drun, trim(fname))

    call ops_fetch_dat_hdf5_file(d_urun, trim(fname))

    call ops_fetch_dat_hdf5_file(d_vrun, trim(fname))

    call ops_fetch_dat_hdf5_file(d_wrun, trim(fname))

!    fname = 'test_dir/erun_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_erun, trim(fname))

    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_yrun(ispec), trim(fname))
        call ops_fetch_dat_hdf5_file(d_rate(ispec), trim(fname))
    END DO

    call ops_par_loop(copy_kernel,"copy",senga_grid,3,rangexyz, &
            ops_arg_dat(d2prun, 1,s3d_000,"real(kind=8)",OPS_WRITE), &
            ops_arg_dat(d_prun, 1,s3d_000,"real(kind=8)",OPS_READ))

    call ops_par_loop(copy_kernel,"copy",senga_grid,3,rangexyz, &
            ops_arg_dat(d2trun, 1,s3d_000,"real(kind=8)",OPS_WRITE), &
            ops_arg_dat(d_trun, 1,s3d_000,"real(kind=8)",OPS_READ))


    call ops_fetch_dat_hdf5_file(d2trun, trim(fname))

    call ops_fetch_dat_hdf5_file(d2prun, trim(fname))

!!VM!!    fname = 'test_dir/drun_timestep'//citime//pnxres
!!VM!!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!!VM!!    call ops_fetch_dat_hdf5_file(d_drun, trim(fname))
!!VM!!
!!VM!!    fname = 'test_dir/urun_timestep'//citime//pnxres
!!VM!!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!!VM!!    call ops_fetch_dat_hdf5_file(d_urun, trim(fname))
!!VM!!
!!VM!!    fname = 'test_dir/vrun_timestep'//citime//pnxres
!!VM!!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!!VM!!    call ops_fetch_dat_hdf5_file(d_vrun, trim(fname))
!!VM!!
!!VM!!    fname = 'test_dir/wrun_timestep'//citime//pnxres
!!VM!!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!!VM!!    call ops_fetch_dat_hdf5_file(d_wrun, trim(fname))
!!VM!!
!!VM!!!    fname = 'test_dir/erun_timestep'//citime//pnxres
!!VM!!!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!!VM!!!    call ops_fetch_dat_hdf5_file(d_erun, trim(fname))
!!VM!!
!!VM!!    fname = 'test_dir/yrun_timestep'//citime//pnxres
!!VM!!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!!VM!!    DO ispec = 1,nspcmx
!!VM!!        call ops_fetch_dat_hdf5_file(d_yrun(ispec), trim(fname))
!!VM!!    END DO
!!VM!!
!!VM!!    call ops_par_loop(copy_kernel,"copy",senga_grid,3,rangexyz, &
!!VM!!            ops_arg_dat(d2prun, 1,s3d_000,"real(kind=8)",OPS_WRITE), &
!!VM!!            ops_arg_dat(d_prun, 1,s3d_000,"real(kind=8)",OPS_READ))
!!VM!!
!!VM!!    call ops_par_loop(copy_kernel,"copy",senga_grid,3,rangexyz, &
!!VM!!            ops_arg_dat(d2trun, 1,s3d_000,"real(kind=8)",OPS_WRITE), &
!!VM!!            ops_arg_dat(d_trun, 1,s3d_000,"real(kind=8)",OPS_READ))
!!VM!!
!!VM!!
!!VM!!    fname = 'test_dir/trun_timestep'//citime//pnxres
!!VM!!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!!VM!!    call ops_fetch_dat_hdf5_file(d2trun, trim(fname))
!!VM!!
!!VM!!    fname = 'test_dir/prun_timestep'//citime//pnxres
!!VM!!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!!VM!!    call ops_fetch_dat_hdf5_file(d2prun, trim(fname))

!    fname = 'test_dir/store1_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_store1, trim(fname))
!
!    fname = 'test_dir/store2_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_store2, trim(fname))
!
!    fname = 'test_dir/store3_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_store3, trim(fname))
!
!    fname = 'test_dir/store4_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_store4, trim(fname))
!
!    fname = 'test_dir/store5_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_store5, trim(fname))
!
!    fname = 'test_dir/store6_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_store6, trim(fname))
!
!    fname = 'test_dir/store7_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_store7, trim(fname))
!
!    fname = 'test_dir/utmp_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_utmp, trim(fname))
!
!    fname = 'test_dir/vtmp_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_vtmp, trim(fname))
!
!    fname = 'test_dir/wtmp_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_wtmp, trim(fname))
!
!    fname = 'test_dir/transp_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_transp, trim(fname))
!
!    fname = 'test_dir/divm_timestep'//citime//pnxres
!    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
!    call ops_fetch_dat_hdf5_file(d_divm, trim(fname))

END SUBROUTINE print_reqdats

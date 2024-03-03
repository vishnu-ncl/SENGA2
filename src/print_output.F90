SUBROUTINE print_output()
    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations

    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

    integer(kind=4) :: ispec,dtime
    character(len=60) :: fname
    character(len=3) :: pnxhdf
    parameter(pnxhdf = '.h5')
    character(len=8) :: citime
    integer(kind=4) :: ic,jc,kc
    
    dtime=INT(itime/ntdump)
    WRITE(citime,'(I8.8)') dtime

    DO kc = 1,nzsize
        DO jc = 1,nysize
            DO ic = 1,nxsize
                prun2(ic,jc,kc)=prun(ic,jc,kc)
                trun2(ic,jc,kc)=trun(ic,jc,kc)
            END DO
        END DO
    END DO

    fname = 'output/timestep'//citime//pnxhdf
    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))

    call ops_fetch_dat_hdf5_file(d_drun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vrun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wrun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_erun, trim(fname))

    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_yrun(ispec), trim(fname))
        call ops_fetch_dat_hdf5_file(d_rrte(ispec), trim(fname))
    END DO

    call ops_fetch_dat_hdf5_file(d2trun, trim(fname))
    call ops_fetch_dat_hdf5_file(d2prun, trim(fname))

END SUBROUTINE print_output

PROGRAM flamin

    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    implicit none

!   *************************************************************************

!   FLAMIN
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   28-DEC-2003:  CREATED
!   08-JAN-2005:  RSC INITIAL 1D LAMINAR FLAME PROFILE

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   SETS INITIAL THERMOCHEMICAL FIELD
!   1D LAMINAR FLAME PROFILE (LEFT OR RIGHT FACING)
!   SPECIAL FOR 21 STEP HYDROGEN MECHAMISM

!   *************************************************************************

!   LOCAL DATA
!   ==========
    integer(kind=4), parameter :: nxglbl=3000, nyglbl=1, nzglbl=1
    integer(kind=4), parameter :: nspcmx=9
    real(kind=8) :: drun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: urun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: vrun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: wrun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: trun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: yrun(nxglbl,nyglbl,nzglbl,nspcmx)
    integer(kind=4) :: ix, ic, ispec, index1d, fnum

    real(kind=8) :: drun_1D(nxglbl*nyglbl*nzglbl), urun_1D(nxglbl*nyglbl*nzglbl), &
                    vrun_1D(nxglbl*nyglbl*nzglbl), wrun_1D(nxglbl*nyglbl*nzglbl), &
                    trun_1D(nxglbl*nyglbl*nzglbl), &
                    yrun1_1D(nxglbl*nyglbl*nzglbl), yrun2_1D(nxglbl*nyglbl*nzglbl), &
                    yrun3_1D(nxglbl*nyglbl*nzglbl), yrun4_1D(nxglbl*nyglbl*nzglbl), &
                    yrun5_1D(nxglbl*nyglbl*nzglbl), yrun6_1D(nxglbl*nyglbl*nzglbl), &
                    yrun7_1D(nxglbl*nyglbl*nzglbl), yrun8_1D(nxglbl*nyglbl*nzglbl), &
                    yrun9_1D(nxglbl*nyglbl*nzglbl)

    character(len=60) :: fname
    character(len=3) :: pnxhdf
    parameter(pnxhdf = '.h5')

!   Declare ops_block and ops_dats
    TYPE(ops_block) :: senga_grid
    TYPE(ops_dat) :: d_drun, d_urun, d_vrun, d_wrun, d_trun
    TYPE(ops_dat) :: d_yrun(nspcmx)
    TYPE(ops_dat) :: d_temp

!   Declare stencil
    TYPE(ops_stencil) :: s3d_000

    integer(kind=4) :: a3d_000(3) = [0,0,0]
    integer(kind=4) :: d_size(3) = [nxglbl,nyglbl,nzglbl]
    integer(kind=4) :: d_base(3) = [1,1,1]
    integer(kind=4) :: d_p(3) = [0,0,0]
    integer(kind=4) :: d_m(3) = [0,0,0]

    real(kind=8), dimension(:), allocatable :: temp_real_null

    character(len=20) :: buf

    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   Init OPS environment
    call ops_init(2)

    call ops_decl_block(3, senga_grid, "SENGA_GRID")

!   Create OPS dats and refer to 1D arrays so that OPS_DATS will be automatically set
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, drun_1D, d_drun, "real(kind=8)", "DRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, urun_1D, d_urun, "real(kind=8)", "URUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, vrun_1D, d_vrun, "real(kind=8)", "VRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, wrun_1D, d_wrun, "real(kind=8)", "WRUN")
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, trun_1D, d_trun, "real(kind=8)", "TRUN")
    
    ispec = 1
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun1_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 2
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun2_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 3
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun3_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 4
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun4_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 5
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun5_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 6
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun6_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 7
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun7_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 8
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun8_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 9
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun9_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_temp, "real(kind=8)", "temp")

    call ops_decl_stencil( 3, 1, a3d_000, s3d_000, "0,0,0")

    call ops_partition(" ")

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_temp, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!   =========================================================================

    fnum = 111
    OPEN(UNIT=fnum,FILE='h2_init.dat',FORM="FORMATTED",STATUS= "OLD")
    DO ic=1,nxglbl
        READ(fnum,*)ix,urun(ic,1,1),trun(ic,1,1),drun(ic,1,1),  &
            yrun(ic,1,1,1),yrun(ic,1,1,2),yrun(ic,1,1,3),yrun(ic,1,1,4),  &
            yrun(ic,1,1,5),yrun(ic,1,1,6),yrun(ic,1,1,7),yrun(ic,1,1,8), yrun(ic,1,1,9)
    END DO
    CLOSE(fnum)

    DO ic = 1,nxglbl
        vrun(ic,1,1) = 0.0_8
        wrun(ic,1,1) = 0.0_8
    END DO

    index1d = 0
    DO ic = 1, nxglbl
        index1d = index1d + 1

        drun_1D(index1d) = drun(ic,1,1)
        urun_1D(index1d) = urun(ic,1,1)
        vrun_1D(index1d) = vrun(ic,1,1)
        wrun_1D(index1d) = wrun(ic,1,1)
        trun_1D(index1d) = trun(ic,1,1)

        yrun1_1D(index1d) = yrun(ic,1,1,1)
        yrun2_1D(index1d) = yrun(ic,1,1,2)
        yrun3_1D(index1d) = yrun(ic,1,1,3)
        yrun4_1D(index1d) = yrun(ic,1,1,4)
        yrun5_1D(index1d) = yrun(ic,1,1,5)
        yrun6_1D(index1d) = yrun(ic,1,1,6)
        yrun7_1D(index1d) = yrun(ic,1,1,7)
        yrun8_1D(index1d) = yrun(ic,1,1,8)
        yrun9_1D(index1d) = yrun(ic,1,1,9)
        
    END DO

    fname = 'flame_init'//pnxhdf

    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))

    call ops_fetch_dat_hdf5_file(d_drun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_urun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_vrun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_wrun, trim(fname))
    call ops_fetch_dat_hdf5_file(d_trun, trim(fname))

    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_yrun(ispec), trim(fname))
    END DO

!   Finalize OPS env
    call ops_exit( )

!   =========================================================================

END PROGRAM flamin

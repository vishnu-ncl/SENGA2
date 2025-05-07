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
    integer(kind=4), parameter :: nxglbl=780, nyglbl=216, nzglbl=216
    integer(kind=4), parameter :: nspcmx=31

    real(kind=8) :: drun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: urun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: vrun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: wrun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: trun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: yrun(nxglbl,nyglbl,nzglbl,nspcmx)

    real(kind=8) :: urunl(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: vrunl(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: wrunl(nxglbl,nyglbl,nzglbl)

    integer(kind=4) :: ic,jc,kc,index1d
    integer(kind=4) :: ispec

    real(kind=8) :: drun_1D(nxglbl*nyglbl*nzglbl), urun_1D(nxglbl*nyglbl*nzglbl), &
                    vrun_1D(nxglbl*nyglbl*nzglbl), wrun_1D(nxglbl*nyglbl*nzglbl), &
                    trun_1D(nxglbl*nyglbl*nzglbl), &
                    yrun1_1D(nxglbl*nyglbl*nzglbl), yrun2_1D(nxglbl*nyglbl*nzglbl), &
                    yrun3_1D(nxglbl*nyglbl*nzglbl), yrun4_1D(nxglbl*nyglbl*nzglbl), &
                    yrun5_1D(nxglbl*nyglbl*nzglbl), yrun6_1D(nxglbl*nyglbl*nzglbl), &
                    yrun7_1D(nxglbl*nyglbl*nzglbl), yrun8_1D(nxglbl*nyglbl*nzglbl), &
                    yrun9_1D(nxglbl*nyglbl*nzglbl), yrun10_1D(nxglbl*nyglbl*nzglbl), &
                    yrun11_1D(nxglbl*nyglbl*nzglbl), yrun12_1D(nxglbl*nyglbl*nzglbl), &
                    yrun13_1D(nxglbl*nyglbl*nzglbl), yrun14_1D(nxglbl*nyglbl*nzglbl), &
                    yrun15_1D(nxglbl*nyglbl*nzglbl), yrun16_1D(nxglbl*nyglbl*nzglbl), &
                    yrun17_1D(nxglbl*nyglbl*nzglbl), yrun18_1D(nxglbl*nyglbl*nzglbl), &
                    yrun19_1D(nxglbl*nyglbl*nzglbl), yrun20_1D(nxglbl*nyglbl*nzglbl), &
                    yrun21_1D(nxglbl*nyglbl*nzglbl), yrun22_1D(nxglbl*nyglbl*nzglbl), &
                    yrun23_1D(nxglbl*nyglbl*nzglbl), yrun24_1D(nxglbl*nyglbl*nzglbl), &
                    yrun25_1D(nxglbl*nyglbl*nzglbl), yrun26_1D(nxglbl*nyglbl*nzglbl), &
                    yrun27_1D(nxglbl*nyglbl*nzglbl), yrun28_1D(nxglbl*nyglbl*nzglbl), &
                    yrun29_1D(nxglbl*nyglbl*nzglbl), yrun30_1D(nxglbl*nyglbl*nzglbl), &
                    yrun31_1D(nxglbl*nyglbl*nzglbl)

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

!     BEGIN
!     =====

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

    ispec = 10
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun10_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 11
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun11_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 12
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun12_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 13
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun13_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 14
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun14_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 15
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun15_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 16
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun16_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 17
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun17_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 18
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun18_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 19
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun19_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 20
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun20_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 21
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun21_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 22
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun22_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 23
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun23_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 24
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun24_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 25
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun25_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 26
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun26_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 27
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun27_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 28
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun28_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 29
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun29_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 30
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun30_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    ispec = 31
    WRITE(buf,"(A4,I2.2)") "YRUN",ispec
    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, yrun31_1D, d_yrun(ispec), "real(kind=8)", trim(buf))

    call ops_decl_dat(senga_grid, 1, d_size, d_base, d_m, d_p, temp_real_null, d_temp, "real(kind=8)", "temp")

    call ops_decl_stencil( 3, 1, a3d_000, s3d_000, "0,0,0")

    call ops_partition(" ")

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_temp, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!   =========================================================================

    OPEN(UNIT=11,FILE='turbin.dat',STATUS='OLD',FORM='UNFORMATTED')
    READ(11)drun,urun,vrun,wrun
    CLOSE(11)

    OPEN(UNIT=10,FILE='lamsol.dat',STATUS='OLD',FORM='FORMATTED')
    DO kc = 1,nzglbl
        DO jc = 1,nyglbl
            DO ic = 1,nxglbl
                READ(10,*)drun(ic,jc,kc), urunl(ic,jc,kc),  &
                            vrunl(ic,jc,kc), wrunl(ic,jc,kc),  &
                            trun(ic,jc,kc), yrun(ic,jc,kc,1),  &
                            yrun(ic,jc,kc,2), yrun(ic,jc,kc,3),  &
                            yrun(ic,jc,kc,4), yrun(ic,jc,kc,5),  &
                            yrun(ic,jc,kc,6), yrun(ic,jc,kc,7),  &
                            yrun(ic,jc,kc,8), yrun(ic,jc,kc,9),  &
                            yrun(ic,jc,kc,10), yrun(ic,jc,kc,11),  &
                            yrun(ic,jc,kc,12), yrun(ic,jc,kc,13),  &
                            yrun(ic,jc,kc,14), yrun(ic,jc,kc,15),  &
                            yrun(ic,jc,kc,16), yrun(ic,jc,kc,17),  &
                            yrun(ic,jc,kc,18), yrun(ic,jc,kc,19),  &
                            yrun(ic,jc,kc,20), yrun(ic,jc,kc,21),  &
                            yrun(ic,jc,kc,22), yrun(ic,jc,kc,23),  &
                            yrun(ic,jc,kc,24), yrun(ic,jc,kc,25),  &
                            yrun(ic,jc,kc,26), yrun(ic,jc,kc,27),  &
                            yrun(ic,jc,kc,28), yrun(ic,jc,kc,29),  &
                            yrun(ic,jc,kc,30), yrun(ic,jc,kc,31)
            END DO
        END DO
    END DO
    CLOSE(10)

    DO kc = 1,nzglbl
        DO jc = 1,nyglbl
            DO ic = 1,nxglbl
                urun(ic,jc,kc)=urun(ic,jc,kc)+urunl(ic,jc,kc)
                vrun(ic,jc,kc)=vrun(ic,jc,kc)+vrunl(ic,jc,kc)
                wrun(ic,jc,kc)=wrun(ic,jc,kc)+wrunl(ic,jc,kc)
            END DO
        END DO
    END DO

    index1d = 0
    DO ic = 1, nxglbl
        index1d = index1d + 1

        drun_1D(index1d) = drun(ic,jc,kc)
        urun_1D(index1d) = urun(ic,jc,kc)
        vrun_1D(index1d) = vrun(ic,jc,kc)
        wrun_1D(index1d) = wrun(ic,jc,kc)
        trun_1D(index1d) = trun(ic,jc,kc)

        yrun1_1D(index1d) = yrun(ic,jc,kc,1)
        yrun2_1D(index1d) = yrun(ic,jc,kc,2)
        yrun3_1D(index1d) = yrun(ic,jc,kc,3)
        yrun4_1D(index1d) = yrun(ic,jc,kc,4)
        yrun5_1D(index1d) = yrun(ic,jc,kc,5)
        yrun6_1D(index1d) = yrun(ic,jc,kc,6)
        yrun7_1D(index1d) = yrun(ic,jc,kc,7)
        yrun8_1D(index1d) = yrun(ic,jc,kc,8)
        yrun9_1D(index1d) = yrun(ic,jc,kc,9)
        yrun10_1D(index1d) = yrun(ic,jc,kc,10)
        yrun11_1D(index1d) = yrun(ic,jc,kc,11)
        yrun12_1D(index1d) = yrun(ic,jc,kc,12)
        yrun13_1D(index1d) = yrun(ic,jc,kc,13)
        yrun14_1D(index1d) = yrun(ic,jc,kc,14)
        yrun15_1D(index1d) = yrun(ic,jc,kc,15)
        yrun16_1D(index1d) = yrun(ic,jc,kc,16)
        yrun17_1D(index1d) = yrun(ic,jc,kc,17)
        yrun18_1D(index1d) = yrun(ic,jc,kc,18)
        yrun19_1D(index1d) = yrun(ic,jc,kc,19)
        yrun20_1D(index1d) = yrun(ic,jc,kc,20)
        yrun21_1D(index1d) = yrun(ic,jc,kc,21)
        yrun22_1D(index1d) = yrun(ic,jc,kc,22)
        yrun23_1D(index1d) = yrun(ic,jc,kc,23)
        yrun24_1D(index1d) = yrun(ic,jc,kc,24)
        yrun25_1D(index1d) = yrun(ic,jc,kc,25)
        yrun26_1D(index1d) = yrun(ic,jc,kc,26)
        yrun27_1D(index1d) = yrun(ic,jc,kc,27)
        yrun28_1D(index1d) = yrun(ic,jc,kc,28)
        yrun29_1D(index1d) = yrun(ic,jc,kc,29)
        yrun30_1D(index1d) = yrun(ic,jc,kc,30)
        yrun31_1D(index1d) = yrun(ic,jc,kc,31)
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


RETURN

END PROGRAM flamin

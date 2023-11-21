SUBROUTINE bctixl

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCTIXL
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   30-MAR-2006:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   INITIALISES TURBULENT INFLOW

!   X-DIRECTION LEFT-HAND END

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   PARAMETERS
!   ==========
    CHARACTER (LEN=4) :: pntixl,pntcxl
    CHARACTER (LEN=4) :: pnxdat
    PARAMETER(pntixl='tixl', pntcxl='tcxl', pnxdat='.dat')

!   LOCAL DATA
!   ==========
    integer(kind=4) :: icproc
    integer(kind=4) :: kxmodd,ixmodd
    integer(kind=4) :: itotxl
    CHARACTER (LEN=6) :: pnproc
    LOGICAL :: fxdump
    integer(kind=4) :: rangexyz(6)

    real(kind=8) ::  store1(1,1,1), store4(1,1,1), store5(1,1,1), store6(1,1,1)

!   BEGIN
!   =====

!   =========================================================================

!   BUILD THE FILENAMES FOR THE INLET TURBULENT VELOCITY FIELD
    WRITE(pnproc,'(I6.6)')iproc
!   RESTART FILE
    fntixl = pntixl//pnproc//pnxdat
!   COLD START FILE
    fntcxl = pntcxl//pnproc//pnxdat

    nctixl = 3

!   =========================================================================

!   INLET COLD START SWITCH
!   PARAMETER I2=0
    IF(nxlprm(2) == 0) THEN

!       =======================================================================

!       INLET COLD START
!       ----------------

!       CHECK AND INITIALISE RESTART FILE
        INQUIRE(FILE=fntixl,EXIST=fxdump)
        IF(.NOT.fxdump) THEN
            OPEN(UNIT=nctixl,FILE=fntixl,STATUS='NEW',FORM='UNFORMATTED')
            CLOSE(nctixl)
        END IF

!       READ THE COLD START INLET TURBULENT VELOCITY FIELD
        OPEN(UNIT=nctixl,FILE=fntcxl,STATUS='OLD',FORM='UNFORMATTED')
            READ(nctixl)store1,store4,store5,store6
        CLOSE(nctixl)

!       SET THE REAL PARTS
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ))

        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ))

        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ))


!       ZERO THE IMAGINARY PARTS
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!       PARTIAL (X-WISE) FOURIER TRANSFORM
        call buftxl

!       MEAN INLET VELOCITY
!       SCANNING PLANE LOCATION AND VELOCITY
!       PARAMETER R1=MEAN VEL, R2=EXTRA SCAN VEL, R3=INITIAL LOCATION
        bvelxl = rxlprm(1)
        svelxl = bvelxl + rxlprm(2)
        elocxl = xgdlen - rxlprm(3)

!       =======================================================================

    ELSE

!       =======================================================================

!       INLET RESTART
!       -------------
!       READ THE RESTART INLET TURBULENT VELOCITY FIELD
        write(*, '(a)') "Using the arrays not allocated by OPS, &
                        Please implement the function in OPS first, bctixl.F90: ID=137"
            STOP
        OPEN(UNIT=nctixl,FILE=fntixl,STATUS='OLD',FORM='UNFORMATTED')
            READ(nctixl)ufxl,vfxl,wfxl,elocxl,svelxl,bvelxl
        CLOSE(nctixl)

!       =======================================================================

    END IF

!   =========================================================================

!   INITIALISE RUNNING LOCATION
    slocxl = elocxl

!   INITIALISE SCALE FACTORS
    tpovxg = two*pi/xgdlen
    scauxl = two/REAL(nxglbl,kind=8)
    scduxl = -two*pi*scauxl*REAL(nxglbl-1,kind=8)/REAL(nxglbl,kind=8)/xgdlen

!   INITIALISE FLAGS AND INDICES FOR INLET PLANE DFT
    fllixl = .false.
    fltrxl = .false.
    istaxl = 1
    istoxl = npmapx(ixproc)

    kminxl = 0
    DO icproc = 0, ixproc-1
        kminxl = kminxl + npmapx(icproc)
    END DO
    kxmodd = MOD(kminxl,2)
    itotxl = kminxl + npmapx(ixproc)

    ixmodd = MOD(istoxl,2)

    IF(kxmodd == 0) THEN

!       EVEN NUMBER OF POINTS TO LH SIDE OF THIS PROCESSOR
!       LH PROCESSOR HAS TRAILING REAL VALUE
!       LOCAL PROCESSOR MUST HAVE LEADING IMAGINARY VALUE
!       EXCEPT FOR FIRST POINT
        IF(kminxl > 0) fllixl = .true.
        kminxl = kminxl/2
        istaxl = 2
        IF(ixmodd == 0) THEN
!           EVEN NUMBER OF POINTS ON LOCAL PROCESSOR
!           LOCAL PROCESSOR MUST HAVE TRAILING REAL VALUE
            istoxl = istoxl - 1
            IF(itotxl == nxglbl) THEN
!               END OF PENCIL: IGNORE TRAILING REAL VALUE
                fltrxl = .false.
            ELSE
                fltrxl = .true.
            END IF
        ELSE
!           ODD NUMBER OF POINTS ON LOCAL PROCESSOR
!           LOCAL PROCESSOR HAS NO TRAILING REAL VALUE
            fltrxl = .false.
        END IF

    ELSE

!       ODD NUMBER OF POINTS TO LH SIDE OF THIS PROCESSOR
!       LH PROCESSOR HAS NO TRAILING REAL VALUE
!       LOCAL PROCESSOR MUST HAVE NO LEADING IMAGINARY VALUE
        fllixl = .false.
        kminxl = kminxl/2 + 1
        istaxl = 1
        IF(ixmodd == 0) THEN
!           EVEN NUMBER OF POINTS ON LOCAL PROCESSOR
!           LOCAL PROCESSOR HAS NO TRAILING REAL VALUE
            fltrxl = .false.
        ELSE
!           ODD NUMBER OF POINTS ON LOCAL PROCESSOR
!           LOCAL PROCESSOR MUST HAVE TRAILING REAL VALUE
            istoxl = istoxl - 1
            IF(itotxl == nxglbl) THEN
!               END OF PENCIL: IGNORE TRAILING REAL VALUE
                fltrxl = .false.
            ELSE
                fltrxl = .true.
            END IF
        END IF

    END IF

!   =========================================================================

END SUBROUTINE bctixl

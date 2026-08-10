SUBROUTINE output

    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
    use OPS_CONSTANTS

    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   OUTPUT
!   ======

!   AUTHOR
!   ------
!   R.S.CANT - CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   28-DEC-2003:  RSC MODIFIED FOR SENGA2
!   11-JUL-2009:  RSC ADD A DUMP FORMAT SWITCH; REPORT THE DUMP
!   29-AUG-2009:  RSC UPDATE NUMBER OF PROCESSORS

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   PROCESSES THE RESULTS

!   FILES: REPORT FILE       - UNIT NCREPT - NAME FNREPT  - FORMATTED
!          DUMP OUTPUT FILES - UNIT NCDMPO - NAMES FNDMPO - UNFORMATTED
!          STATISTICS FILE   - UNIT NCSTAT - NAME FNSTAT  - FORMATTED

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------

!   LOCAL parameters
!   ================
!   DIAGNOSTICS
    character(len=4) :: pnxres
    parameter(pnxres = '.res')
    character(len=2) :: pnxryf
    parameter(pnxryf = 'yf')
    integer(kind=4) :: ncdiag
    parameter(ncdiag = 11)
    integer(kind=4) :: iddump

!   LOCAL DATA
!   ==========
!   DIAGNOSTICS
    real(kind=8) :: deltag, deltagx, deltagy, deltagz, fornow
!    real(kind=8) :: ttemp(nxsize,nysize,nzsize)
!    real(kind=8) :: ptemp(nxsize,nysize,nzsize)
!    real(kind=8) :: ytemp(nspec,nxsize,nysize,nzsize)
    real(kind=8) :: tkeg,enstrg

    integer(kind=4) :: ispec
    integer(kind=4) :: ic, jc, kc
    integer(kind=4) :: ix
    integer(kind=4) :: rangexyz(6)
    character(len=21) :: fndiag

!   RSC UPDATE NUMBER OF PROCESSORS
    character(len=6) :: pnproc
    character(len=11) :: strqty
    character(len=2) :: strspc
    logical :: bcflag, file_exist

    character(len=60) :: fname
    character(len=4) :: proc
    character(len=5) :: ipdump

    character(len=25) :: fndump
    character(len=16) :: pndump
    character(len=3) :: pnxhdf
    parameter(pndump = 'output/dmpi_dats', pnxhdf = '.h5')

!   BEGIN
!   =====

!   =========================================================================

!   REPORT OUTPUT
!   =============
    IF(MOD(itime,ntrept) == 0) THEN

!       REPORT ON PROCESSOR NO.1 ONLY
!       ------
        IF (ops_is_root() == 1) THEN

            OPEN(UNIT=ncrept,FILE=fnrept,STATUS='OLD',FORM='FORMATTED')

!           GO TO EOF
            1000      CONTINUE
            READ(ncrept,9000,END=1010)
            GO TO 1000
            1010      BACKSPACE(ncrept)

            WRITE(ncrept,9100)itime
            WRITE(ncrept,9110)etime,tstep
            CLOSE(ncrept)

        END IF

!       =======================================================================

!       DIAGNOSTICS
!ops        jc = MAX(nyglbl/2,1)
!ops        kc = MAX(nzglbl/2,1)

!ops        WRITE(pnproc,'(I6.6)')iproc

!       GLOBAL INDEXING
!ops        deltag = xgdlen/(REAL(nxglbl-1))

!        STRQTY = 'output/pres'
!        FNDIAG = STRQTY//PNPROC//PNXRES
!        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C       GO TO EOF
!8000    CONTINUE
!          READ(NCDIAG,9000,END=8001)
!          GOTO 8000
!8001    BACKSPACE(NCDIAG)
!        DO IC = ISTAL,ISTOL
!          IX = IGOFST + IC
!          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,PRUN(IC,JC,KC)
!        ENDDO
!        WRITE(NCDIAG,*)
!        CLOSE(NCDIAG)

!        STRQTY = 'output/uvel'
!        FNDIAG = STRQTY//PNPROC//PNXRES
!        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C       GO TO EOF
!8010    CONTINUE
!          READ(NCDIAG,9000,END=8011)
!          GOTO 8010
!8011    BACKSPACE(NCDIAG)
!        DO IC = ISTAL,ISTOL
!          IX = IGOFST + IC
!          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,
!     +                      URUN(IC,JC,KC)/DRUN(IC,JC,KC)
!        ENDDO
!        WRITE(NCDIAG,*)
!        CLOSE(NCDIAG)

!        STRQTY = 'output/temp'
!        FNDIAG = STRQTY//PNPROC//PNXRES
!        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C       GO TO EOF
!8020    CONTINUE
!          READ(NCDIAG,9000,END=8021)
!          GOTO 8020
!8021    BACKSPACE(NCDIAG)
!        DO IC = ISTAL,ISTOL
!          IX = IGOFST + IC
!          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,TRUN(IC,JC,KC)
!        ENDDO
!        WRITE(NCDIAG,*)
!        CLOSE(NCDIAG)

!        STRQTY = 'output/dens'
!        FNDIAG = STRQTY//PNPROC//PNXRES
!        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C       GO TO EOF
!8030    CONTINUE
!          READ(NCDIAG,9000,END=8031)
!          GOTO 8030
!8031    BACKSPACE(NCDIAG)
!        DO IC = ISTAL,ISTOL
!          IX = IGOFST + IC
!          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,DRUN(IC,JC,KC)
!        ENDDO
!        WRITE(NCDIAG,*)
!        CLOSE(NCDIAG)

!        STRQTY = 'output/ener'
!        FNDIAG = STRQTY//PNPROC//PNXRES
!        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C       GO TO EOF
!8040    CONTINUE
!          READ(NCDIAG,9000,END=8041)
!          GOTO 8040
!8041    BACKSPACE(NCDIAG)
!        DO IC = ISTAL,ISTOL
!          IX = IGOFST + IC
!          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,
!     +                      ERUN(IC,JC,KC)/DRUN(IC,JC,KC)
!        ENDDO
!        WRITE(NCDIAG,*)
!        CLOSE(NCDIAG)

!        DO ISPEC = 1, NSPEC

!          WRITE(STRSPC,'(I2.2)')ISPEC
!          STRQTY = 'output/'//PNXRYF//STRSPC
!          FNDIAG = STRQTY//PNPROC//PNXRES
!          OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!C         GO TO EOF
!8050      CONTINUE
!            READ(NCDIAG,9000,END=8051)
!            GOTO 8050
!8051      BACKSPACE(NCDIAG)
!          DO IC = ISTAL,ISTOL
!            IX = IGOFST + IC
!              FORNOW = YRUN(IC,JC,KC,ISPEC)/DRUN(IC,JC,KC)
!              FORNOW = MAX(FORNOW,1.0D-30)
!            WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,FORNOW
!C     +                        YRUN(IC,JC,KC,ISPEC)/DRUN(IC,JC,KC)
!          ENDDO
!          WRITE(NCDIAG,*)
!          CLOSE(NCDIAG)

!        ENDDO

!C        ISPEC = NSPEC
!C        WRITE(STRSPC,'(I2.2)')ISPEC
!C        STRQTY = PNXRYF//STRSPC
!C        FNDIAG = STRQTY//PNPROC//PNXRES
!C        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
!CC       GO TO EOF
!C8050    CONTINUE
!C          READ(NCDIAG,9000,END=8051)
!C          GOTO 8050
!C8051    BACKSPACE(NCDIAG)
!C        DO IC = ISTAL,ISTOL
!C          IX = IGOFST + IC
!C          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,
!C     +                      YRUN(IC,JC,KC,ISPEC)/DRUN(IC,JC,KC)
!C        ENDDO
!C        WRITE(NCDIAG,*)
!C        CLOSE(NCDIAG)

    END IF

!   =========================================================================
!   UMOD START
!   DATA OUTPUT FOR POST-PROCESSING
!   ======
!ops    IF(MOD(itime,ntdump) == 0) THEN
!ops        iddump=itime/ntdump
!ops        WRITE(ipdump,'(I5.5)') iddump
!ops        WRITE(proc,'(I4.4)') iproc

!ops        fname = 'output/out'//ipdump//proc//pnxres


!ops        DO kc = 1,nzsize
!ops            DO jc = 1,nysize
!ops                DO ic = 1,nxsize
!ops                    DO ispec =1,nspec
!ops                        ytemp(ispec,ic,jc,kc)=yrun(ispec,ic,jc,kc)/drun(ic,jc,kc)
!ops                    END DO
!ops                    ttemp(ic,jc,kc)=trun(ic,jc,kc)
!ops                    ptemp(ic,jc,kc)=prun(ic,jc,kc)
!ops                END DO
!ops            END DO
!ops        END DO

!ops        OPEN(UNIT=16,FILE=trim(fname),FORM='UNFORMATTED',STATUS='NEW')
!ops        WRITE(16)drun,urun/drun,vrun/drun,wrun/drun,erun/drun,ttemp,  &
!ops                 ptemp,ytemp,rrte,etime
!ops        CLOSE(16)
!ops    END IF

!   UMOD END
!   =====

!   FULL DUMP OUTPUT
!   ================
!    IF( MOD(itime,ntdump) == 0 .and. (.not. (((itime == ntime1) .or. (itime == 0)) .and. ncdmpi == 1)) ) THEN
    IF( MOD(itime,ntdump) == 0 ) THEN

!       CARRY OUT A FULL DUMP
!       ---------------------
!       Datasets dumped for output visualization as well as for restart purpose
        call print_output()

!       USE THE DUMP FILE INDICATED BY IDFLAG
!       RSC 11-JUL-2009 ADD A DUMP FORMAT SWITCH
        idflag = MOD(INT(itime/ntdump), 2) + 1
        IF( ndofmt == 0 ) THEN

            IF (nxlprm(1) == 4) THEN
                IF(ops_is_root() == 1)THEN
                    OPEN(UNIT=ncdmpo,FILE='output/intran.dat',STATUS='unknown',  &
                         FORM='FORMATTED')
                    WRITE(ncdmpo,*)intran
                    CLOSE(ncdmpo)
                END IF

                fname = 'output/inflow'//pnxhdf
                call ops_fetch_block_hdf5_file(senga_grid, trim(fname))
                call ops_fetch_dat_hdf5_file(d_uinf2, trim(fname))
                call ops_fetch_dat_hdf5_file(d_vinf2, trim(fname))
                call ops_fetch_dat_hdf5_file(d_winf2, trim(fname))

                DO ispec = 1,nspcmx
                    call ops_fetch_dat_hdf5_file(d_yinf2(ispec), trim(fname))
                END DO
            
            END IF

!           UNFORMATTED DUMP OUTPUT
            OPEN(UNIT=ncdmpo, FILE=fndmpo(idflag), STATUS='OLD', FORM='UNFORMATTED')

            IF (ops_is_root() == 1) THEN
                WRITE(*,*) "Writing run information to file(unformatted): ", trim(fndmpo(idflag)), "  idflag: ", idflag
            END IF

            REWIND(ncdmpo)
            WRITE(ncdmpo)nxglbl,nyglbl,nzglbl,nspec,&
                         etime,tstep,errold,errldr
            CLOSE(ncdmpo)
        ELSE
!           FORMATTED DUMP OUTPUT
            OPEN(UNIT=ncdmpo,FILE=fndmpo(idflag),STATUS='OLD', FORM='FORMATTED')

            IF (ops_is_root() == 1) THEN
                WRITE(*,*) "Writing run information to file(formatted): ", trim(fndmpo(idflag)), "  idflag: ", idflag
            END IF

            REWIND(ncdmpo)
            WRITE(ncdmpo,*)nxglbl,nyglbl,nzglbl,nspec
            WRITE(ncdmpo,*)etime,tstep,errold,errldr
            CLOSE(ncdmpo)
        END IF

!       REPORT THE DUMP
!       RSC 11-JUL-2009
        IF (ops_is_root() == 1) THEN

            OPEN(UNIT=ncrept,FILE=fnrept,STATUS='OLD',FORM='FORMATTED')
            3000      CONTINUE
            READ(ncrept,9000,END=3010)
            GO TO 3000
            3010      BACKSPACE(ncrept)
            WRITE(ncrept,9120)fndmpo(idflag)
            CLOSE(ncrept)

        END IF

        IF (ops_is_root() == 1) THEN
            INQUIRE(FILE="output/filed_time.dat",EXIST=file_exist)
            IF ( file_exist ) THEN
                OPEN(UNIT=1011,FILE="output/filed_time.dat",STATUS='OLD',POSITION='APPEND',FORM='FORMATTED')
            ELSE
                OPEN(UNIT=1011,FILE="output/filed_time.dat",STATUS='NEW',FORM='FORMATTED')
            END IF
            WRITE(1011,*) INT(itime/ntdump), etime
            CLOSE(1011)
        END IF

    END IF

!   =========================================================================

!   DUMP BC INFORMATION AS REQUIRED
!   ===============================
!ops    IF(MOD(itime,ntdump) == 0) THEN

!ops        bcflag = (nsbcxl == nsbci2).OR.(nsbcxl == nsbci3)
!ops        bcflag = bcflag.AND.(nxlprm(1) == 3)

!ops        IF(bcflag) THEN

!           DUMP THE INLET TURBULENT VELOCITY FIELD
!ops            OPEN(UNIT=nctixl,FILE=fntixl,STATUS='OLD', FORM='UNFORMATTED')
!ops            REWIND(nctixl)
!ops            WRITE(nctixl)ufxl,vfxl,wfxl,slocxl,svelxl,bvelxl
!ops            CLOSE(nctixl)

!ops        END IF

!ops    END IF

!   =========================================================================

!   TIME STEP HISTORY
    IF (ops_is_root() == 1) THEN
        WRITE(*,'(I7,1PE12.4,I5)')itime,tstep,inderr
    END IF

!   =========================================================================

!   STATISTICS ON THE FLY
!   =====================

!   STATISTICS MASTER SWITCH
!   ------------------------
    IF(ntstat >= 0) THEN

!       =========================================================================

!       OUTPUT STATISTICS
!       =================

!       STATISTICS ON ONE PROCESSOR ONLY
!       ----------
        IF(iproc == 0) THEN

            IF(MOD(itime,ntstat) == 0) THEN

                OPEN(UNIT=ncstat,FILE=fnstat,STATUS='OLD',FORM='FORMATTED')

!               GO TO EOF
                2000        CONTINUE
                READ(ncstat,9200,END=2010)
                GO TO 2000
                2010        BACKSPACE(ncstat)

                WRITE(ncstat,9100)itime

                CLOSE(ncstat)

            END IF

        END IF

!       RESET STORAGE INDEX
        itstat = 0

!       =========================================================================

!   STATISTICS MASTER SWITCH
    END IF

!   =========================================================================

9000  FORMAT(a)
9100  FORMAT('Time step number: ',i7)
9110  FORMAT('Elapsed time: ',1PE12.4,';',2X,'next time step:',1PE12.4)
9120  FORMAT('Dump completed: ',a)
9200  FORMAT(i5,/ 5X,6(1PE12.4)/  &
    5X,6(1PE12.4)/ 5X,6(1PE12.4)/  &
    5X,4(1PE12.4)/ 5X,3(1PE12.4)/  &
    5X,3(1PE12.4)/ 5X,3(1PE12.4)/  &
    5X,3(1PE12.4)/ 5X,3(1PE12.4)/  &
    5X,2(1PE12.4))
9300  FORMAT(2(1PE15.7))

END SUBROUTINE output

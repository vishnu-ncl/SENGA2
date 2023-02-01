SUBROUTINE output

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
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
#ifdef HDF5
    use hdf5io
#endif
!   -------------------------------------------------------------------------

!   LOCAL PARAMETERS
!   ================
!   DIAGNOSTICS
    CHARACTER (LEN=4) :: pnxres
    PARAMETER(pnxres = '.res')
    CHARACTER (LEN=2) :: pnxryf
    PARAMETER(pnxryf = 'yf')
    INTEGER :: ncdiag
    PARAMETER(ncdiag = 11)
    INTEGER :: iddump

!   LOCAL DATA
!   ==========
!   DIAGNOSTICS
    REAL(KIND=8) :: deltag,fornow
    REAL(KIND=8) :: ttemp(nxsize,nysize,nzsize)
    REAL(KIND=8) :: ptemp(nxsize,nysize,nzsize)
    REAL(KIND=8) :: ytemp(nspec,nxsize,nysize,nzsize)

    INTEGER :: ispec
    INTEGER :: ic,jc,kc
    INTEGER :: igofst,ix
    INTEGER :: rangexyz(6)
    CHARACTER (LEN=21) :: fndiag

!   RSC UPDATE NUMBER OF PROCESSORS
    CHARACTER (LEN=6) :: pnproc
    CHARACTER (LEN=11) :: strqty
    CHARACTER (LEN=2) :: strspc
    LOGICAL :: bcflag

    CHARACTER (LEN=60) :: fname
    CHARACTER (LEN=4) :: proc
    CHARACTER (LEN=5) :: ipdump

!   BEGIN
!   =====

!   =========================================================================

!   REPORT OUTPUT
!   =============
    IF(MOD(itime,ntrept) == 0) THEN
  
!       REPORT ON PROCESSOR NO.1 ONLY
!       ------
        IF(iproc == 0) THEN
    
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
        jc = MAX(nyglbl/2,1)
        kc = MAX(nzglbl/2,1)
  
        WRITE(pnproc,'(I6.6)')iproc
  
!       GLOBAL INDEXING
        deltag = xgdlen/(REAL(nxglbl-1))
  
        igofst = 0
        DO ic = 0, ixproc-1
            igofst = igofst + npmapx(ic)
        END DO
  
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
    IF(MOD(itime,ntdump) == 0) THEN
        iddump=itime/ntdump
        WRITE(ipdump,'(I5.5)') iddump
        WRITE(proc,'(I4.4)') iproc

        fname = 'output/out'//ipdump//proc//pnxres


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
    END IF

!   UMOD END
!   =====
!   FULL DUMP OUTPUT
!   ================
    IF(MOD(itime,ntdump) == 0) THEN
  
!       CARRY OUT A FULL DUMP
!       ---------------------
!       USE THE DUMP FILE INDICATED BY IDFLAG
!       RSC 11-JUL-2009 ADD A DUMP FORMAT SWITCH
#ifndef HDF5
        IF(ndofmt == 0) THEN
    
!           UNFORMATTED DUMP OUTPUT
!ops            OPEN(UNIT=ncdmpo,FILE=fndmpo(idflag+1),STATUS='OLD', FORM='UNFORMATTED')
!ops            REWIND(ncdmpo)
!ops            WRITE(ncdmpo)nxsize,nysize,nzsize,nspec, drun,urun,vrun,wrun,erun,yrun,  &
!ops                         etime,tstep,errold,errldr
!ops            CLOSE(ncdmpo)
    
        ELSE
    
!           FORMATTED DUMP OUTPUT
!ops            OPEN(UNIT=ncdmpo,FILE=fndmpo(idflag+1),STATUS='OLD', FORM='FORMATTED')
!ops            REWIND(ncdmpo)
!ops            WRITE(ncdmpo,*)nxsize,nysize,nzsize,nspec
!ops            DO kc = 1,nzsize
!ops                DO jc = 1,nysize
!ops                    DO ic = 1,nxsize
!ops                        WRITE(ncdmpo,*)drun(ic,jc,kc),  &
!ops                            urun(ic,jc,kc),vrun(ic,jc,kc),wrun(ic,jc,kc), erun(ic,jc,kc),  &
!ops                            (yrun(ispec,ic,jc,kc),ispec=1,nspec)
!ops                    END DO
!ops                END DO
!ops            END DO
!ops            WRITE(ncdmpo,*)etime,tstep,errold,errldr
!ops            CLOSE(ncdmpo)
    
        END IF
  
#else
    CALL write_h5_dumpfile
#endif

!       REPORT THE DUMP
!       RSC 11-JUL-2009
        IF(iproc == 0) THEN
  
            OPEN(UNIT=ncrept,FILE=fnrept,STATUS='OLD',FORM='FORMATTED')
            3000      CONTINUE
            READ(ncrept,9000,END=3010)
            GO TO 3000
            3010      BACKSPACE(ncrept)
            WRITE(ncrept,9120)fndmpo(idflag+1)
            CLOSE(ncrept)

        END IF

!       RESET THE DUMP FLAG
        idflag = MOD(idflag+1,2)

    END IF

!   =========================================================================

!   DUMP BC INFORMATION AS REQUIRED
!   ===============================
    IF(MOD(itime,ntdump) == 0) THEN
  
        bcflag = (nsbcxl == nsbci2).OR.(nsbcxl == nsbci3)
        bcflag = bcflag.AND.(nxlprm(1) == 3)
  
        IF(bcflag) THEN
    
!           DUMP THE INLET TURBULENT VELOCITY FIELD
!ops            OPEN(UNIT=nctixl,FILE=fntixl,STATUS='OLD', FORM='UNFORMATTED')
!ops            REWIND(nctixl)
!ops            WRITE(nctixl)ufxl,vfxl,wfxl,slocxl,svelxl,bvelxl
!ops            CLOSE(nctixl)
    
        END IF
  
    END IF

!   =========================================================================

!   TIME STEP HISTORY
    IF(iproc == 0) THEN
        WRITE(*,'(I7,1PE12.4,I5)')itime,tstep,inderr
    END IF

!    IF(iproc == 0) THEN
        rangexyz = (/255,255,nyglbl,nyglbl,nzglbl,nzglbl/)
        call ops_par_loop(math_kernel_print_drhs, "print single value", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                        ops_arg_gbl(itime, 1, "integer", OPS_READ))

        rangexyz = (/256,256,nyglbl,nyglbl,nzglbl,nzglbl/)
        call ops_par_loop(math_kernel_print_erhs, "print single value", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_READ), &
                        ops_arg_gbl(itime, 1, "integer", OPS_READ))

        rangexyz = (/257,257,nyglbl,nyglbl,nzglbl,nzglbl/)
        call ops_par_loop(math_kernel_print_urhs, "print single value", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_READ), &
                        ops_arg_gbl(itime, 1, "integer", OPS_READ))
!    END IF

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

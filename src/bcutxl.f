      SUBROUTINE BCUTXL
 
C     *************************************************************************
C
C     BCUTXL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     30-DEC-2003:  CREATED
C     04-JAN-2007:  RSC REVISE PARALLEL RECEIVES
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
C     AND THEIR TIME DERIVATIVES
C
C     X-DIRECTION LEFT-HAND END
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
CKA   FIX INFLOW BUG, BTIME IS DEFINED IN COM_SENGA2.H
CKA      DOUBLE PRECISION BTIME
      DOUBLE PRECISION FORNOW,ARGMNT,ARGVAL,REALKX 
      DOUBLE PRECISION COSVAL,SINVAL,COSTHT,SINTHT
      DOUBLE PRECISION PCOUNT
      INTEGER IC,JC,KC
      INTEGER IIC,IIM,KX,KXBASE
      INTEGER ICPROC,NCOUNT,IRPROC,IRTAG


C     BEGIN
C     =====

C     =========================================================================

CKA   THIS WAS MOVED TO BOUNDT & BOUNTT TO FIX INFLOW SCANNING LOCATION
C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)
CKA      BTIME = ETIME + RKTIM(IRKSTP)

C     =========================================================================

C     CONSTANT U-VELOCITY
C     PARAMETER I1=1, R1=U-VELOCITY
      IF(NXLPRM(1).EQ.1)THEN

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRUXL(JC,KC) = RXLPRM(1)
C            FORNOW = REAL(JC)*DELTAY/(HALF*YGDLEN)
C            STRUXL(JC,KC) = RXLPRM(1)*TANH(FORNOW)
            STRVXL(JC,KC) = ZERO
            STRWXL(JC,KC) = ZERO

            DUDTXL(JC,KC) = ZERO
            DVDTXL(JC,KC) = ZERO
            DWDTXL(JC,KC) = ZERO

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     SINUSOIDAL U-VELOCITY
C     PARAMETER I1=2, R1=AMPLITUDE, R2=PERIOD
      IF(NXLPRM(1).EQ.2)THEN

        FORNOW = TWO*PI/RXLPRM(2)
        ARGMNT = FORNOW*BTIME

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRUXL(JC,KC) = RXLPRM(1)*SIN(ARGMNT)
            STRVXL(JC,KC) = ZERO
            STRWXL(JC,KC) = ZERO

            DUDTXL(JC,KC) = FORNOW*RXLPRM(1)*COS(ARGMNT)
            DVDTXL(JC,KC) = ZERO
            DWDTXL(JC,KC) = ZERO

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     TURBULENT VELOCITY FIELD
C     PARAMETER I1=3
      IF(NXLPRM(1).EQ.3)THEN

C       INTERPOLATE STORED TURBULENT VELOCITY FIELD ONTO INLET PLANE
C       DO THE INTERPOLATION BY DFT: LOCAL-PROCESSOR CONTRIBUTION

C       -----------------------------------------------------------------------

C       UPDATE THE SCANNING PLANE LOCATION
        SLOCXL = ELOCXL - SVELXL*BTIME
        IF(SLOCXL.LT.ZERO)SLOCXL = XGDLEN + SLOCXL
CKA     FIX INFLOW
CKA        IF(IRKSTP.EQ.NRKSTP)ELOCXL = SLOCXL
        IF(FUPELC)ELOCXL = SLOCXL

C       INITIALISE THE PHASE ANGLE TERMS
        ARGMNT = TPOVXG*SLOCXL
        COSTHT = COS(ARGMNT)
        SINTHT = SIN(ARGMNT)
        KXBASE = KMINXL

C       ZERO THE LOCAL-PROCESSOR CONTRIBUTION TO THE DFT
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRUXL(JC,KC) = ZERO
            STRVXL(JC,KC) = ZERO
            STRWXL(JC,KC) = ZERO

            DUDTXL(JC,KC) = ZERO
            DVDTXL(JC,KC) = ZERO
            DWDTXL(JC,KC) = ZERO

          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       SPECIAL CASE OF LEADING IMAGINARY TERM
        IF(FLLIXL)THEN

          KX = KXBASE
          REALKX = REAL(KX)
          ARGVAL = ARGMNT*REALKX
          COSVAL = COS(ARGVAL)
          SINVAL = SIN(ARGVAL)
          IIC = 1

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRUXL(JC,KC) = STRUXL(JC,KC)
     +                      + UFXL(IIC,JC,KC)*SINVAL
              STRVXL(JC,KC) = STRVXL(JC,KC)
     +                      + VFXL(IIC,JC,KC)*SINVAL
              STRWXL(JC,KC) = STRWXL(JC,KC)
     +                      + WFXL(IIC,JC,KC)*SINVAL

              DUDTXL(JC,KC) = DUDTXL(JC,KC)
     +                      - REALKX*UFXL(IIC,JC,KC)*COSVAL
              DVDTXL(JC,KC) = DVDTXL(JC,KC)
     +                      - REALKX*VFXL(IIC,JC,KC)*COSVAL
              DWDTXL(JC,KC) = DWDTXL(JC,KC)
     +                      - REALKX*WFXL(IIC,JC,KC)*COSVAL

            ENDDO
          ENDDO

          KXBASE = KXBASE + 1

        ENDIF

C       -----------------------------------------------------------------------

C       STANDARD LOCAL CONTRIBUTION

C       ZEROTH WAVENUMBER
        IF(KXBASE.EQ.0)THEN

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              KX = KXBASE
              REALKX = REAL(KX)
              ARGVAL = ARGMNT*REALKX
              COSVAL = COS(ARGVAL)
              SINVAL = SIN(ARGVAL)
              IIM = 1

              STRUXL(JC,KC) = STRUXL(JC,KC)
     +                      + HALF*UFXL(IIM,JC,KC)*COSVAL
              STRVXL(JC,KC) = STRVXL(JC,KC)
     +                      + HALF*VFXL(IIM,JC,KC)*COSVAL
              STRWXL(JC,KC) = STRWXL(JC,KC)
     +                      + HALF*WFXL(IIM,JC,KC)*COSVAL

            ENDDO
          ENDDO

          KXBASE = KXBASE + 1

        ENDIF

C       ALL OTHER WAVENUMBERS
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            KX = KXBASE
            REALKX = REAL(KX)
            ARGVAL = ARGMNT*REALKX
            COSVAL = COS(ARGVAL)
            SINVAL = SIN(ARGVAL)

            DO IC = ISTAXL,ISTOXL,2

              IIM = IC
              IIC = IC+1

              STRUXL(JC,KC) = STRUXL(JC,KC)
     +                      + UFXL(IIM,JC,KC)*COSVAL
     +                      + UFXL(IIC,JC,KC)*SINVAL
              STRVXL(JC,KC) = STRVXL(JC,KC)
     +                      + VFXL(IIM,JC,KC)*COSVAL
     +                      + VFXL(IIC,JC,KC)*SINVAL
              STRWXL(JC,KC) = STRWXL(JC,KC)
     +                      + WFXL(IIM,JC,KC)*COSVAL
     +                      + WFXL(IIC,JC,KC)*SINVAL

              DUDTXL(JC,KC) = DUDTXL(JC,KC)
     +                      + REALKX*(UFXL(IIM,JC,KC)*SINVAL
     +                              - UFXL(IIC,JC,KC)*COSVAL)
              DVDTXL(JC,KC) = DVDTXL(JC,KC)
     +                      + REALKX*(VFXL(IIM,JC,KC)*SINVAL
     +                              - VFXL(IIC,JC,KC)*COSVAL)
              DWDTXL(JC,KC) = DWDTXL(JC,KC)
     +                      + REALKX*(WFXL(IIM,JC,KC)*SINVAL
     +                              - WFXL(IIC,JC,KC)*COSVAL)

              KX = KX + 1
              REALKX = REAL(KX)
              FORNOW = COSVAL
              COSVAL = COSTHT*COSVAL - SINTHT*SINVAL
              SINVAL = SINTHT*FORNOW + COSTHT*SINVAL

            ENDDO

          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       SPECIAL CASE OF TRAILING REAL TERM
        IF(FLTRXL)THEN

          KX = KXBASE + ISTOXL/2
          REALKX = REAL(KX)
          ARGVAL = ARGMNT*REALKX
          COSVAL = COS(ARGVAL)
          SINVAL = SIN(ARGVAL)
          IIM = ISTOXL + 1

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRUXL(JC,KC) = STRUXL(JC,KC)
     +                      + UFXL(IIM,JC,KC)*COSVAL
              STRVXL(JC,KC) = STRVXL(JC,KC)
     +                      + VFXL(IIM,JC,KC)*COSVAL
              STRWXL(JC,KC) = STRWXL(JC,KC)
     +                      + WFXL(IIM,JC,KC)*COSVAL

              DUDTXL(JC,KC) = DUDTXL(JC,KC)
     +                      + REALKX*UFXL(IIM,JC,KC)*SINVAL
              DVDTXL(JC,KC) = DVDTXL(JC,KC)
     +                      + REALKX*VFXL(IIM,JC,KC)*SINVAL
              DWDTXL(JC,KC) = DWDTXL(JC,KC)
     +                      + REALKX*WFXL(IIM,JC,KC)*SINVAL

            ENDDO
          ENDDO

        ENDIF

C       -----------------------------------------------------------------------

C       PARALLEL TRANSFER
C       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
        IF(IXPROC.EQ.0)THEN

C         LEFTMOST PROCESSOR IN X
C         RECEIVE FROM ALL OTHER PROCESSORS IN X
          DO ICPROC = 1,NXPRM1

            IRPROC = NPROCX(ICPROC)
            IRTAG = IRPROC*NPROC+IPROC
            CALL P_RECV(PCOUNT,1,IRPROC,IRTAG)
            CALL P_RECV(PARRAY,NPARAY,IRPROC,IRTAG)

            NCOUNT = 0
            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                NCOUNT = NCOUNT + 1
                STRUXL(JC,KC) = STRUXL(JC,KC) + PARRAY(NCOUNT)
                NCOUNT = NCOUNT + 1
                STRVXL(JC,KC) = STRVXL(JC,KC) + PARRAY(NCOUNT)
                NCOUNT = NCOUNT + 1
                STRWXL(JC,KC) = STRWXL(JC,KC) + PARRAY(NCOUNT)
                NCOUNT = NCOUNT + 1
                DUDTXL(JC,KC) = DUDTXL(JC,KC) + PARRAY(NCOUNT)
                NCOUNT = NCOUNT + 1
                DVDTXL(JC,KC) = DVDTXL(JC,KC) + PARRAY(NCOUNT)
                NCOUNT = NCOUNT + 1
                DWDTXL(JC,KC) = DWDTXL(JC,KC) + PARRAY(NCOUNT)

              ENDDO
            ENDDO
          
          ENDDO

C         SCALING OF DFT
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

C             VELOCITIES
              STRUXL(JC,KC) = STRUXL(JC,KC)*SCAUXL
              STRVXL(JC,KC) = STRVXL(JC,KC)*SCAUXL
              STRWXL(JC,KC) = STRWXL(JC,KC)*SCAUXL

C             DERIVATIVES
              DUDTXL(JC,KC) = DUDTXL(JC,KC)*SCDUXL
              DVDTXL(JC,KC) = DVDTXL(JC,KC)*SCDUXL
              DWDTXL(JC,KC) = DWDTXL(JC,KC)*SCDUXL

C             ADD MEAN VELOCITY
              STRUXL(JC,KC) = STRUXL(JC,KC) + BVELXL

C             CONVERT SPATIAL TO TEMPORAL DERIVATIVES
              DUDTXL(JC,KC) = DUDTXL(JC,KC)*SVELXL
              DVDTXL(JC,KC) = DVDTXL(JC,KC)*SVELXL
              DWDTXL(JC,KC) = DWDTXL(JC,KC)*SVELXL

            ENDDO
          ENDDO

        ELSE 

C         NOT THE LEFTMOST PROCESSOR IN X
C         SEND TO LEFTMOST PROCESSOR IN X
          NCOUNT = 0
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              NCOUNT = NCOUNT + 1
              PARRAY(NCOUNT) = STRUXL(JC,KC)
              NCOUNT = NCOUNT + 1
              PARRAY(NCOUNT) = STRVXL(JC,KC)
              NCOUNT = NCOUNT + 1
              PARRAY(NCOUNT) = STRWXL(JC,KC)
              NCOUNT = NCOUNT + 1
              PARRAY(NCOUNT) = DUDTXL(JC,KC)
              NCOUNT = NCOUNT + 1
              PARRAY(NCOUNT) = DVDTXL(JC,KC)
              NCOUNT = NCOUNT + 1
              PARRAY(NCOUNT) = DWDTXL(JC,KC)

            ENDDO
          ENDDO

          PCOUNT = REAL(NCOUNT)
          IRPROC = NPROCX(0)
          IRTAG = IPROC*NPROC+IRPROC
          CALL P_SEND(PCOUNT,1,1,IRPROC,IRTAG)
          CALL P_SEND(PARRAY,NPARAY,NCOUNT,IRPROC,IRTAG)

        ENDIF

      ENDIF

C     =========================================================================


      RETURN
      END

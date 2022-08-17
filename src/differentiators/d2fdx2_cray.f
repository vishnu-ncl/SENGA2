      SUBROUTINE D2FDX2(FUNCTN,FDERIV)
 
C     *************************************************************************
C
C     D2FDX2
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     06-APR-2003:  RSC MODIFIED FOR SENGA2
 
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES SECOND X-DERIVATIVE OF SPECIFIED FUNCTION
C     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
C     EXPLICIT 8TH,6TH,4TH,4TH ORDER END CONDITIONS
C     WITH CRAY ARCHER PERFORMANCE IMPROVEMENTS
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     ARGUMENTS
C     =========
      DOUBLE PRECISION FUNCTN(NXBIGL:NXBIGR,NYBIGL:NYBIGR,NZBIGL:NZBIGR)
      DOUBLE PRECISION FDERIV(NXSIZE,NYSIZE,NZSIZE)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION FDIFAP,FDIFBP,FDIFCP,FDIFDP,FDIFEP
      DOUBLE PRECISION FDIFAM,FDIFBM,FDIFCM,FDIFDM,FDIFEM
      INTEGER IC,JC,KC
      INTEGER ISTART,IFINIS


C     BEGIN
C     =====

C     =========================================================================

C     END CONDITIONS
C     ==============

      ISTART = ISTAL
      IFINIS = ISTOL
      IF(NENDXL.EQ.NBOUND)ISTART = ISTAP5
      IF(NENDXR.EQ.NBOUND)IFINIS = ISTOM5

C     =========================================================================

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL

C         =====================================================================

C         INTERIOR SCHEME
C         ===============
          DO IC = ISTART,IFINIS

C           TENTH ORDER EXPLICIT DIFFERENCES
            FDIFAP = FUNCTN(IC+1,JC,KC) - FUNCTN(IC,JC,KC) 
            FDIFAM = FUNCTN(IC-1,JC,KC) - FUNCTN(IC,JC,KC)
            FDIFBP = FUNCTN(IC+2,JC,KC) - FUNCTN(IC,JC,KC) 
            FDIFBM = FUNCTN(IC-2,JC,KC) - FUNCTN(IC,JC,KC) 
            FDIFCP = FUNCTN(IC+3,JC,KC) - FUNCTN(IC,JC,KC) 
            FDIFCM = FUNCTN(IC-3,JC,KC) - FUNCTN(IC,JC,KC) 
            FDIFDP = FUNCTN(IC+4,JC,KC) - FUNCTN(IC,JC,KC) 
            FDIFDM = FUNCTN(IC-4,JC,KC) - FUNCTN(IC,JC,KC) 
            FDIFEP = FUNCTN(IC+5,JC,KC) - FUNCTN(IC,JC,KC) 
            FDIFEM = FUNCTN(IC-5,JC,KC) - FUNCTN(IC,JC,KC) 

            FDERIV(IC,JC,KC) = ACOFSX*(FDIFAP+FDIFAM)
     +                       + BCOFSX*(FDIFBP+FDIFBM)
     +                       + CCOFSX*(FDIFCP+FDIFCM)
     +                       + DCOFSX*(FDIFDP+FDIFDM)
     +                       + ECOFSX*(FDIFEP+FDIFEM)

            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDLX2

          ENDDO

C         =====================================================================

C         LH END
C         ======
          IF(NENDXL.EQ.NBOUND)THEN

C           EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT

C           LH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(ISTAP1,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDIFBP = FUNCTN(ISTAP2,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDIFCP = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDIFDP = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDIFEP = FUNCTN(ISTAP5,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAL,JC,KC) = (ACFS1X*FDIFAP
     +                          +  BCFS1X*FDIFBP
     +                          +  CCFS1X*FDIFCP
     +                          +  DCFS1X*FDIFDP
     +                          +  ECFS1X*FDIFEP)*OVDLX2

C           LH POINT PLUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(ISTAL,JC,KC)  - FUNCTN(ISTAP1,JC,KC) 
            FDIFBP = FUNCTN(ISTAP2,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFCP = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFDP = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFEP = FUNCTN(ISTAP5,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDERIV(ISTAP1,JC,KC) = (ACFS2X*FDIFAP
     +                           +  BCFS2X*FDIFBP
     +                           +  CCFS2X*FDIFCP
     +                           +  DCFS2X*FDIFDP
     +                           +  ECFS2X*FDIFEP)*OVDLX2

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAP2,JC,KC) 
            FDIFAM = FUNCTN(ISTAP2,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFBP = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP2,JC,KC) 
            FDIFBM = FUNCTN(ISTAP2,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAP2,JC,KC) = (ACFS3X*(FDIFAP-FDIFAM)
     +                           +  BCFS3X*(FDIFBP-FDIFBM))*OVDLX2

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP3,JC,KC) 
            FDIFAM = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAP2,JC,KC) 
            FDIFBP = FUNCTN(ISTAP5,JC,KC) - FUNCTN(ISTAP3,JC,KC) 
            FDIFBM = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFCP = FUNCTN(ISTAP6,JC,KC) - FUNCTN(ISTAP3,JC,KC) 
            FDIFCM = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAP3,JC,KC) = (ACFS4X*(FDIFAP-FDIFAM)
     +                           +  BCFS4X*(FDIFBP-FDIFBM)
     +                           +  CCFS4X*(FDIFCP-FDIFCM))*OVDLX2

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTAP5,JC,KC) - FUNCTN(ISTAP4,JC,KC) 
            FDIFAM = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP3,JC,KC) 
            FDIFBP = FUNCTN(ISTAP6,JC,KC) - FUNCTN(ISTAP4,JC,KC) 
            FDIFBM = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP2,JC,KC) 
            FDIFCP = FUNCTN(ISTAP7,JC,KC) - FUNCTN(ISTAP4,JC,KC) 
            FDIFCM = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFDP = FUNCTN(ISTAP8,JC,KC) - FUNCTN(ISTAP4,JC,KC) 
            FDIFDM = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAP4,JC,KC) = (ACFS5X*(FDIFAP-FDIFAM)
     +                           +  BCFS5X*(FDIFBP-FDIFBM)
     +                           +  CCFS5X*(FDIFCP-FDIFCM)
     +                           +  DCFS5X*(FDIFDP-FDIFDM))*OVDLX2

          ENDIF 

C         =====================================================================

C         RH END
C         ======
          IF(NENDXR.EQ.NBOUND)THEN

C           EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTOM3,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDIFAM = FUNCTN(ISTOM4,JC,KC) - FUNCTN(ISTOM5,JC,KC) 
            FDIFBP = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDIFBM = FUNCTN(ISTOM4,JC,KC) - FUNCTN(ISTOM6,JC,KC) 
            FDIFCP = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDIFCM = FUNCTN(ISTOM4,JC,KC) - FUNCTN(ISTOM7,JC,KC) 
            FDIFDP = FUNCTN(ISTOL,JC,KC)  - FUNCTN(ISTOM4,JC,KC) 
            FDIFDM = FUNCTN(ISTOM4,JC,KC) - FUNCTN(ISTOM8,JC,KC) 
            FDERIV(ISTOM4,JC,KC) = (ACFS5X*(FDIFAP-FDIFAM)
     +                           +  BCFS5X*(FDIFBP-FDIFBM)
     +                           +  CCFS5X*(FDIFCP-FDIFCM)
     +                           +  DCFS5X*(FDIFDP-FDIFDM))*OVDLX2
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM3,JC,KC) 
            FDIFAM = FUNCTN(ISTOM3,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDIFBP = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM3,JC,KC) 
            FDIFBM = FUNCTN(ISTOM3,JC,KC) - FUNCTN(ISTOM5,JC,KC) 
            FDIFCP = FUNCTN(ISTOL,JC,KC)  - FUNCTN(ISTOM3,JC,KC) 
            FDIFCM = FUNCTN(ISTOM3,JC,KC) - FUNCTN(ISTOM6,JC,KC) 
            FDERIV(ISTOM3,JC,KC) = (ACFS4X*(FDIFAP-FDIFAM)
     +                           +  BCFS4X*(FDIFBP-FDIFBM)
     +                           +  CCFS4X*(FDIFCP-FDIFCM))*OVDLX2
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM2,JC,KC) 
            FDIFAM = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM3,JC,KC) 
            FDIFBP = FUNCTN(ISTOL,JC,KC)  - FUNCTN(ISTOM2,JC,KC) 
            FDIFBM = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDERIV(ISTOM2,JC,KC) = (ACFS3X*(FDIFAP-FDIFAM)
     +                           +  BCFS3X*(FDIFBP-FDIFBM))*OVDLX2
      
C           RH POINT MINUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(ISTOL,JC,KC)  - FUNCTN(ISTOM1,JC,KC) 
            FDIFBP = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM1,JC,KC) 
            FDIFCP = FUNCTN(ISTOM3,JC,KC) - FUNCTN(ISTOM1,JC,KC) 
            FDIFDP = FUNCTN(ISTOM4,JC,KC) - FUNCTN(ISTOM1,JC,KC) 
            FDIFEP = FUNCTN(ISTOM5,JC,KC) - FUNCTN(ISTOM1,JC,KC) 
            FDERIV(ISTOM1,JC,KC) = (ACFS2X*FDIFAP
     +                           +  BCFS2X*FDIFBP
     +                           +  CCFS2X*FDIFCP
     +                           +  DCFS2X*FDIFDP
     +                           +  ECFS2X*FDIFEP)*OVDLX2
      
C           RH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(ISTOM1,JC,KC)  - FUNCTN(ISTOL,JC,KC) 
            FDIFBP = FUNCTN(ISTOM2,JC,KC)  - FUNCTN(ISTOL,JC,KC) 
            FDIFCP = FUNCTN(ISTOM3,JC,KC)  - FUNCTN(ISTOL,JC,KC) 
            FDIFDP = FUNCTN(ISTOM4,JC,KC)  - FUNCTN(ISTOL,JC,KC) 
            FDIFEP = FUNCTN(ISTOM5,JC,KC)  - FUNCTN(ISTOL,JC,KC) 
            FDERIV(ISTOL,JC,KC) = (ACFS1X*FDIFAP
     +                          +  BCFS1X*FDIFBP
     +                          +  CCFS1X*FDIFCP
     +                          +  DCFS1X*FDIFDP
     +                          +  ECFS1X*FDIFEP)*OVDLX2

          ENDIF

C         =====================================================================

        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END

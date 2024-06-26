      SUBROUTINE FFTNRA(NX,NI,IFORW,NRADIX)

C     *************************************************************************
C
C     FFTNRA
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     23-MAY-1999:  CREATED
C
C     DESCRIPTION
C     -----------
C     CARRIES OUT AN FFT IN 1D
C     USING TEMPERTON SELF-SORTING PRIME FACTOR ALGORITHM
C     ARBITRARY RADIX
C     WITH ROTATIONS
C     INITIALISATION ONLY
C
C     REFERENCES
C     ----------
C     TEMPERTON C.: SIAM J SCI STAT COMP 13,3,676-686, 1992
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      DOUBLE PRECISION ONE,EIGHT
      PARAMETER(ONE=1.0D0, EIGHT=8.0D0)


C     GLOBAL DATA
C     ===========
C     FFACOM-------------------------------------------------------------------

      INTEGER NFTMAX
      PARAMETER(NFTMAX=1024)

      INTEGER NRDXMX,NRDXM1
      PARAMETER(NRDXMX=11,NRDXM1=NRDXMX-1)

      DOUBLE PRECISION COSTMA(NFTMAX),SINTMA(NFTMAX)
      DOUBLE PRECISION COEFRA(NRDXMX),COEFIA(NRDXMX)

      INTEGER NXFRAC,MPOWRA,MPHLFA,NIA,INCRMA,NRADXA,NRADM1

      COMMON/FFACOM/COSTMA,SINTMA,
     +              COEFRA,COEFIA,
     +              NXFRAC,MPOWRA,MPHLFA,NIA,INCRMA,NRADXA,NRADM1

C     FFACOM-------------------------------------------------------------------


C     ARGUMENTS
C     =========
      INTEGER NX,NI,IFORW,NRADIX


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ARGMNT,TWOPIN,TWOPIR,REALNX,RLOGR,RRADIX
      INTEGER IADD,IX,II,IROTN,IROTR


C     BEGIN
C     =====

C     USEFUL CONSTANTS
      RRADIX = REAL(NRADIX)
      REALNX = REAL(NI)
      RLOGR  = LOG(RRADIX)
      TWOPIN = REAL(IFORW)*EIGHT*ATAN(ONE)
      TWOPIR = TWOPIN/RRADIX
      TWOPIN = TWOPIN/REALNX
      NXFRAC = NX/NRADIX
      MPOWRA = NINT(LOG(REALNX)/RLOGR)
      IADD = MOD(MPOWRA,2)
      MPHLFA = MPOWRA/2+IADD
      NRADM1 = NRADIX-1
      INCRMA = NX/NI
      IROTN = MOD(INCRMA,NI)
      IROTR = MOD(IROTN,NRADIX)
      TWOPIR = TWOPIR*REAL(IROTR)
      NRADXA = NRADIX
      NIA = NI


C     PRECOMPUTE THE RADIX-N COEFFICIENTS
      DO IX = 1,NRADM1
        ARGMNT = TWOPIR*REAL(IX)
        COEFRA(IX) = COS(ARGMNT)
        COEFIA(IX) = SIN(ARGMNT)
      ENDDO

C     PRECOMPUTE THE TRIG FACTORS
      II = 0
      DO IX = 1,NRADM1*(NXFRAC-1)+1
        ARGMNT = TWOPIN*REAL(II)
        COSTMA(IX) = COS(ARGMNT)
        SINTMA(IX) = SIN(ARGMNT)
        II = II+IROTN
        IF(II.GT.NI)II = II-NI
      ENDDO


      RETURN
      END

      SUBROUTINE BCDTXR
 
C     *************************************************************************
C
C     BCDTXR
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     30-DEC-2003:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR DENSITY
C     AND ITS TIME DERIVATIVE
C
C     X-DIRECTION RIGHT-HAND END
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      USE com_senga
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER JC,KC


C     BEGIN
C     =====

C     =========================================================================

C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

C     =========================================================================

C     EVALUATE AND RETURN STRDXR,DDDTXR

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL

C         SET DENSITY TO CONSTANT (INITIAL) VALUE
          STRDXR(JC,KC) = DRIN

C         SET DENSITY TIME DERIVATIVE TO ZERO
          DDDTXR(JC,KC) = ZERO

        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END

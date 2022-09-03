      SUBROUTINE BCTTZR
 
C     *************************************************************************
C
C     BCTTZR
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     26-OCT-2013:  CREATED
C     09-MAY-2015:  RSC MODIFIED FOR ISOTHERMAL WALL
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR TEMPERATURE
C     AND ITS TIME DERIVATIVE
C
C     Z-DIRECTION RIGHT-HAND END
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      USE com_senga
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER IC,JC


C     BEGIN
C     =====

C     =========================================================================

C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

C     =========================================================================

C     EVALUATE AND RETURN STRTZR,DTDTZR
      DO JC = JSTAL,JSTOL
        DO IC = ISTAL,ISTOL

          STRTZR(IC,JC) = TRIN

          DTDTZR(IC,JC) = ZERO

        ENDDO
      ENDDO

C     =========================================================================

C     ISOTHERMAL WALL
      IF(NSBCZR.EQ.NSBCW2)THEN

        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STRTZR(IC,JC) = RZRPRM(1)

            DTDTZR(IC,JC) = ZERO

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================


      RETURN
      END

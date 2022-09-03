      SUBROUTINE BCUTZL
 
C     *************************************************************************
C
C     BCUTZL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     26-OCT-2013:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
C     AND THEIR TIME DERIVATIVES
C
C     Z-DIRECTION LEFT-HAND END
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

C     CONSTANT W-VELOCITY
C     PARAMETER I1=1, R1=W-VELOCITY
      IF(NZLPRM(1).EQ.1)THEN

        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STRUZL(IC,JC) = ZERO
            STRVZL(IC,JC) = ZERO
            STRWZL(IC,JC) = RZLPRM(1)

            DUDTZL(IC,JC) = ZERO
            DVDTZL(IC,JC) = ZERO
            DWDTZL(IC,JC) = ZERO

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================


      RETURN
      END

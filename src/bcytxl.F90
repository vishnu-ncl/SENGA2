SUBROUTINE bcytxl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2024-11-08  Time: 05:06:48

!     *************************************************************************

!     BCYTXL
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     30-DEC-2003:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR MASS FRACTIONS
!     AND THEIR TIME DERIVATIVES

!     X-DIRECTION LEFT-HAND END

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
INTEGER :: jc,kc
INTEGER :: ispec
DOUBLE PRECISION :: toty


!     BEGIN
!     =====

!     =========================================================================

!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!     =========================================================================

!     EVALUATE AND RETURN STRYXL,DYDTXL
DO ispec = 1,nspec
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
!           SET MASS FRACTIONS TO CONSTANT (INITIAL) VALUES
      stryxl(jc,kc,ispec) = yrin(ispec)
      
!           SET MASS FRACTION TIME DERIVATIVES TO ZERO
      dydtxl(jc,kc,ispec) = zero
      
    END DO
  END DO
  
END DO

!     VM: SYNTHETIC SCALAR INFLOW
!     VM: NXLPRM(2)=1 IMPLIES THAT THE SCALAR SYTHETIC DIGITAL FILTERING
!     IS ON
IF ((nxlprm(2)==1).AND.(nxlprm(1)==4).AND.(ngbcxl==12))THEN
  DO ispec=1,nspec
    DO kc=kstal,kstol
      DO jc=jstal,jstol
        stryxl(jc,kc,ispec)=yrin(ispec)+yinf2(jc,kc,ispec)
        IF(stryxl(jc,kc,ispec) > 1.0D0) THEN
          yinf2(jc,kc,ispec)=1.0D0-yrin(ispec)
          stryxl(jc,kc,ispec)=1.0D0
        END IF
        IF(stryxl(jc,kc,ispec) < 0.0D0) THEN
          yinf2(jc,kc,ispec)=yrin(ispec)-0.0D0
          stryxl(jc,kc,ispec)=0.0D0
        END IF
        dydtxl(jc,kc,ispec)=(yinf2(jc,kc,ispec)- yinf1(jc,kc,ispec))/tstep
      END DO
    END DO
  END DO
  DO kc=kstal,kstol
    DO jc=jstal,jstol
      toty=0.0D0
      DO ispec=1,nspec-1
        toty = toty+stryxl(jc,kc,ispec)
      END DO
      stryxl(jc,kc,nspec)=1.0D0-toty
      stryxl(jc,kc,nspec)=MAX(0.0,stryxl(jc,kc,nspec))
      stryxl(jc,kc,nspec)=MIN(1.0,stryxl(jc,kc,nspec))
      dydtxl(jc,kc,nspec)=(yinf2(jc,kc,nspec) -yinf1(jc,kc,nspec))/tstep
    END DO
  END DO
END IF

!     =========================================================================


RETURN
END SUBROUTINE bcytxl

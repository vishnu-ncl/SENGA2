SUBROUTINE bcytxl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:15:08

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
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
INTEGER :: jc,kc
INTEGER :: ispec


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
      stryxl(ispec,1,jc,kc) = yrin(ispec)
      
!           SET MASS FRACTION TIME DERIVATIVES TO ZERO
      dydtxl(ispec,1,jc,kc) = zero
      
    END DO
  END DO
  
END DO

!     =========================================================================


RETURN
END SUBROUTINE bcytxl

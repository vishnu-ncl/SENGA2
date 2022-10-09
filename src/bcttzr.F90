SUBROUTINE bcttzr
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-13  Time: 20:59:13

!     *************************************************************************

!     BCTTZR
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     26-OCT-2013:  CREATED
!     09-MAY-2015:  RSC MODIFIED FOR ISOTHERMAL WALL

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR TEMPERATURE
!     AND ITS TIME DERIVATIVE

!     Z-DIRECTION RIGHT-HAND END

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
INTEGER :: ic,jc


!     BEGIN
!     =====

!     =========================================================================

!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!     =========================================================================

!     EVALUATE AND RETURN STRTZR,DTDTZR
DO jc = jstal,jstol
  DO ic = istal,istol
    
    strtzr(ic,jc,1) = trin
    
    dtdtzr(ic,jc) = zero
    
  END DO
END DO

!     =========================================================================

!     ISOTHERMAL WALL
IF(nsbczr == nsbcw2)THEN
  
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      strtzr(ic,jc,1) = rzrprm(1)
      
      dtdtzr(ic,jc) = zero
      
    END DO
  END DO
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE bcttzr

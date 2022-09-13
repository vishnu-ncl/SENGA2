SUBROUTINE bcttxl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-13  Time: 20:58:53

!     *************************************************************************

!     BCTTXL
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     30-DEC-2003:  CREATED
!     09-MAY-2015:  RSC MODIFIED FOR ISOTHERMAL WALL

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR TEMPERATURE
!     AND ITS TIME DERIVATIVE

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


!     BEGIN
!     =====

!     =========================================================================

!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!     =========================================================================

!     EVALUATE AND RETURN STRTXL,DTDTXL
DO kc = kstal,kstol
  DO jc = jstal,jstol
    
    strtxl(jc,kc) = trin
    
    dtdtxl(jc,kc) = zero
    
  END DO
END DO

!     =========================================================================

!     ISOTHERMAL WALL
IF(nsbcxl == nsbcw2)THEN
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      strtxl(jc,kc) = rxlprm(1)
      
      dtdtxl(jc,kc) = zero
      
    END DO
  END DO
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE bcttxl

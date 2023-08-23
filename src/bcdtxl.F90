SUBROUTINE bcdtxl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-13  Time: 20:58:19

!     *************************************************************************

!     BCDTXL
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
!     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR DENSITY
!     AND ITS TIME DERIVATIVE

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


!     BEGIN
!     =====

!     =========================================================================

!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!     =========================================================================

!     EVALUATE AND RETURN STRDXL,DDDTXL

DO kc = kstal,kstol
  DO jc = jstal,jstol
    
!         SET DENSITY TO CONSTANT (INITIAL) VALUE
    strdxl(jc,kc) = drin
    
!         SET DENSITY TIME DERIVATIVE TO ZERO
    dddtxl(jc,kc) = zero
    
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE bcdtxl

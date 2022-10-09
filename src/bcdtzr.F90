SUBROUTINE bcdtzr
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-13  Time: 20:58:46

!     *************************************************************************

!     BCDTZR
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     26-OCT-2013:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR DENSITY
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

!     EVALUATE AND RETURN STRDZR,DDDTZR

DO jc = jstal,jstol
  DO ic = istal,istol
    
!         SET DENSITY TO CONSTANT (INITIAL) VALUE
    strdzr(ic,jc,1) = drin
    
!         SET DENSITY TIME DERIVATIVE TO ZERO
    dddtzr(ic,jc) = zero
    
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE bcdtzr

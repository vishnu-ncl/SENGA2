SUBROUTINE bcdtyr
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-13  Time: 20:58:32

!     *************************************************************************

!     BCDTYR
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

!     Y-DIRECTION RIGHT-HAND END

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
INTEGER :: ic,kc


!     BEGIN
!     =====

!     =========================================================================

!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!     =========================================================================

!     EVALUATE AND RETURN STRDYR,DDDTYR

DO kc = kstal,kstol
  DO ic = istal,istol
    
!         SET DENSITY TO CONSTANT (INITIAL) VALUE
    strdyr(ic,kc) = drin
    
!         SET DENSITY TIME DERIVATIVE TO ZERO
    dddtyr(ic,kc) = zero
    
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE bcdtyr

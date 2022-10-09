SUBROUTINE bcdtyl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-13  Time: 20:58:28

!     *************************************************************************

!     BCDTYL
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

!     Y-DIRECTION LEFT-HAND END

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

!     EVALUATE AND RETURN STRDYL,DDDTYL

DO kc = kstal,kstol
  DO ic = istal,istol
    
!         SET DENSITY TO CONSTANT (INITIAL) VALUE
    strdyl(ic,1,kc) = drin
    
!         SET DENSITY TIME DERIVATIVE TO ZERO
    dddtyl(ic,kc) = zero
    
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE bcdtyl

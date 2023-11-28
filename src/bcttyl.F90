SUBROUTINE bcttyl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-13  Time: 20:59:02

!     *************************************************************************

!     BCTTYL
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

!     Y-DIRECTION LEFT-HAND END

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

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

!     EVALUATE AND RETURN STRTYL,DTDTYL
DO kc = kstal,kstol
  DO ic = istal,istol
    
    strtyl(ic,kc) = trin
    
    dtdtyl(ic,kc) = zero
    
  END DO
END DO

!     =========================================================================

!     ISOTHERMAL WALL
IF(nsbcyl == nsbcw2)THEN
  
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      strtyl(ic,kc) = rylprm(1)
      
      dtdtyl(ic,kc) = zero
      
    END DO
  END DO
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE bcttyl

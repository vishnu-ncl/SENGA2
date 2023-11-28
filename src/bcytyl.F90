SUBROUTINE bcytyl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:15:16

!     *************************************************************************

!     BCYTYL
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
!     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR MASS FRACTIONS
!     AND THEIR TIME DERIVATIVES

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
INTEGER :: ispec


!     BEGIN
!     =====

!     =========================================================================

!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!     =========================================================================

!     EVALUATE AND RETURN STRYYL,DYDTYL
DO ispec = 1,nspec
  
  DO kc = kstal,kstol
    DO ic = istal,istol
      
!           SET MASS FRACTIONS TO CONSTANT (INITIAL) VALUES
      stryyl(ic,kc,ispec) = yrin(ispec)
      
!           SET MASS FRACTION TIME DERIVATIVES TO ZERO
      dydtyl(ic,kc,ispec) = zero
      
    END DO
  END DO
  
END DO

!     =========================================================================


RETURN
END SUBROUTINE bcytyl

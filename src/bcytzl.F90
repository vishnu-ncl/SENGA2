SUBROUTINE bcytzl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:24:27

!     *************************************************************************

!     BCYTZL
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

!     Z-DIRECTION LEFT-HAND END

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
INTEGER :: ispec


!     BEGIN
!     =====

!     =========================================================================

!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!     =========================================================================

!     EVALUATE AND RETURN STRYZL,DYDTZL
DO ispec = 1,nspec
  
  DO jc = jstal,jstol
    DO ic = istal,istol
      
!           SET MASS FRACTIONS TO CONSTANT (INITIAL) VALUES
      stryzl(ic,jc,ispec) = yrin(ispec)
      
!           SET MASS FRACTION TIME DERIVATIVES TO ZERO
      dydtzl(ic,jc,ispec) = zero
      
    END DO
  END DO
  
END DO

!     =========================================================================


RETURN
END SUBROUTINE bcytzl

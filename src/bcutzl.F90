SUBROUTINE bcutzl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:15:01

!     *************************************************************************

!     BCUTZL
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
!     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
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


!     BEGIN
!     =====

!     =========================================================================

!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!     =========================================================================

!     CONSTANT W-VELOCITY
!     PARAMETER I1=1, R1=W-VELOCITY
IF(nzlprm(1) == 1)THEN
  
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      struzl(ic,jc,1) = zero
      strvzl(ic,jc,1) = zero
      strwzl(ic,jc,1) = rzlprm(1)
      
      dudtzl(ic,jc) = zero
      dvdtzl(ic,jc) = zero
      dwdtzl(ic,jc) = zero
      
    END DO
  END DO
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE bcutzl

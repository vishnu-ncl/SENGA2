SUBROUTINE bcutyr
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:14:57

!     *************************************************************************

!     BCUTYR
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

!     CONSTANT V-VELOCITY
!     PARAMETER I1=1, R1=V-VELOCITY
IF(nyrprm(1) == 1)THEN
  
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      struyr(ic,kc) = zero
      strvyr(ic,kc) = ryrprm(1)
      strwyr(ic,kc) = zero
      
      dudtyr(ic,kc) = zero
      dvdtyr(ic,kc) = zero
      dwdtyr(ic,kc) = zero
      
    END DO
  END DO
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE bcutyr

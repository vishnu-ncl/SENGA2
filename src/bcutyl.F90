SUBROUTINE bcutyl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:14:52

!     *************************************************************************

!     BCUTYL
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

!     CONSTANT V-VELOCITY
!     PARAMETER I1=1, R1=V-VELOCITY
IF(nylprm(1) == 1)THEN
  
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      struyl(ic,kc) = zero
      strvyl(ic,kc) = rylprm(1)
      strwyl(ic,kc) = zero
      
      dudtyl(ic,kc) = zero
      dvdtyl(ic,kc) = zero
      dwdtyl(ic,kc) = zero
      
    END DO
  END DO
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE bcutyl

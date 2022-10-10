SUBROUTINE bcutxr
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:14:49

!     *************************************************************************

!     BCUTXR
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
!     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
!     AND THEIR TIME DERIVATIVES

!     X-DIRECTION RIGHT-HAND END

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
!KA   FIX INFLOW BC
!KA      REAL(KIND=dp) BTIME
REAL(KIND=dp) :: fornow,argmnt
INTEGER :: jc,kc


!     BEGIN
!     =====

!     =========================================================================

!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)
!KA   FIX INFLOW BC
!KA      BTIME = ETIME + RKTIM(IRKSTP)

!     =========================================================================

!     CONSTANT U-VELOCITY
!     PARAMETER I1=1, R1=U-VELOCITY
IF(nxrprm(1) == 1)THEN
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      struxr(1,jc,kc) = rxrprm(1)
      strvxr(1,jc,kc) = zero
      strwxr(1,jc,kc) = zero
      
      dudtxr(jc,kc) = zero
      dvdtxr(jc,kc) = zero
      dwdtxr(jc,kc) = zero
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     SINUSOIDAL U-VELOCITY
!     PARAMETER I1=2, R1=AMPLITUDE, R2=PERIOD
IF(nxrprm(1) == 2)THEN
  
  fornow = two*pi/rxrprm(2)
  argmnt = fornow*btime
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      struxr(1,jc,kc) = rxrprm(1)*SIN(argmnt)
      strvxr(1,jc,kc) = zero
      strwxr(1,jc,kc) = zero
      
      dudtxr(jc,kc) = fornow*rxrprm(1)*COS(argmnt)
      dvdtxr(jc,kc) = zero
      dwdtxr(jc,kc) = zero
      
    END DO
  END DO
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE bcutxr

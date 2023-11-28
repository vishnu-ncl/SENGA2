SUBROUTINE integf(functn,argmin,argmax,answer)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:26:02

!     *************************************************************************

!     INTEGF
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     20-APR-1997:  CREATED
!     24-MAY-2003:  UPDATED FOR SENGA2
!     15-SEP-2012:  RSC IMPROVE INITIAL SAMPLING OF FUNCTION

!     DESCRIPTION
!     -----------
!     INTEGRATES THE SPECIFIED FUNCTION BETWEEN THE SPECIFIED LIMITS
!     USING SUCCESSIVELY-REFINED SIMPSON'S RULE

!     REFERENCE
!     ---------
!     NUMERICAL RECIPES 1st Ed, CUP (1986), pp110-113

!     *************************************************************************


!     PARAMETERS
!     ==========

REAL(kind=8), INTENT(IN)             :: argmin
REAL(kind=8), INTENT(IN OUT)         :: argmax
REAL(kind=8), INTENT(OUT)            :: answer


REAL(kind=8), PARAMETER :: tolint=0.000001_8
INTEGER, PARAMETER :: intmax=20
REAL(kind=8), PARAMETER :: vsmall = 0.000000000000000000000000000001_8
REAL(kind=8), PARAMETER :: zero=0.0_8
REAL(kind=8), PARAMETER :: half=0.5_8
REAL(kind=8), PARAMETER :: three=3.0_8
REAL(kind=8), PARAMETER :: four=4.0_8


!     ARGUMENTS
!     =========

EXTERNAL functn


!     LOCAL DATA
!     ==========
REAL(kind=8) :: oldtmp,oldans,tmpans,argtmp,deltrg,addpts
INTEGER :: icount,ic,interv


!     BEGIN
!     =====

!     =========================================================================

!     RSC 15-SEP-2012 IMPROVE INITIAL SAMPLING OF FUNCTION
!      INTERV = 1
interv = 4
icount = 1
answer = vsmall
tmpans = half*(argmax-argmin)*(functn(argmin)+functn(argmax))

!     =========================================================================

!     MAIN LOOP TO REFINE THE INTEGRAL
!     --------------------------------

1000  CONTINUE

oldans = answer
oldtmp = tmpans

deltrg = (argmax-argmin)/REAL(interv,kind=8)
argtmp = argmin-half*deltrg

addpts = zero
DO ic = 1, interv
  argtmp = argtmp + deltrg
  addpts = addpts + functn(argtmp)
END DO

tmpans = half*(tmpans+addpts*deltrg)
answer = (four*tmpans-oldtmp)/three

IF(ABS(answer-oldans) > (tolint*ABS(oldans)))THEN
  icount = icount + 1
  IF(icount <= intmax)THEN
    interv = 2*interv
    GO TO 1000
  ELSE
!           FINISH EVEN IF INTEGRAL IS NOT CONVERGED
!           RSC 15-SEP-2012 MINOR BUG FIX
    WRITE(6,*)'Warning: INTEGF: integral not converged'
    WRITE(6,*)'at iteration:',intmax
    WRITE(6,*)'with values',answer,oldans
  END IF
END IF

!     END OF MAIN LOOP
!     ----------------

!     =========================================================================


RETURN
END SUBROUTINE integf

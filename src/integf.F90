SUBROUTINE integf(functn,argmin,argmax,answer)

!   *************************************************************************

!   INTEGF
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   20-APR-1997:  CREATED
!   24-MAY-2003:  UPDATED FOR SENGA2
!   15-SEP-2012:  RSC IMPROVE INITIAL SAMPLING OF FUNCTION

!   DESCRIPTION
!   -----------
!   INTEGRATES THE SPECIFIED FUNCTION BETWEEN THE SPECIFIED LIMITS
!   USING SUCCESSIVELY-REFINED SIMPSON'S RULE

!   REFERENCE
!   ---------
!   NUMERICAL RECIPES 1st Ed, CUP (1986), pp110-113

!   *************************************************************************

!   PARAMETERS
!   ==========
    real(kind=8), intent(in)             :: argmin
    real(kind=8), intent(in out)         :: argmax
    real(kind=8), intent(out)            :: answer


    real(kind=8), parameter :: tolint=0.000001_8
    integer(kind=4), parameter :: intmax=20
    real(kind=8), parameter :: vsmall = 0.000000000000000000000000000001_8
    real(kind=8), parameter :: zero=0.0_8
    real(kind=8), parameter :: half=0.5_8
    real(kind=8), parameter :: three=3.0_8
    real(kind=8), parameter :: four=4.0_8

!   ARGUMENTS
!   =========
    real(kind=8) :: functn
    EXTERNAL functn

!   LOCAL DATA
!   ==========
    real(kind=8) :: oldtmp,oldans,tmpans,argtmp,deltrg,addpts
    integer(kind=4) :: icount,ic,interv

!   BEGIN
!   =====

!   =========================================================================

!   RSC 15-SEP-2012 IMPROVE INITIAL SAMPLING OF FUNCTION
!   INTERV = 1
    interv = 4
    icount = 1
    answer = vsmall
    tmpans = half*(argmax-argmin)*(functn(argmin)+functn(argmax))

!   =========================================================================

!   MAIN LOOP TO REFINE THE INTEGRAL
!   --------------------------------

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

    IF(ABS(answer-oldans) > (tolint*ABS(oldans))) THEN
        icount = icount + 1
        IF(icount <= intmax) THEN
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

!   END OF MAIN LOOP
!   ----------------

!   =========================================================================

END SUBROUTINE integf

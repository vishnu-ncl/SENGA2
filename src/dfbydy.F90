SUBROUTINE dfbydy(functn,fderiv)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-03  Time: 21:30:59

!     *************************************************************************

!     DFBYDY
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     11-APR-2003:  RSC MODIFIED FOR SENGA2
!     10-OCT-2004:  RSC NULL VERSION

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES FIRST Y-DERIVATIVE OF SPECIFIED FUNCTION

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     ARGUMENTS
!     =========

real(kind=8), INTENT(IN OUT)         :: functn(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
real(kind=8), INTENT(OUT)            :: fderiv(nxsize,nysize,nzsize)




!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc


!     BEGIN
!     =====

!     =========================================================================

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fderiv(ic,jc,kc) = zero
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE dfbydy

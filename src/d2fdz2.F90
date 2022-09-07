SUBROUTINE d2fdz2(functn,fderiv)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-05  Time: 11:29:16

!     *************************************************************************

!     D2FDZ2
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
!     EVALUATES SECOND Z-DERIVATIVE OF SPECIFIED FUNCTION
!     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
!     EXPLICIT 8TH,6TH,4TH,2ND ORDER END CONDITIONS

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!     -------------------------------------------------------------------------


!     ARGUMENTS
!     =========

REAL(KIND=8), INTENT(IN OUT)         :: functn(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
REAL(KIND=8), INTENT(OUT)            :: fderiv(nxsize,nysize,nzsize)




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
END SUBROUTINE d2fdz2

SUBROUTINE d2fdyz(functn,fderiv)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-05  Time: 11:29:11

!     *************************************************************************

!     D2FDYZ
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     15-MAY-2003:  RSC MODIFIED FOR SENGA2
!     10-OCT-2004:  RSC NULL VERSION
!     31-DEC-2006:  RSC BUG FIX INDICES

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES SECOND YZ-DERIVATIVE OF SPECIFIED FUNCTION

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

!     RSC 31-DEC-2006 BUG FIX INDICES
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fderiv(ic,jc,kc) = zero
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE d2fdyz

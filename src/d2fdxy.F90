SUBROUTINE d2fdxy(functn,fderiv)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-05  Time: 11:29:02

!     *************************************************************************

!     D2FDXY
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     15-APR-2003:  RSC MODIFIED FOR SENGA2
!     10-OCT-2004:  RSC NULL VERSION

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES SECOND XY-DERIVATIVE OF SPECIFIED FUNCTION

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     ARGUMENTS
!     =========

REAL(KIND=dp), INTENT(IN OUT)         :: functn(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
REAL(KIND=dp), INTENT(OUT)            :: fderiv(nxsize,nysize,nzsize)




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
END SUBROUTINE d2fdxy

SUBROUTINE d2fdyz(functn,fderiv)

! Code converted using TO_F90 by Alan Miller
! Date: 2022-11-09  Time: 14:51:28

!     *************************************************************************

!     D2FDYZ
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     15-APR-2003:  RSC MODIFIED FOR SENGA2

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES SECOND YZ-DERIVATIVE OF SPECIFIED FUNCTION
!     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
!     EXPLICIT 8TH,6TH,4TH,4TH,4TH COMPATIBLE ORDER END CONDITIONS

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------

!     ARGUMENTS
!     =========

REAL(kind=8), INTENT(IN)             :: functn(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
REAL(kind=8), INTENT(OUT)            :: fderiv(nxsize,nysize,nzsize)




!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc

!     BEGIN
!     =====

!     =========================================================================

!     RSC 31-DEC-2006 BUG FIX INDICES
      do kc = kstal,kstol
        do jc = jstal,jstol
          do ic = istal,istol
             
            fderiv(ic,jc,kc) = zero

          enddo
        enddo
      enddo

!     =========================================================================


RETURN
END SUBROUTINE d2fdyz

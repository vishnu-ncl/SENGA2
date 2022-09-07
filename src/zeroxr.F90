SUBROUTINE zeroxr(farray)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-05  Time: 10:55:42

!     *************************************************************************

!     ZEROXR
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     06-JUL-2003:  RSC MODIFIED FOR SENGA2

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     ZEROS THE X-WISE RIGHT END ELEMENTS OF THE ARRAY

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!     -------------------------------------------------------------------------


!     ARGUMENTS
!     =========

DOUBLE PRECISION, INTENT(OUT)            :: farray(nxsize,nysize,nzsize)



!     LOCAL DATA
!     ==========
INTEGER :: jc,kc


!     BEGIN
!     =====

!     =========================================================================

DO kc = kstal,kstol
  DO jc = jstal,jstol
    
    farray(istol,jc,kc) = zero
    
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE zeroxr

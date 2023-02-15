SUBROUTINE zerozl(farray)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-05  Time: 10:55:51

!     *************************************************************************

!     ZEROZL
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
!     ZEROS THE Z-WISE LEFT END ELEMENTS OF THE ARRAY

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
INTEGER :: ic,jc


!     BEGIN
!     =====

!     =========================================================================

DO jc = jstal,jstol
  DO ic = istal,istol
    
    farray(ic,jc,kstal) = zero
    
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE zerozl

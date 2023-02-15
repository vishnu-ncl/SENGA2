SUBROUTINE zeroyr(farray)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-05  Time: 10:55:48

!     *************************************************************************

!     ZEROYR
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
!     ZEROS THE Y-WISE RIGHT END ELEMENTS OF THE ARRAY

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
INTEGER :: ic,kc


!     BEGIN
!     =====

!     =========================================================================

DO kc = kstal,kstol
  DO ic = istal,istol
    
    farray(ic,jstol,kc) = zero
    
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE zeroyr

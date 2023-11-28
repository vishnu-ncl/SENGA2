SUBROUTINE cxhalo(bigarr,jndexl,jndexr,kndexl,kndexr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:15:48

!     *************************************************************************

!     CXHALO
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPERTMENT

!     CHANGE RECORD
!     -------------
!     14-MAY-2003:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     CARRIES OUT HALO EXCHANGE FOR PERIODIC BCS IN X DIRECTION

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     ARGUMENTS
!     =========

REAL(kind=8), INTENT(OUT)            :: bigarr(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
INTEGER, INTENT(IN)                      :: jndexl
INTEGER, INTENT(IN)                      :: jndexr
INTEGER, INTENT(IN)                      :: kndexl
INTEGER, INTENT(IN)                      :: kndexr




!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc
INTEGER :: is


!     BEGIN
!     =====

!     =========================================================================

!     RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
DO kc = kndexl,kndexr
  DO jc = jndexl,jndexr
    
    is = istali - 1
    DO ic = istaro,istoro
      
      is = is + 1
      bigarr(ic,jc,kc) = bigarr(is,jc,kc)
      
    END DO
  END DO
END DO

!     =========================================================================

!     LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
DO kc = kndexl,kndexr
  DO jc = jndexl,jndexr
    
    is = istari - 1
    DO ic = istalo,istolo
      
      is = is + 1
      bigarr(ic,jc,kc) = bigarr(is,jc,kc)
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE cxhalo

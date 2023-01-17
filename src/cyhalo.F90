SUBROUTINE cyhalo(bigarr,indexl,indexr,kndexl,kndexr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:15:58

!     *************************************************************************

!     CYHALO
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPERTMENT

!     CHANGE RECORD
!     -------------
!     16-MAY-2003:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     CARRIES OUT HALO EXCHANGE FOR PERIODIC BCS IN Y DIRECTION

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     ARGUMENTS
!     =========

REAL(KIND=dp), INTENT(OUT)            :: bigarr(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
INTEGER, INTENT(IN)                      :: indexl
INTEGER, INTENT(IN)                      :: indexr
INTEGER, INTENT(IN)                      :: kndexl
INTEGER, INTENT(IN)                      :: kndexr




!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc
INTEGER :: js


!     BEGIN
!     =====

!     =========================================================================

!     RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
DO kc = kndexl,kndexr
  
  js = 0
  DO jc = nysize+1,nysize+nhaloy
    
    js = js + 1
    
    DO ic = indexl,indexr
      
      bigarr(ic,jc,kc) = bigarr(ic,js,kc)
      
    END DO
    
  END DO
  
END DO

!     =========================================================================

!     LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
DO kc = kndexl,kndexr
  
  js = (nysize+1-nhaloy) - 1
  DO jc = 1-nhaloy,0
    
    js = js + 1
    
    DO ic = indexl,indexr
      
      bigarr(ic,jc,kc) = bigarr(ic,js,kc)
      
    END DO
    
  END DO
  
END DO

!     =========================================================================


RETURN
END SUBROUTINE cyhalo

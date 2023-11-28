SUBROUTINE cyhals(bigarr,indexl,indexr,kndexl,kndexr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:16:04

!     *************************************************************************

!     CYHALS
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
!     FOR SPECIES MASS FRACTIONS

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     ARGUMENTS
!     =========

REAL(kind=8), INTENT(OUT)            :: bigarr(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr,nspcmx)
INTEGER, INTENT(IN)                      :: indexl
INTEGER, INTENT(IN)                      :: indexr
INTEGER, INTENT(IN)                      :: kndexl
INTEGER, INTENT(IN)                      :: kndexr




!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc
INTEGER :: js,ispec


!     BEGIN
!     =====

!     =========================================================================

!     RUN THROUGH ALL SPECIES
!     -----------------------
DO ispec = 1,nspec
  
!       =======================================================================
  
!       RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
  DO kc = kndexl,kndexr
    
    js = jstali - 1
    DO jc = jstaro,jstoro
      
      js = js + 1
      
      DO ic = indexl,indexr
        
        bigarr(ic,jc,kc,ispec) = bigarr(ic,js,kc,ispec)
        
      END DO
      
    END DO
    
  END DO
  
!       =======================================================================
  
!       LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
  DO kc = kndexl,kndexr
    
    js = jstari - 1
    DO jc = jstalo,jstolo
      
      js = js + 1
      
      DO ic = indexl,indexr
        
        bigarr(ic,jc,kc,ispec) = bigarr(ic,js,kc,ispec)
        
      END DO
      
    END DO
    
  END DO
  
!       =======================================================================
  
END DO
!     END OF SPECIES LOOP

!     =========================================================================


RETURN
END SUBROUTINE cyhals

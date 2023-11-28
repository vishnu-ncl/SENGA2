SUBROUTINE czhals(bigarr,indexl,indexr,jndexl,jndexr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:16:17

!     *************************************************************************

!     CZHALS
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
!     CARRIES OUT HALO EXCHANGE FOR PERIODIC BCS IN Z DIRECTION
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
INTEGER, INTENT(IN)                      :: jndexl
INTEGER, INTENT(IN)                      :: jndexr




!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc
INTEGER :: ks,ispec


!     BEGIN
!     =====

!     =========================================================================

!     RUN THROUGH ALL SPECIES
!     -----------------------
DO ispec = 1,nspec
  
!       =======================================================================
  
!       RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
  ks = kstali - 1
  DO kc = kstaro,kstoro
    
    ks = ks + 1
    
    DO jc = jndexl,jndexr
      DO ic = indexl,indexr
        
        bigarr(ic,jc,kc,ispec) = bigarr(ic,jc,ks,ispec)
        
      END DO
    END DO
    
  END DO
  
!       =======================================================================
  
!       LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
  ks = kstari - 1
  DO kc = kstalo,kstolo
    
    ks = ks + 1
    
    DO jc = jndexl,jndexr
      DO ic = indexl,indexr
        
        bigarr(ic,jc,kc,ispec) = bigarr(ic,jc,ks,ispec)
        
      END DO
    END DO
    
  END DO
  
!       =======================================================================
  
END DO
!     END OF SPECIES LOOP

!     =========================================================================


RETURN
END SUBROUTINE czhals

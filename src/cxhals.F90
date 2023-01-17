SUBROUTINE cxhals(bigarr,jndexl,jndexr,kndexl,kndexr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:15:53

!     *************************************************************************

!     CXHALS
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
!     CARRIES OUT HALO EXCHANGE FOR PERIODIC BCS IN X DIRECTION
!     FOR SPECIES MASS FRACTIONS

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     ARGUMENTS
!     =========

REAL(KIND=dp), INTENT(OUT)            :: bigarr(nspcmx,nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
INTEGER, INTENT(IN)                      :: jndexl
INTEGER, INTENT(IN)                      :: jndexr
INTEGER, INTENT(IN)                      :: kndexl
INTEGER, INTENT(IN)                      :: kndexr




!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc
INTEGER :: is,ispec


!     BEGIN
!     =====

!     =========================================================================

!     RUN THROUGH ALL SPECIES
!     -----------------------
DO ispec = 1,nspec
  
!       =======================================================================
  
!       RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
  DO kc = kndexl,kndexr
    DO jc = jndexl,jndexr
      
      is = 0
      DO ic = nxsize+1,nxsize+nhalox
        
        is = is + 1
        bigarr(ispec,ic,jc,kc) = bigarr(ispec,is,jc,kc)
        
      END DO
    END DO
  END DO
  
!       =======================================================================
  
!       LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
  DO kc = kndexl,kndexr
    DO jc = jndexl,jndexr
      
      is = (nxsize+1-nhalox) - 1
      DO ic = 1-nhalox,0
        
        is = is + 1
        bigarr(ispec,ic,jc,kc) = bigarr(ispec,is,jc,kc)
        
      END DO
    END DO
  END DO
  
!       =======================================================================
  
END DO
!     END OF SPECIES LOOP

!     =========================================================================


RETURN
END SUBROUTINE cxhals

SUBROUTINE radcal
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:26:10

!     *************************************************************************

!     RADCAL
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     14-JUL-2013:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     RADIATION TREATMENT
!     USING OPTICALLY THIN ASSUMPTION: Ju et al: JFM 342, 315-334, 1997.
!     AFTER TOM DUNSTAN 2012

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
real(kind=8):: plspec,fornow
INTEGER :: ic,jc,kc,ispec,jspec,icp


!     BEGIN
!     =====

!     =========================================================================

!     BUILD THE PLANCK MEAN ABSORPTION COEFFICIENT OF THE MIXTURE
!     -----------------------------------------------------------

!     INITIALISE THE ACCUMULATOR
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store1(ic,jc,kc) = zero
      
    END DO
  END DO
END DO

!     -------------------------------------------------------------------------

!     RUN THROUGH ALL RADIATING SPECIES
DO jspec = 1, nsprad
  
!       PLANCK MEAN ABSORPTION COEFFICIENT OF EACH SPECIES
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        fornow = trun(ic,jc,kc)
        plspec = akprad(nkprad(jspec),jspec)
        DO icp = nkprm1(jspec),1,-1
          plspec = plspec*fornow + akprad(icp,jspec)
        END DO
        store2(ic,jc,kc) = plspec
        
      END DO
    END DO
  END DO
  
!       SPECIES ID
  ispec = nsprid(jspec)
  
!       ADD THE SPECIES CONTRIBUTION
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        fornow = yrhs(ic,jc,kc,ispec)*rgspec(ispec)*trun(ic,jc,kc)
        store1(ic,jc,kc) = store1(ic,jc,kc) + store2(ic,jc,kc)*fornow
        
      END DO
    END DO
  END DO
  
END DO

!     =========================================================================

!     INCLUDE THE RADIATION TERM IN THE ENERGY EQUATION

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fornow = trun(ic,jc,kc)
      fornow = fornow*fornow*fornow*fornow
      
      erhs(ic,jc,kc) = erhs(ic,jc,kc)  &
          - foursb*store1(ic,jc,kc)*(fornow - trfrth)
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE radcal

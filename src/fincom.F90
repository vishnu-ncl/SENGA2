SUBROUTINE fincom
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:25:43

!     *************************************************************************

!     FINCOM
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     15-JAN-2003:  CREATED
!     08-AUG-2012:  RSC EVALUATE ALL SPECIES

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     COMPUTES FINAL SOLUTION VALUES IN ERK SCHEME
!     BY DOING A LINEAR COMBINATION OF LEFT- AND RIGHT-HAND SIDES

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc,ispec

!     -------------------------------------------------------------------------

!     BEGIN
!     =====

!     =========================================================================

!     FINAL ERK SUBSTEP
!     =================

!     -------------------------------------------------------------------------
!     NOTE: ALL ERK ERROR ARRAYS ARE INITIALISED TO ZERO IN SUBROUTINE ADAPTT
!     -------------------------------------------------------------------------

!     DENSITY
!     ----------
DO kc = kstald,kstold
  DO jc = jstald,jstold
    DO ic = istald,istold
      
      derr(ic,jc,kc) = derr(ic,jc,kc) + rkerr(nrkstp)*drhs(ic,jc,kc)
      
      drun(ic,jc,kc) = drun(ic,jc,kc) + rklhs(nrkstp)*drhs(ic,jc,kc)
      drhs(ic,jc,kc) = drun(ic,jc,kc)
      
    END DO
  END DO
END DO

!     -------------------------------------------------------------------------

!     U VELOCITY
!     ----------
DO kc = kstalu,kstolu
  DO jc = jstalu,jstolu
    DO ic = istalu,istolu
      
      uerr(ic,jc,kc) = uerr(ic,jc,kc) + rkerr(nrkstp)*urhs(ic,jc,kc)
      
      urun(ic,jc,kc) = urun(ic,jc,kc) + rklhs(nrkstp)*urhs(ic,jc,kc)
      urhs(ic,jc,kc) = urun(ic,jc,kc)
      
    END DO
  END DO
END DO

!     -------------------------------------------------------------------------

!     V-VELOCITY
!     ----------
DO kc = kstalv,kstolv
  DO jc = jstalv,jstolv
    DO ic = istalv,istolv
      
      verr(ic,jc,kc) = verr(ic,jc,kc) + rkerr(nrkstp)*vrhs(ic,jc,kc)
      
      vrun(ic,jc,kc) = vrun(ic,jc,kc) + rklhs(nrkstp)*vrhs(ic,jc,kc)
      vrhs(ic,jc,kc) = vrun(ic,jc,kc)
      
    END DO
  END DO
END DO

!     -------------------------------------------------------------------------

!     W-VELOCITY
!     ----------
DO kc = kstalw,kstolw
  DO jc = jstalw,jstolw
    DO ic = istalw,istolw
      
      werr(ic,jc,kc) = werr(ic,jc,kc) + rkerr(nrkstp)*wrhs(ic,jc,kc)
      
      wrun(ic,jc,kc) = wrun(ic,jc,kc) + rklhs(nrkstp)*wrhs(ic,jc,kc)
      wrhs(ic,jc,kc) = wrun(ic,jc,kc)
      
    END DO
  END DO
END DO

!     -------------------------------------------------------------------------

!     STAGNATION INTERNAL ENERGY
!     --------------------------
DO kc = kstale,kstole
  DO jc = jstale,jstole
    DO ic = istale,istole
      
      eerr(ic,jc,kc) = eerr(ic,jc,kc) + rkerr(nrkstp)*erhs(ic,jc,kc)
      
      erun(ic,jc,kc) = erun(ic,jc,kc) + rklhs(nrkstp)*erhs(ic,jc,kc)
      erhs(ic,jc,kc) = erun(ic,jc,kc)
      
    END DO
  END DO
END DO

!     -------------------------------------------------------------------------

!     SPECIES MASS FRACTIONS
!     ----------------------
!     RSC 08-AUG-2012 EVALUATE ALL SPECIES
!      DO ISPEC = 1,NSPM1
DO ispec = 1,nspec
  
  DO kc = kstaly,kstoly
    DO jc = jstaly,jstoly
      DO ic = istaly,istoly
        
        yerr(ic,jc,kc,ispec) = yerr(ic,jc,kc,ispec)  &
            + rkerr(nrkstp)*yrhs(ic,jc,kc,ispec)
        
        yrun(ic,jc,kc,ispec) = yrun(ic,jc,kc,ispec)  &
            + rklhs(nrkstp)*yrhs(ic,jc,kc,ispec)
        yrhs(ic,jc,kc,ispec) = yrun(ic,jc,kc,ispec)
        
      END DO
    END DO
  END DO
  
END DO

!     -------------------------------------------------------------------------

!C     NTH SPECIES
!      DO KC = KSTALY,KSTOLY
!        DO JC = JSTALY,JSTOLY
!          DO IC = ISTALY,ISTOLY

!            YRUN(IC,JC,KC,NSPEC) = ZERO

!          ENDDO
!        ENDDO
!      ENDDO

!      DO ISPEC = 1,NSPM1
!        DO KC = KSTALY,KSTOLY
!          DO JC = JSTALY,JSTOLY
!            DO IC = ISTALY,ISTOLY

!              YRUN(IC,JC,KC,NSPEC) = YRUN(IC,JC,KC,NSPEC)
!     +                             + YRUN(IC,JC,KC,ISPEC)

!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO

!      DO KC = KSTALY,KSTOLY
!        DO JC = JSTALY,JSTOLY
!          DO IC = ISTALY,ISTOLY

!            YRUN(IC,JC,KC,NSPEC)
!     +        = DRUN(IC,JC,KC)*(ONE-YRUN(IC,JC,KC,NSPEC)/DRUN(IC,JC,KC))

!            YRHS(IC,JC,KC,NSPEC) = YRUN(IC,JC,KC,NSPEC)

!          ENDDO
!        ENDDO
!      ENDDO

!     =========================================================================


RETURN
END SUBROUTINE fincom

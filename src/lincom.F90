SUBROUTINE lincom
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:26:06

!     *************************************************************************

!     LINCOM
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
!     COMPUTES INTERMEDIATE SOLUTION VALUES IN ERK SCHEME
!     BY DOING LINEAR COMBINATIONS OF LEFT- AND RIGHT-HAND SIDES

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
real(kind=dp) :: fornow
INTEGER :: ic,jc,kc,ispec


!     BEGIN
!     =====

!     =========================================================================

!     ERK SUBSTEP
!     ===========

!     -------------------------------------------------------------------------
!     NOTE: ALL ERK ERROR ARRAYS ARE INITIALISED TO ZERO IN SUBROUTINE ADAPTT
!     -------------------------------------------------------------------------

!     DENSITY
!     -------
DO kc = kstald,kstold
  DO jc = jstald,jstold
    DO ic = istald,istold
      
      derr(ic,jc,kc) = derr(ic,jc,kc) + rkerr(irkstp)*drhs(ic,jc,kc)
      
      fornow = drun(ic,jc,kc)
      drun(ic,jc,kc) = fornow + rklhs(irkstp)*drhs(ic,jc,kc)
      drhs(ic,jc,kc) = fornow + rkrhs(irkstp)*drhs(ic,jc,kc)
      
    END DO
  END DO
END DO

!     -------------------------------------------------------------------------

!     U-VELOCITY
!     ----------
DO kc = kstalu,kstolu
  DO jc = jstalu,jstolu
    DO ic = istalu,istolu
      
      uerr(ic,jc,kc) = uerr(ic,jc,kc) + rkerr(irkstp)*urhs(ic,jc,kc)
      
      fornow = urun(ic,jc,kc)
      urun(ic,jc,kc) = fornow + rklhs(irkstp)*urhs(ic,jc,kc)
      urhs(ic,jc,kc) = fornow + rkrhs(irkstp)*urhs(ic,jc,kc)
      
    END DO
  END DO
END DO

!     -------------------------------------------------------------------------

!     V-VELOCITY
!     ----------
DO kc = kstalv,kstolv
  DO jc = jstalv,jstolv
    DO ic = istalv,istolv
      
      verr(ic,jc,kc) = verr(ic,jc,kc) + rkerr(irkstp)*vrhs(ic,jc,kc)
      
      fornow = vrun(ic,jc,kc)
      vrun(ic,jc,kc) = fornow + rklhs(irkstp)*vrhs(ic,jc,kc)
      vrhs(ic,jc,kc) = fornow + rkrhs(irkstp)*vrhs(ic,jc,kc)
      
    END DO
  END DO
END DO

!     -------------------------------------------------------------------------

!     W-VELOCITY
!     ----------
DO kc = kstalw,kstolw
  DO jc = jstalw,jstolw
    DO ic = istalw,istolw
      
      werr(ic,jc,kc) = werr(ic,jc,kc) + rkerr(irkstp)*wrhs(ic,jc,kc)
      
      fornow = wrun(ic,jc,kc)
      wrun(ic,jc,kc) = fornow + rklhs(irkstp)*wrhs(ic,jc,kc)
      wrhs(ic,jc,kc) = fornow + rkrhs(irkstp)*wrhs(ic,jc,kc)
      
    END DO
  END DO
END DO

!     -------------------------------------------------------------------------

!     STAGNATION INTERNAL ENERGY
!     --------------------------
DO kc = kstale,kstole
  DO jc = jstale,jstole
    DO ic = istale,istole
      
      eerr(ic,jc,kc) = eerr(ic,jc,kc) + rkerr(irkstp)*erhs(ic,jc,kc)
      
      fornow = erun(ic,jc,kc)
      erun(ic,jc,kc) = fornow + rklhs(irkstp)*erhs(ic,jc,kc)
      erhs(ic,jc,kc) = fornow + rkrhs(irkstp)*erhs(ic,jc,kc)
      
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
            + rkerr(irkstp)*yrhs(ic,jc,kc,ispec)
        
        fornow = yrun(ic,jc,kc,ispec)
        yrun(ic,jc,kc,ispec) = fornow + rklhs(irkstp)*yrhs(ic,jc,kc,ispec)
        yrhs(ic,jc,kc,ispec) = fornow + rkrhs(irkstp)*yrhs(ic,jc,kc,ispec)
        
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
!            YRHS(IC,JC,KC,NSPEC) = ZERO

!          ENDDO
!        ENDDO
!      ENDDO

!      DO ISPEC = 1,NSPM1
!        DO KC = KSTALY,KSTOLY
!          DO JC = JSTALY,JSTOLY
!            DO IC = ISTALY,ISTOLY

!              YRUN(IC,JC,KC,NSPEC) = YRUN(IC,JC,KC,NSPEC)
!     +                             + YRUN(IC,JC,KC,ISPEC)
!              YRHS(IC,JC,KC,NSPEC) = YRHS(IC,JC,KC,NSPEC)
!     +                             + YRHS(IC,JC,KC,ISPEC)

!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO

!      DO KC = KSTALY,KSTOLY
!        DO JC = JSTALY,JSTOLY
!          DO IC = ISTALY,ISTOLY

!            YRUN(IC,JC,KC,NSPEC)
!     +        = DRUN(IC,JC,KC)*(ONE-YRUN(IC,JC,KC,NSPEC)/DRUN(IC,JC,KC))

!            YRHS(IC,JC,KC,NSPEC)
!     +        = DRHS(IC,JC,KC)*(ONE-YRHS(IC,JC,KC,NSPEC)/DRHS(IC,JC,KC))

!          ENDDO
!        ENDDO
!      ENDDO

!     =========================================================================


RETURN
END SUBROUTINE lincom

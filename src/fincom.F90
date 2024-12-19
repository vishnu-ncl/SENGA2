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
DOUBLE PRECISION :: temp1,temp2,temp3,temp4
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

!     VM & NC: GRADIENT OF SPECIES AT WALL EQUAL TO ZERO
IF((nsbcxl == nsbcw2).OR.(nsbcxl == nsbcw1))THEN
  DO ispec=1,nspec
    DO kc=kstal,kstol
      DO jc=jstal,jstol
        temp1=yrhs(istal+1,jc,kc,ispec)/drhs(istal+1,jc,kc)
        temp2=yrhs(istal+2,jc,kc,ispec)/drhs(istal+2,jc,kc)
        temp3=yrhs(istal+3,jc,kc,ispec)/drhs(istal+3,jc,kc)
        temp4=yrhs(istal+4,jc,kc,ispec)/drhs(istal+4,jc,kc)
        yrun(istal,jc,kc,ispec)=(12.0/25.0)*(4.0*temp1-3.0*temp2+  &
            (4.0/3.0)*temp3-(1.0/4.0*temp4))
        yrun(istal,jc,kc,ispec)=MIN(1.0,yrun(istal,jc,kc,ispec))
        yrun(istal,jc,kc,ispec)=MAX(0.0,yrun(istal,jc,kc,ispec))
        yrun(istal,jc,kc,ispec)=drhs(istal,jc,kc) *yrun(istal,jc,kc,ispec)
        yrhs(istal,jc,kc,ispec)=yrun(istal,jc,kc,ispec)
      END DO
    END DO
  END DO
END IF

IF((nsbcxr == nsbcw2).OR.(nsbcxr == nsbcw1))THEN
  DO ispec=1,nspec
    DO kc=kstal,kstol
      DO jc=jstal,jstol
        temp1=yrhs(istol-1,jc,kc,ispec)/drhs(istol-1,jc,kc)
        temp2=yrhs(istol-2,jc,kc,ispec)/drhs(istol-2,jc,kc)
        temp3=yrhs(istol-3,jc,kc,ispec)/drhs(istol-3,jc,kc)
        temp4=yrhs(istol-4,jc,kc,ispec)/drhs(istol-4,jc,kc)
        yrun(istol,jc,kc,ispec)=(12.0/25.0)*(4.0*temp1-3.0*temp2+  &
            (4.0/3.0)*temp3-(1.0/4.0*temp4))
        yrun(istol,jc,kc,ispec)=MIN(1.0,yrun(istol,jc,kc,ispec))
        yrun(istol,jc,kc,ispec)=MAX(0.0,yrun(istol,jc,kc,ispec))
        yrun(istol,jc,kc,ispec)=drhs(istol,jc,kc) *yrun(istol,jc,kc,ispec)
        yrhs(istol,jc,kc,ispec)=yrun(istol,jc,kc,ispec)
      END DO
    END DO
  END DO
END IF

IF((nsbcyl == nsbcw2).OR.(nsbcyl == nsbcw1))THEN
  DO ispec=1,nspec
    DO kc=kstal,kstol
      DO ic=istal,istol
        temp1=yrhs(ic,jstal+1,kc,ispec)/drhs(ic,jstal+1,kc)
        temp2=yrhs(ic,jstal+2,kc,ispec)/drhs(ic,jstal+2,kc)
        temp3=yrhs(ic,jstal+3,kc,ispec)/drhs(ic,jstal+3,kc)
        temp4=yrhs(ic,jstal+4,kc,ispec)/drhs(ic,jstal+4,kc)
        yrun(ic,jstal,kc,ispec)=(12.0/25.0)*(4.0*temp1-3.0*temp2+  &
            (4.0/3.0)*temp3-(1.0/4.0*temp4))
        yrun(ic,jstal,kc,ispec)=MIN(1.0,yrun(ic,jstal,kc,ispec))
        yrun(ic,jstal,kc,ispec)=MAX(0.0,yrun(ic,jstal,kc,ispec))
        yrun(ic,jstal,kc,ispec)=drhs(ic,jstal,kc) *yrun(ic,jstal,kc,ispec)
        yrhs(ic,jstal,kc,ispec)=yrun(ic,jstal,kc,ispec)
      END DO
    END DO
  END DO
END IF

IF((nsbcyr == nsbcw2).OR.(nsbcyr == nsbcw1))THEN
  DO ispec=1,nspec
    DO kc=kstal,kstol
      DO ic=istal,istol
        temp1=yrhs(ic,jstol-1,kc,ispec)/drhs(ic,jstol-1,kc)
        temp2=yrhs(ic,jstol-2,kc,ispec)/drhs(ic,jstol-2,kc)
        temp3=yrhs(ic,jstol-3,kc,ispec)/drhs(ic,jstol-3,kc)
        temp4=yrhs(ic,jstol-4,kc,ispec)/drhs(ic,jstol-4,kc)
        yrun(ic,jstol,kc,ispec)=(12.0/25.0)*(4.0*temp1-3.0*temp2+  &
            (4.0/3.0)*temp3-(1.0/4.0*temp4))
        yrun(ic,jstol,kc,ispec)=MIN(1.0,yrun(ic,jstol,kc,ispec))
        yrun(ic,jstol,kc,ispec)=MAX(0.0,yrun(ic,jstol,kc,ispec))
        yrun(ic,jstol,kc,ispec)=drhs(ic,jstol,kc) *yrun(ic,jstol,kc,ispec)
        yrhs(ic,jstol,kc,ispec)=yrun(ic,jstol,kc,ispec)
      END DO
    END DO
  END DO
END IF

IF((nsbczl == nsbcw2).OR.(nsbczl == nsbcw1))THEN
  DO ispec=1,nspec
    DO jc=jstal,jstol
      DO ic=istal,istol
        temp1=yrhs(ic,jc,kstal+1,ispec)/drhs(ic,jc,kstal+1)
        temp2=yrhs(ic,jc,kstal+2,ispec)/drhs(ic,jc,kstal+2)
        temp3=yrhs(ic,jc,kstal+3,ispec)/drhs(ic,jc,kstal+3)
        temp4=yrhs(ic,jc,kstal+4,ispec)/drhs(ic,jc,kstal+4)
        yrun(ic,jc,kstal,ispec)=(12.0/25.0)*(4.0*temp1-3.0*temp2+  &
            (4.0/3.0)*temp3-(1.0/4.0*temp4))
        yrun(ic,jc,kstal,ispec)=MIN(1.0,yrun(ic,jc,kstal,ispec))
        yrun(ic,jc,kstal,ispec)=MAX(0.0,yrun(ic,jc,kstal,ispec))
        yrun(ic,jc,kstal,ispec)=drhs(ic,jc,kstal) *yrun(ic,jc,kstal,ispec)
        yrhs(ic,jc,kstal,ispec)=yrun(ic,jc,kstal,ispec)
      END DO
    END DO
  END DO
END IF

IF((nsbczr == nsbcw2).OR.(nsbczr == nsbcw1))THEN
  DO ispec=1,nspec
    DO jc=jstal,jstol
      DO ic=istal,istol
        temp1=yrhs(ic,jc,kstol-1,ispec)/drhs(ic,jc,kstol-1)
        temp2=yrhs(ic,jc,kstol-2,ispec)/drhs(ic,jc,kstol-2)
        temp3=yrhs(ic,jc,kstol-3,ispec)/drhs(ic,jc,kstol-3)
        temp4=yrhs(ic,jc,kstol-4,ispec)/drhs(ic,jc,kstol-4)
        yrun(ic,jc,kstol,ispec)=(12.0/25.0)*(4.0*temp1-3.0*temp2+  &
            (4.0/3.0)*temp3-(1.0/4.0*temp4))
        yrun(ic,jc,kstol,ispec)=MIN(1.0,yrun(ic,jc,kstol,ispec))
        yrun(ic,jc,kstol,ispec)=MAX(0.0,yrun(ic,jc,kstol,ispec))
        yrun(ic,jc,kstol,ispec)=drhs(ic,jc,kstol) *yrun(ic,jc,kstol,ispec)
        yrhs(ic,jc,kstol,ispec)=yrun(ic,jc,kstol,ispec)
      END DO
    END DO
  END DO
END IF

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

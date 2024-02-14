SUBROUTINE flamin
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:25:46

!     *************************************************************************

!     FLAMIN
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     28-DEC-2003:  CREATED
!     08-JAN-2005:  RSC INITIAL 1D LAMINAR FLAME PROFILE

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     SETS INITIAL THERMOCHEMICAL FIELD
!     1D LAMINAR FLAME PROFILE (LEFT OR RIGHT FACING)
!     SPECIAL FOR 21 STEP HYDROGEN MECHAMISM

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!     -------------------------------------------------------------------------


!     PARAMETERS
!     ==========
!     ESTIMATED FLAME LOCATION AND THICKNESS
real(kind=8) :: clocat,cthick
PARAMETER(clocat = 0.0025_8, cthick = 0.0005_8)

!C     PINCH OF HYDROGEN ATOM
!      real(kind=8) HPINCH,HLOCAT,HTHICK
!      PARAMETER(HPINCH = 1.0D-10, HLOCAT = 2.5D-3, HTHICK = 1.0D-4)
!C     PINCH OF HYDROGEN MOLECULE
!      real(kind=8) H2PNCH,H2LOCT,H2THCK
!      PARAMETER(H2PNCH = 1.0D-6, H2LOCT = 2.5D-3, H2THCK = 2.5D-4)


!     FUNCTION
!     ========
real(kind=8) :: erfunc
EXTERNAL erfunc


!     LOCAL DATA
!     ==========
real(kind=8) :: crin(1:nxsize)
real(kind=8) :: yrinr(nspcmx),yrinp(nspcmx)
real(kind=8) :: trinr,trinp
real(kind=8) :: deltag,xcoord,argmnt
real(kind=8) :: flxmas
INTEGER :: icproc
INTEGER :: igofst
INTEGER :: ix
INTEGER :: ic,jc,kc
INTEGER :: ispec


!     BEGIN
!     =====

!     =========================================================================

!     SPECIFY INITIAL THERMOCHEMICAL FIELD HERE
!     =========================================


!     SET PRODUCT TEMPERATURE
!     -----------------------
!     REACTANT TEMPERATURE SET IN CONTROL FILE
trinr = trin
!      TRINP = 2330.96554
trinp = 2200.0_8


!     SET SPECIES MASS FRACTIONS
!     --------------------------
!     OVERRIDE MASS FRACTION VALUES SET IN CONTROL FILE

!     REACTANTS
yrinr(1) = 0.0199886_8     !2.8312571D-2
yrinr(2) = 0.2286239_8     !2.26500566D-1
DO ispec = 3,nspm1
  yrinr(ispec) = zero
END DO

yrinr(nspec) = zero
DO ispec = 1,nspm1
  yrinr(nspec) = yrinr(nspec) + yrinr(ispec)
END DO
yrinr(nspec) = one - yrinr(nspec)

!     PRODUCTS
yrinp(1) = zero
yrinp(2) = 0.0685323_8     !ZERO
yrinp(3) = 0.1798974_8     !2.54716981D-1
DO ispec = 4,nspm1
  yrinp(ispec) = zero
END DO
yrinp(nspec) = zero
DO ispec = 1,nspm1
  yrinp(nspec) = yrinp(nspec) + yrinp(ispec)
END DO
yrinp(nspec) = one - yrinp(nspec)

!     WRITE TO REPORT FILE
IF(iproc == 0)THEN
  
!        OPEN(UNIT=NCREPT,FILE=FNREPT,STATUS='OLD',FORM='FORMATTED')
  
!C       GO TO EOF
!1000    CONTINUE
!          READ(NCREPT,9000,END=1010)
!          GOTO 1000
!1010    BACKSPACE(NCREPT)
  
  WRITE(ncrept,*)
  WRITE(ncrept,*)'FLAMIN: reactant mass fractions:'
  DO ispec = 1,nspec
    WRITE(ncrept,'(I5,1PE15.7)')ispec,yrinr(ispec)
  END DO
  WRITE(ncrept,*)
  
  WRITE(ncrept,*)'FLAMIN: product mass fractions:'
  DO ispec = 1,nspec
    WRITE(ncrept,'(I5,1PE15.7)')ispec,yrinp(ispec)
  END DO
  WRITE(ncrept,*)
  
  WRITE(ncrept,*)'FLAMIN: reactant and product temperatures:'
  WRITE(ncrept,'(2(1PE15.7))')trinr,trinp
  WRITE(ncrept,*)
  
!        CLOSE(NCREPT)
  
END IF


!     GLOBAL INDEXING
!     ---------------
deltag = xgdlen/(REAL(nxglbl-1,kind=8))

igofst = 0
DO icproc = 0, ixproc-1
  igofst = igofst + npmapx(icproc)
END DO


!     SET REACTION PROGRESS VARIABLE PROFILE
!     --------------------------------------
!     SIMPLE 1D LEFT-FACING ERROR FUNCTION PROFILE
DO ic = istal,istol
  
  ix = igofst + ic
  xcoord = REAL(ix-1,kind=8)*deltag
  argmnt = (xcoord-clocat)/cthick
  crin(ic) = half*(one+erfunc(argmnt))
  
END DO

!C     SIMPLE 1D RIGHT-FACING ERROR FUNCTION PROFILE
!      DO IC = ISTAL,ISTOL

!        IX = IGOFST + IC
!        XCOORD = REAL(IX-1)*DELTAG
!        ARGMNT = (XCOORD-CLOCAT)/CTHICK
!        CRIN(IC) = HALF*(ONE+ERFUNC(-ARGMNT))

!      ENDDO


!     SET SPECIES MASS FRACTION PROFILES
!     ----------------------------------
DO ispec = 1, nspm1
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        yrun(ic,jc,kc,ispec) = yrinr(ispec)  &
            + crin(ic)*(yrinp(ispec) - yrinr(ispec))
        
      END DO
    END DO
  END DO
  
END DO

!C     SG 25-STEP MECHANISM
!C     ADD A PINCH OF HYDROGEN ATOM
!C     AND HYDROGEN MOLECULE
!      DO KC = KSTAL,KSTOL
!        DO JC = JSTAL,JSTOL
!          DO IC = ISTAL,ISTOL

!            IX = IGOFST + IC
!            XCOORD = REAL(IX-1)*DELTAG
!            ARGMNT = (XCOORD-H2LOCT)/H2THCK
!            YRUN(IC,JC,KC,5) = H2PNCH*EXP(-ARGMNT*ARGMNT)
!            ARGMNT = (XCOORD-HLOCAT)/HTHICK
!            YRUN(IC,JC,KC,8) = HPINCH*EXP(-ARGMNT*ARGMNT)

!          ENDDO
!        ENDDO
!      ENDDO

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      yrun(ic,jc,kc,nspec) = zero
      
    END DO
  END DO
END DO

DO ispec = 1, nspm1
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        yrun(ic,jc,kc,nspec) = yrun(ic,jc,kc,nspec) + yrun(ic,jc,kc,ispec)
        
      END DO
    END DO
  END DO
END DO

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      yrun(ic,jc,kc,nspec) = one - yrun(ic,jc,kc,nspec)
      
    END DO
  END DO
END DO


!     SET TEMPERATURE PROFILE
!     -----------------------
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      trun(ic,jc,kc) = trinr +  crin(ic)*(trinp - trinr)
      
    END DO
  END DO
END DO


!     SET DENSITY PROFILE ASSUMING CONSTANT PRESSURE
!     -------------------
!     PRESSURE SET IN CONTROL FILE

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      store1(ic,jc,kc) = zero
      
    END DO
  END DO
END DO

DO ispec = 1,nspec
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        store1(ic,jc,kc) = store1(ic,jc,kc)  &
            + rgspec(ispec)*yrun(ic,jc,kc,ispec)
        
      END DO
    END DO
  END DO
END DO

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      drun(ic,jc,kc) = prin/(store1(ic,jc,kc)*trun(ic,jc,kc))
      
    END DO
  END DO
END DO


!     SET VELOCITY PROFILE ASSUMING CONSTANT MASS FLUX
!     --------------------
!     INITIAL (INLET) VEOCITY SET IN CONTROL FILE
flxmas = drin*urin
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      urun(ic,jc,kc) = flxmas/drun(ic,jc,kc)
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN

9000  FORMAT(a)

END SUBROUTINE flamin

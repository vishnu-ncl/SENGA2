SUBROUTINE flamin

!   *************************************************************************

!   FLAMIN
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   28-DEC-2003:  CREATED
!   08-JAN-2005:  RSC INITIAL 1D LAMINAR FLAME PROFILE

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   SETS INITIAL THERMOCHEMICAL FIELD
!   1D LAMINAR FLAME PROFILE (LEFT OR RIGHT FACING)
!   SPECIAL FOR 21 STEP HYDROGEN MECHAMISM

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   PARAMETERS
!   ==========
!   ESTIMATED FLAME LOCATION AND THICKNESS
    DOUBLE PRECISION :: clocat,cthick
    PARAMETER(clocat = 2.5D-3, cthick = 5.0D-4)

!   PINCH OF HYDROGEN ATOM
!   DOUBLE PRECISION HPINCH,HLOCAT,HTHICK
!   PARAMETER(HPINCH = 1.0D-10, HLOCAT = 2.5D-3, HTHICK = 1.0D-4)
!   PINCH OF HYDROGEN MOLECULE
!   DOUBLE PRECISION H2PNCH,H2LOCT,H2THCK
!   PARAMETER(H2PNCH = 1.0D-6, H2LOCT = 2.5D-3, H2THCK = 2.5D-4)

!   LOCAL DATA
!   ==========
DOUBLE PRECISION :: crin(1:nxsize)
DOUBLE PRECISION :: yrinr(nspcmx),yrinp(nspcmx)
DOUBLE PRECISION :: trinr,trinp
DOUBLE PRECISION :: deltag,xcoord,argmnt
DOUBLE PRECISION :: flxmas
DOUBLE PRECISION :: drunr(nxglbl,nyglbl,nzglbl)
DOUBLE PRECISION :: urunr(nxglbl,nyglbl,nzglbl)
DOUBLE PRECISION :: trunr(nxglbl,nyglbl,nzglbl)
DOUBLE PRECISION :: yrunr(nxglbl,nyglbl,nzglbl,nspec)
INTEGER :: icproc
INTEGER :: igofst
INTEGER :: ix
INTEGER :: ic,jc,kc
INTEGER :: ispec


!   BEGIN
!   =====

!   =========================================================================

    num = 111
    OPEN(UNIT=num,FILE='input/h2_init.dat',FORM="FORMATTED",STATUS= "OLD")
    DO i=1,nxglbl
        READ(num,*)ix,urunr(i,1,1),trunr(i,1,1),drunr(i,1,1),  &
            yrunr(i,1,1,1),yrunr(i,1,1,2),yrunr(i,1,1,3),yrunr(i,1,1,4),  &
            yrunr(i,1,1,5),yrunr(i,1,1,6),yrunr(i,1,1,7),yrunr(i,1,1,8), yrunr(i,1,1,9)
    END DO
    CLOSE(num)

    igofst = 0
    DO icproc = 0, ixproc-1
        igofst = igofst + npmapx(icproc)
    END DO

    DO ic = istal,istol
        ix = igofst + ic
        urun(ic,1,1) = urunr(ix,1,1)
        vrun(ic,1,1) = 0.0D0
        wrun(ic,1,1) = 0.0D0
        drun(ic,1,1) = drunr(ix,1,1)
        trun(ic,1,1) = trunr(ix,1,1)
        yrun(ic,1,1,1) = yrunr(ix,1,1,1)
        yrun(ic,1,1,2) = yrunr(ix,1,1,2)
        yrun(ic,1,1,3) = yrunr(ix,1,1,3)
        yrun(ic,1,1,4) = yrunr(ix,1,1,4)
        yrun(ic,1,1,5) = yrunr(ix,1,1,5)
        yrun(ic,1,1,6) = yrunr(ix,1,1,6)
        yrun(ic,1,1,7) = yrunr(ix,1,1,7)
        yrun(ic,1,1,8) = yrunr(ix,1,1,8)
        yrun(ic,1,1,9) = yrunr(ix,1,1,9)
    END DO

!   =========================================================================

END SUBROUTINE flamin

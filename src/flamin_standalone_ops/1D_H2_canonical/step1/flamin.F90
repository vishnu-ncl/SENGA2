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

!   LOCAL DATA
!   ==========
    integer(kind=4), parameter :: nxglbl=3000, nyglbl=1, nzglbl=1
    integer(kind=4), parameter :: nspcmx=9
    real(kind=8) :: drun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: urun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: vrun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: wrun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: trun(nxglbl,nyglbl,nzglbl)
    real(kind=8) :: yrun(nxglbl,nyglbl,nzglbl,nspec)
    integer(kind=4) :: ic

!   BEGIN
!   =====

!   =========================================================================

    num = 111
    OPEN(UNIT=num,FILE='input/h2_init.dat',FORM="FORMATTED",STATUS= "OLD")
    DO i=1,nxglbl
        READ(num,*)ix,urun(i,1,1),trun(i,1,1),drun(i,1,1),  &
            yrun(i,1,1,1),yrun(i,1,1,2),yrun(i,1,1,3),yrun(i,1,1,4),  &
            yrun(i,1,1,5),yrun(i,1,1,6),yrun(i,1,1,7),yrun(i,1,1,8), yrun(i,1,1,9)
    END DO
    CLOSE(num)

    DO ic = istal,istol
        vrun(ic,1,1) = 0.0D0
        wrun(ic,1,1) = 0.0D0
    END DO

!   =========================================================================

END SUBROUTINE flamin

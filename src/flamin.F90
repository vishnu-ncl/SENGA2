SUBROUTINE flamin
 
! Code converted using TO_F90 by Alan Miller
! Date: 2025-01-24  Time: 09:58:37

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
!     SPECIAL FOR 25_STEP MECHAMISM

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!     -------------------------------------------------------------------------


!     PARAMETERS
!     ==========
!     ESTIMATED FLAME LOCATION AND THICKNESS
REAL(kind=8) :: clocat,cthick
PARAMETER(clocat = 0.0025_8, cthick = 0.0005_8)

!     FUNCTION
!     ========
REAL(kind=8) :: erfunc
EXTERNAL erfunc

!     LOCAL DATA
!     ==========
REAL(kind=8) :: crin(1:nxsize)
REAL(kind=8) :: yrinr(nspcmx),yrinp(nspcmx)
REAL(kind=8) :: trinr,trinp,fornow,rd
REAL(kind=8) :: deltagx,xcoord,argmntx
REAL(kind=8) :: deltagy,ycoord,argmnty
REAL(kind=8) :: deltagz,zcoord,argmntz
REAL(kind=8) :: u0,flxmas,rglocl
REAL(kind=8) :: angfrx,angfry,angfrz
REAL(kind=8) :: xrgmnt,yrgmnt,zrgmnt
REAL(kind=8) :: psi(nxsize,nysize,nzsize)

REAL(kind=8) :: cpsi,cfo
PARAMETER(cpsi=3.0_8,cfo=0.0556_8)
REAL(kind=8) :: rpsi

REAL(kind=8) :: xx
REAL(kind=8) :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11

INTEGER :: nl
PARAMETER(nl=513)
REAL :: xr(nl),tr(nl)
DOUBLE PRECISION :: xl(nl),theta(nl)

INTEGER :: icproc
INTEGER :: igofstx
INTEGER :: igofsty
INTEGER :: igofstz
INTEGER :: ix,iy,iz
INTEGER :: ic,jc,kc
INTEGER :: ispec,il

!     BEGIN
!     =====

!     =========================================================================

!     SPECIFY INITIAL THERMOCHEMICAL FIELD HERE
!     =========================================


!     SET PRODUCT TEMPERATURE
!     -----------------------
!     REACTANT TEMPERATURE SET IN CONTROL FILE
trinr = trin
u0=4.0_8
rpsi=1.0_8/8.0_8
!     GLOBAL INDEXING
!     ---------------
deltagx = xgdlen/(REAL(nxglbl-1,kind=8))
deltagy = ygdlen/(REAL(nyglbl-1,kind=8))
deltagz = zgdlen/(REAL(nzglbl-1,kind=8))

igofstx = 0
DO icproc = 0, ixproc-1
  igofstx = igofstx + npmapx(icproc)
END DO

igofsty = 0
DO icproc = 0, iyproc-1
  igofsty = igofsty + npmapy(icproc)
END DO

igofstz = 0
DO icproc = 0, izproc-1
  igofstz = igofstz + npmapz(icproc)
END DO

!     SET THE VELOCITY PROFILE FOR TGV
!     --------------------------------
angfrx = 8.0_8*ATAN(1.0_8)/xgdlen
angfry = 8.0_8*ATAN(1.0_8)/ygdlen
angfrz = 8.0_8*ATAN(1.0_8)/zgdlen
rpsi=rpsi*xgdlen

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      ix = igofstx + ic
      iy = igofsty + jc
      iz = igofstz + kc
      
      xcoord = REAL(ix-1,kind=8)*deltagx
      ycoord = REAL(iy-1,kind=8)*deltagy
      zcoord = REAL(iz-1,kind=8)*deltagz
      
      xrgmnt = angfrx*xcoord
      yrgmnt = angfry*ycoord
      zrgmnt = angfrz*zcoord
      
!       SET TAYLOR-GREEN VORTEX VELOCITY FIELD
      urun(ic,jc,kc)=u0*SIN(xrgmnt)*COS(yrgmnt)*COS(zrgmnt)
      vrun(ic,jc,kc)=-u0*COS(xrgmnt)*SIN(yrgmnt)*COS(zrgmnt)
      wrun(ic,jc,kc)=0.0_8
!       SET PRESSURE PROFILE ASSUMING CONSTANT DENSITY
      xrgmnt = 2.0_8*xrgmnt
      yrgmnt = 2.0_8*yrgmnt
      zrgmnt = 2.0_8*zrgmnt
      
      rd=ABS(xcoord-0.5_8*xgdlen)

      fornow=cpsi*(rd-rpsi)/rpsi
      fornow=(EXP(fornow)-EXP(-fornow))/ (EXP(fornow)+EXP(-fornow))
      psi(ic,jc,kc)=0.5_8*(1.0_8+fornow)
      
      
    END DO
  END DO
END DO

!     SET SPECIES MASS FRACTION PROFILES
!     ----------------------------------
DO ispec = 1, nspec
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        yrun(ic,jc,kc,ispec) = 0.0_8
        
      END DO
    END DO
  END DO
END DO

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      yrun(ic,jc,kc,1) = cfo*(1.0_8-psi(ic,jc,kc))
      yrun(ic,jc,kc,2) = 0.233_8*psi(ic,jc,kc)
      yrun(ic,jc,kc,nspec)=1.0_8-yrun(ic,jc,kc,1) -yrun(ic,jc,kc,2)
      
    END DO
  END DO
END DO

OPEN(UNIT=16,FILE='t_initial.dat',STATUS='unknown')
DO il=1,nl
  READ(16,*)xr(il),tr(il)
END DO
CLOSE(16)

DO il=1,nl
  xl(il)=REAL(xr(il))
  theta(il)=REAL(tr(il))
END DO

!     SET TEMPERATURE PROFILE
!     -----------------------
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
!      trun(ic,jc,kc) = trin
!            NEW ADDITION FOR REACTING CASE

      ix = igofstx + ic

      xcoord = REAL(ix-1)*deltagx

    DO il=1,nl-1

        IF((xl(il) <= xcoord).AND.(xl(il+1) > xcoord))THEN

          fornow=(xcoord-xl(il))/(xl(il+1)-xl(il))

          trun(ic,jc,kc)=theta(il)+fornow*(theta(il+1)-theta(il))

        END IF

      END DO

!    trun(ic,jc,kc) = max(trin,trun(ic,jc,kc))

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

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      DO ispec = 1,nspec
        
        store1(ic,jc,kc) = store1(ic,jc,kc) + rgspec(ispec)*yrun(ic,jc,kc,ispec)
        
      END DO
    END DO
  END DO
END DO

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      drun(ic,jc,kc) = prin/(store1(ic,jc,kc)*trin)
 
    END DO
  END DO
END DO

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol

    ix = igofstx + ic
    iy = igofsty + jc
    iz = igofstz + kc

    xcoord = REAL(ix-1)*deltagx
    ycoord = REAL(iy-1)*deltagy
    zcoord = REAL(iz-1)*deltagz

    xrgmnt = angfrx*xcoord
    yrgmnt = angfry*ycoord
    zrgmnt = angfrz*zcoord

!       SET PRESSURE PROFILE ASSUMING CONSTANT DENSITY
    xrgmnt = 2.0_8*xrgmnt
    yrgmnt = 2.0_8*yrgmnt
    zrgmnt = 2.0_8*zrgmnt

    prun(ic,jc,kc)=prin+((drun(ic,jc,kc)*u0*u0)/16.0_8)*  &
        (COS(xrgmnt)+COS(yrgmnt))* (COS(zrgmnt)+2.0_8)


    END DO
  END DO
END DO

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol

      drun(ic,jc,kc) = prun(ic,jc,kc)/(store1(ic,jc,kc)*trun(ic,jc,kc))

    END DO
  END DO
END DO


!CC     SET VELOCITY PROFILE ASSUMING CONSTANT MASS FLUX
!CC     --------------------
!CC     INITIAL (INLET) VEOCITY SET IN CONTROL FILE
!C      FLXMAS = DRIN*URIN
!C      DO KC = KSTAL,KSTOL
!C        DO JC = JSTAL,JSTOL
!C          DO IC = ISTAL,ISTOL
!C
!C            URUN(IC,JC,KC) = FLXMAS/DRUN(IC,JC,KC)
!C
!C          ENDDO
!C        ENDDO
!C      ENDDO
!C

!     =========================================================================

RETURN

END SUBROUTINE flamin

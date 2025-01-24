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
REAL(kind=8) :: trinr,trinp,fornow
REAL(kind=8) :: deltagx,xcoord,argmntx,rd
REAL(kind=8) :: deltagy,ycoord,argmnty
REAL(kind=8) :: deltagz,zcoord,argmntz
REAL(kind=8) :: u0,flxmas,rglocl
REAL(kind=8) :: angfrx,angfry,angfrz
REAL(kind=8) :: xrgmnt,yrgmnt,zrgmnt
REAL(kind=8) :: psi(nxsize,nysize,nzsize)

REAL(kind=8) :: cpsi,cfo
PARAMETER(cpsi=3.0_8,cfo=0.0556_8)
REAL(kind=8) :: rpsi


REAL(kind=8) :: xx,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11
INTEGER :: icproc
INTEGER :: igofstx
INTEGER :: igofsty
INTEGER :: igofstz
INTEGER :: ix,iy,iz
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
u0=34.789806_8
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
      
      prun(ic,jc,kc)=prin+((drin*u0*u0)/16.0_8)* (COS(xrgmnt)+COS(yrgmnt))*  &
          (COS(zrgmnt)+2.0_8)
      
      
      rd=ABS(xcoord-0.5_8*xgdlen)
      rpsi=rpsi*xgdlen
      fornow=cpsi*(rd-rpsi)/rpsi
      fornow=(EXP(fornow)-EXP(-fornow))/ (EXP(fornow)+EXP(-fornow))
      psi(ic,jc,kc)=0.5_8*(1.0_8+fornow)
      
      
      xx=ABS(xcoord-0.5*xgdlen)/xgdlen
      
      p1 = 307.0597_8
      p2 = -1200.3_8
      p3 = 246980.0_8
      p4 = -14803000.0_8
      p5 = 458540000.0_8
      p6 = -7629300000.0_8
      p7 = 70736000000.0_8
      p8 = 336430000000.0_8
      p9 = 637340000000.0_8
      
      IF((xx >= 0.0_8).AND.(xx <= 0.129_8)) THEN
        
        trun(ic,jc,kc) = p1 + p2*xx &
                       + p3*xx**2.0_8 &
                       + p4*xx**3.0_8 &
                       + p5*xx**4.0_8 &
                       + p6*xx**5.0_8 &
                       + p7*xx**6.0_8 &
                       + p8*xx**7.0_8 &
                       + p9*xx**8.0_8
        
      END IF
      
      p1 = -15833000.0_8
      p2 = 733200000.0_8
      p3 = -15154000000.0_8
      p4 = 184100000000.0_8
      p5 = -1456500000000.0_8
      p6 = 7842400000000.0_8
      p7 = -29107000000000.0_8
      p8 = 73543000000000.0_8
      p9 = -121070000000000.0_8
      p10 = 117300000000000.0_8
      p11 = -50789000000000.0_8      

      IF((xx >= 0.129_8).AND.(xx <= 0.3057_8)) THEN
        
        trun(ic,jc,kc) = p1 + p2*xx + p3*xx**2.0_8 &
                       + p4*xx**3.0_8 &
                       + p5*xx**4.0_8 &
                       + p6*xx**5.0_8 &
                       + p7*xx**6.0_8 &
                       + p8*xx**7.0_8 &
                       + p9*xx**8.0_8 &
                       + p10*xx**9.0_8 &
                       + p11*xx**10.0_8
        
      END IF
      
      p1 = -50428000.0
      p2 = 1273000000.0
      p3 = -14416000000.0
      p4 = 96442000000.0
      p5 = -422100000000.0
      p6 = 1262900000000.0
      p7 = -2616200000000.0
      p8 = 3705300000000.0
      p9 = -3433700000000.0
      p10 = 1880200000000.0
      p11 = -461970000000.0
      
      IF(xx > 0.3057_8) THEN
        
        fornow = p1 + p2*xx + p3*xx**2.0_8 &
               + p4*xx**3.0_8 &
               + p5*xx**4.0_8 &
               + p6*xx**5.0_8 &
               + p7*xx**6.0_8 &
               + p8*xx**7.0_8 &
               + p9*xx**8.0_8 &
               + p10*xx**9.0_8 &
               + p11*xx**10.0_8
        
        trun(ic,jc,kc)=MAX(300.0_8,fornow)
        
      END IF
      
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


!     SET TEMPERATURE PROFILE
!     -----------------------
!DO kc = kstal,kstol
!  DO jc = jstal,jstol
!    DO ic = istal,istol
      
!      trun(ic,jc,kc) = trin
      
!    END DO
!  END DO
!END DO

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
        
        store1(ic,jc,kc) = store1(ic,jc,kc) + rgspec(ispec)*yrun(ic,jc,kc,ispec)
        
      END DO
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

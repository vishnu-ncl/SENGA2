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
real(kind=8):: clocat,cthick
PARAMETER(clocat = 0.0025_8, cthick = 0.0005_8)

!C     PINCH OF HYDROGEN ATOM
!      real(kind=8)HPINCH,HLOCAT,HTHICK
!      PARAMETER(HPINCH = 1.0D-10, HLOCAT = 2.5D-3, HTHICK = 1.0D-4)
!C     PINCH OF HYDROGEN MOLECULE
!      real(kind=8)H2PNCH,H2LOCT,H2THCK
!      PARAMETER(H2PNCH = 1.0D-6, H2LOCT = 2.5D-3, H2THCK = 2.5D-4)


!     FUNCTION
!     ========
real(kind=8):: erfunc
EXTERNAL erfunc


!     LOCAL DATA
!     ==========
real(kind=8):: deltagx,deltagy,deltagz,rglocl,u0,angfrx,angfry,angfrz
real(kind=8):: crin(1:nxsize)
real(kind=8):: yrinr(nspcmx),yrinp(nspcmx)
real(kind=8):: trinr,trinp
real(kind=8):: deltag,xcoord,argmnt,ycoord,zcoord
real(kind=8):: flxmas
real(kind=8):: xrgmnt,yrgmnt,zrgmnt
INTEGER :: icproc
INTEGER :: igofst,igofsty,igofstx,igofstz
INTEGER :: ix,iy,iz
INTEGER :: ic,jc,kc
INTEGER :: ispec,jspec


!     BEGIN
!     =====

!     =========================================================================

!     SPECIFY INITIAL THERMOCHEMICAL FIELD HERE
!     =========================================


!     SET PRODUCT TEMPERATURE
!     -----------------------
!     REACTANT TEMPERATURE SET IN CONTROL FILE
trinr = trin
u0 = 34.789806_8

!Mixture gas constant
rglocl = zero
do jspec = 1,nspec
  rglocl = rglocl + rgspec(jspec)*yrin(jspec)
end do
!Times (constant) density
rglocl  = drin*rglocl

!Global indexing
!---------------
deltagx = xgdlen/(real(nxglbl-1,kind=8))
deltagy = ygdlen/(real(nyglbl-1,kind=8))
deltagz = zgdlen/(real(nzglbl-1,kind=8))

igofstx = 0
do icproc = 0, ixproc-1
  igofstx = igofstx + npmapx(icproc)
end do

igofsty = 0
do icproc = 0, iyproc-1
  igofsty = igofsty + npmapy(icproc)
end do

igofstz = 0
do icproc = 0, izproc-1
  igofstz = igofstz + npmapz(icproc)
end do

!SET THE VELOCITY PROFILE FOR TGV
!--------------------------------
angfrx = 8.0_8*atan(1.0_8)/xgdlen
angfry = 8.0_8*atan(1.0_8)/ygdlen
angfrz = 8.0_8*atan(1.0_8)/zgdlen

do kc = kstal,kstol
 do jc = jstal,jstol
  do ic = istal,istol

   ix = igofstx + ic
   iy = igofsty + jc
   iz = igofstz + kc

   xcoord = real(ix-1,kind=8)*deltagx
   ycoord = real(iy-1,kind=8)*deltagy
   zcoord = real(iz-1,kind=8)*deltagz

   xrgmnt = angfrx*xcoord
   yrgmnt = angfry*ycoord
   zrgmnt = angfrz*zcoord

!       set taylor-green vortex velocity field
   urun(ic,jc,kc)=u0*dsin(xrgmnt)*dcos(yrgmnt)*dcos(zrgmnt)
   vrun(ic,jc,kc)=-u0*dcos(xrgmnt)*dsin(yrgmnt)*dcos(zrgmnt)
   wrun(ic,jc,kc)=0.0_8
!       set pressure profile assuming constant density
   xrgmnt = 2.0_8*xrgmnt
   yrgmnt = 2.0_8*yrgmnt
   zrgmnt = 2.0_8*zrgmnt

   prun(ic,jc,kc)=prin+((drin*u0*u0)/16.0_8)*(dcos(xrgmnt)+dcos(yrgmnt))*(dcos(zrgmnt)+2.0_8)
  end do
 end do
end do

!set temperature profile assuming constant density

do kc = kstal,kstol
  do jc = jstal,jstol
    do ic = istal,istol

      trun(ic,jc,kc) = prun(ic,jc,kc)/rglocl

    end do
  end do
end do


!set constant density

do kc = kstal,kstol
  do jc = jstal,jstol
    do ic = istal,istol

      drun(ic,jc,kc) = drin

    end do
  end do
end do


RETURN

9000  FORMAT(a)

END SUBROUTINE flamin

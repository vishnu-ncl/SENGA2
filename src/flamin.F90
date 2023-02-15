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
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     PARAMETERS
!     ==========
!     ESTIMATED FLAME LOCATION AND THICKNESS
real(kind=dp) :: clocat,cthick
PARAMETER(clocat = 0.0025_dp, cthick = 0.0005_dp)

!C     PINCH OF HYDROGEN ATOM
!      real(kind=dp) HPINCH,HLOCAT,HTHICK
!      PARAMETER(HPINCH = 1.0D-10, HLOCAT = 2.5D-3, HTHICK = 1.0D-4)
!C     PINCH OF HYDROGEN MOLECULE
!      real(kind=dp) H2PNCH,H2LOCT,H2THCK
!      PARAMETER(H2PNCH = 1.0D-6, H2LOCT = 2.5D-3, H2THCK = 2.5D-4)


!     FUNCTION
!     ========
real(kind=dp) :: erfunc
EXTERNAL erfunc


!     LOCAL DATA
!     ==========
real(kind=dp) :: crin(1:nxsize)
real(kind=dp) :: yrinr(nspcmx),yrinp(nspcmx)
real(kind=dp) :: trinr,trinp
real(kind=dp) :: deltag,xcoord,argmnt
real(kind=dp) :: flxmas
real(kind=dp) :: xrgmnt,yrgmnt,zrgmnt
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
u0 = 34.789806_dp

!Mixture gas constant
rglocl = zero
do jspec = 1,nspec
  rglocl = rglocl + rgspec(jspec)*yrin(jspec)
end do
!Times (constant) density
rglocl  = drin*rglocl

!Global indexing
!---------------
deltagx = xgdlen/(real(nxglbl-1,dp))
deltagy = ygdlen/(real(nyglbl-1,dp))
deltagz = zgdlen/(real(nzglbl-1,dp))

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
angfrx = 8.0_dp*atan(1.0_dp)/xgdlen
angfry = 8.0_dp*atan(1.0_dp)/ygdlen
angfrz = 8.0_dp*atan(1.0_dp)/zgdlen

do kc = kstal,kstol
 do jc = jstal,jstol
  do ic = istal,istol

   ix = igofstx + ic
   iy = igofsty + jc
   iz = igofstz + kc

   xcoord = real(ix-1,dp)*deltagx
   ycoord = real(iy-1,dp)*deltagy
   zcoord = real(iz-1,dp)*deltagz

   xrgmnt = angfrx*xcoord
   yrgmnt = angfry*ycoord
   zrgmnt = angfrz*zcoord

!       set taylor-green vortex velocity field
   urun(ic,jc,kc)=u0*dsin(xrgmnt)*dcos(yrgmnt)*dcos(zrgmnt)
   vrun(ic,jc,kc)=-u0*dcos(xrgmnt)*dsin(yrgmnt)*dcos(zrgmnt)
   wrun(ic,jc,kc)=0.0_dp
!       set pressure profile assuming constant density
   xrgmnt = 2.0_dp*xrgmnt
   yrgmnt = 2.0_dp*yrgmnt
   zrgmnt = 2.0_dp*zrgmnt

   prun(ic,jc,kc)=prin+((drin*u0*u0)/16.0_dp)*(dcos(xrgmnt)+dcos(yrgmnt))*(dcos(zrgmnt)+2.0_dp)
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

SUBROUTINE pardom
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-04  Time: 20:08:33

!     *************************************************************************

!     PARDOM
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     01-OCT-1996:  CREATED
!     09-MAR-2003:  RSC UPDATED FOR SENGA2
!     27-SEP-2009:  RSC UPDATE INDEXING FOR FILTERS

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     CARRIES OUT PARALLEL DOMAIN DECOMPOSITION

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ===========
INTEGER :: nextra
INTEGER :: icproc


!     BEGIN
!     =====

!     =========================================================================

!     INITIALISE PARALLEL MESSAGE PASSING
!     ===================================
!     INITIALISE
CALL p_init

!     GET THE NUMBER OF PROCESSORS NPROC
!     AND THE LOCAL PROCESSOR ID IPROC: IPROC=0,NPROC-1
CALL p_size(nproc,iproc)

!     =========================================================================

!     GET XYZ COORDINATES OF LOCAL PROCESSOR
!     ======================================
!     IPROC = IXPROC+IYPROC*NXPROC+IZPROC*NXPROC*NYPROC
ixproc = MOD(iproc,nxproc)
iyproc = iproc/nxproc
izproc = iyproc/nyproc
iyproc = MOD(iyproc,nyproc)

!     =========================================================================

!     ASSIGN GRID POINTS TO PROCESSORS
!     ================================
!     STATIC LOAD-BALANCING:
!     LOAD IS AUTOMATICALLY BALANCED IF NXGLBL IS EXACTLY DIVISIBLE BY NXPROC
!     OTHERWISE ANY REMAINDER IS SPLIT EQUALLY BETWEEN THE HIGHER-RANKED
!     PROCESSORS IN EACH DIRECTION
nxnode = nxglbl/nxproc
nextra = MOD(nxglbl,nxproc)
nextra = nxproc-nextra
nxprm1 = nxproc-1
DO icproc = 0,nxprm1
  npmapx(icproc) = nxnode
  IF(icproc >= nextra)npmapx(icproc) = nxnode+1
END DO
nxnode = npmapx(ixproc)
nxnbig = nxnode+2*nhalox

nynode = nyglbl/nyproc
nextra = MOD(nyglbl,nyproc)
nextra = nyproc-nextra
nyprm1 = nyproc-1
DO icproc = 0,nyprm1
  npmapy(icproc) = nynode
  IF(icproc >= nextra)npmapy(icproc) = nynode+1
END DO
nynode = npmapy(iyproc)
nynbig = nynode+2*nhaloy

nznode = nzglbl/nzproc
nextra = MOD(nzglbl,nzproc)
nextra = nzproc-nextra
nzprm1 = nzproc-1
DO icproc = 0,nzprm1
  npmapz(icproc) = nznode
  IF(icproc >= nextra)npmapz(icproc) = nznode+1
END DO
nznode = npmapz(izproc)
nznbig = nznode+2*nhaloz

!     =========================================================================

!     MAP PROCESSORS ALONG CURRENT ROW
!     ================================
DO icproc = 0,nxprm1
  nprocx(icproc) = icproc+nxproc*(iyproc+nyproc*izproc)
END DO
DO icproc = 0,nyprm1
  nprocy(icproc) = ixproc+nxproc*(icproc+nyproc*izproc)
END DO
DO icproc = 0,nzprm1
  nprocz(icproc) = ixproc+nxproc*(iyproc+nyproc*icproc)
END DO

!     =========================================================================

!     IDENTIFY NEAREST NEIGHBOURS FOR MESSAGE-PASSING
!     ===============================================
ixprom = ixproc-1
IF(ixprom < 0)ixprom = nxproc-1
ixprop = ixproc+1
IF(ixprop >= nxproc)ixprop = 0
ixprom = ixprom+nxproc*(iyproc+izproc*nyproc)
ixprop = ixprop+nxproc*(iyproc+izproc*nyproc)

iyprom = iyproc-1
IF(iyprom < 0)iyprom = nyproc-1
iyprop = iyproc+1
IF(iyprop >= nyproc)iyprop = 0
iyprom = ixproc+nxproc*(iyprom+izproc*nyproc)
iyprop = ixproc+nxproc*(iyprop+izproc*nyproc)

izprom = izproc-1
IF(izprom < 0)izprom = nzproc-1
izprop = izproc+1
IF(izprop >= nzproc)izprop = 0
izprom = ixproc+nxproc*(iyproc+izprom*nyproc)
izprop = ixproc+nxproc*(iyproc+izprop*nyproc)

!     =========================================================================

!     SET FLAGS AND TAGS FOR MESSAGE-PASSING
!     ======================================
!     LOCAL PROCESSOR EVEN/ODD FLAGS
proddx = .true.
IF(MOD(ixproc,2) == 0)proddx = .false.
proddy = .true.
IF(MOD(iyproc,2) == 0)proddy = .false.
proddz = .true.
IF(MOD(izproc,2) == 0)proddz = .false.

!     SET MESSAGE TAGS
!     EACH TAG UNIQUELY IDENTIFIES THE SENDING AND RECEIVING PROCESS
!     FOR EACH PROCESSOR FOR EACH DIRECTION THERE ARE FOUR TAGS
!     CORRESPONDING TO SEND AND RECEIVE LEFT AND RIGHT
itgxsl = iproc*nproc + ixprom
itgxrl = ixprom*nproc + iproc
itgxsr = iproc*nproc + ixprop
itgxrr = ixprop*nproc + iproc
itgysl = iproc*nproc + iyprom
itgyrl = iyprom*nproc + iproc
itgysr = iproc*nproc + iyprop
itgyrr = iyprop*nproc + iproc
itgzsl = iproc*nproc + izprom
itgzrl = izprom*nproc + iproc
itgzsr = iproc*nproc + izprop
itgzrr = izprop*nproc + iproc

!     =========================================================================

!     SET LIMITS ON SPATIAL COUNTERS
!     ==============================
!     ON THE LOCAL PROCESSOR

!     STANDARD SIZE ARRAYS
istal = 1
istol = nxnode
jstal = 1
jstol = nynode
kstal = 1
kstol = nznode

istap1 = istal+1
istom1 = istol-1
jstap1 = jstal+1
jstom1 = jstol-1
kstap1 = kstal+1
kstom1 = kstol-1

istap2 = istal+2
istom2 = istol-2
jstap2 = jstal+2
jstom2 = jstol-2
kstap2 = kstal+2
kstom2 = kstol-2

istap3 = istal+3
istom3 = istol-3
jstap3 = jstal+3
jstom3 = jstol-3
kstap3 = kstal+3
kstom3 = kstol-3

istap4 = istal+4
istom4 = istol-4
jstap4 = jstal+4
jstom4 = jstol-4
kstap4 = kstal+4
kstom4 = kstol-4

istap5 = istal+5
istom5 = istol-5
jstap5 = jstal+5
jstom5 = jstol-5
kstap5 = kstal+5
kstom5 = kstol-5

istap6 = istal+6
istom6 = istol-6
jstap6 = jstal+6
jstom6 = jstol-6
kstap6 = kstal+6
kstom6 = kstol-6

istap7 = istal+7
istom7 = istol-7
jstap7 = jstal+7
jstom7 = jstol-7
kstap7 = kstal+7
kstom7 = kstol-7

istap8 = istal+8
istom8 = istol-8
jstap8 = jstal+8
jstom8 = jstol-8
kstap8 = kstal+8
kstom8 = kstol-8

istap9 = istal+9
istom9 = istol-9
jstap9 = jstal+9
jstom9 = jstol-9
kstap9 = kstal+9
kstom9 = kstol-9

istapa = istal+10
istoma = istol-10
jstapa = jstal+10
jstoma = jstol-10
kstapa = kstal+10
kstoma = kstol-10

istapb = istal+11
istomb = istol-11
jstapb = jstal+11
jstomb = jstol-11
kstapb = kstal+11
kstomb = kstol-11

!     BIGGER SIZE ARRAYS
istab = 1-nhalox
istob = nxnode+nhalox
jstab = 1-nhaloy
jstob = nynode+nhaloy
kstab = 1-nhaloz
kstob = nznode+nhaloz

!     HALO DATA
istali = 1
istoli = nhalox
jstali = 1
jstoli = nhaloy
kstali = 1
kstoli = nhaloz
istari = nxnode+1-nhalox
istori = nxnode
jstari = nynode+1-nhaloy
jstori = nynode
kstari = nznode+1-nhaloz
kstori = nznode
istalo = 1-nhalox
istolo = 0
jstalo = 1-nhaloy
jstolo = 0
kstalo = 1-nhaloz
kstolo = 0
istaro = nxnode+1
istoro = nxnode+nhalox
jstaro = nynode+1
jstoro = nynode+nhaloy
kstaro = nznode+1
kstoro = nznode+nhaloz

!     =========================================================================


RETURN
END SUBROUTINE pardom

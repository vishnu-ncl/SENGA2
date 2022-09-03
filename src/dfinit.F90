SUBROUTINE dfinit
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-02  Time: 15:43:38

!     *************************************************************************

!     DFINIT
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     28-MAR-2003:  RSC MODIFIED FOR SENGA2
!     10-NOV-2013:  RSC MODIFIED FOR MESH STRETCHING

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     INITIALISES SPATIAL DIFFERENTIATORS
!     10TH ORDER EXPLICIT DIFFERENCING
!     8TH,6TH,4TH,4TH ORDER EXPLICIT BOUNDARY SCHEMES

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
       use com_senga
!     -------------------------------------------------------------------------


!     BEGIN
!     =====

!     =========================================================================

!     SPATIAL STEP SIZES
!     ==================
deltax = xgdlen/DBLE(nxglbl-1)
ovdelx = one/deltax
ovdlx2 = ovdelx*ovdelx
deltay = ygdlen/DBLE(nyglbl-1)
ovdely = one/deltay
ovdly2 = ovdely*ovdely
deltaz = zgdlen/DBLE(nzglbl-1)
ovdelz = one/deltaz
ovdlz2 = ovdelz*ovdelz

!     =========================================================================

!     FIRST DERIVATIVES
!     =================

!     INTERIOR SCHEME
!     ---------------

!     TENTH ORDER EXPLICIT CENTRED DIFFERENCES
acoeff = 5.0_8/3.0_8
bcoeff = -20.0_8/21.0_8
ccoeff = 5.0_8/14.0_8
dcoeff = -5.0_8/63.0_8
ecoeff = 1.0_8/126.0_8

acoffx = acoeff/2.0_8
bcoffx = bcoeff/4.0_8
ccoffx = ccoeff/6.0_8
dcoffx = dcoeff/8.0_8
ecoffx = ecoeff/10.0_8

acoffy = acoeff/2.0_8
bcoffy = bcoeff/4.0_8
ccoffy = ccoeff/6.0_8
dcoffy = dcoeff/8.0_8
ecoffy = ecoeff/10.0_8

acoffz = acoeff/2.0_8
bcoffz = bcoeff/4.0_8
ccoffz = ccoeff/6.0_8
dcoffz = dcoeff/8.0_8
ecoffz = ecoeff/10.0_8

!     BOUNDARY TREATMENT
!     ------------------

!     FIRST POINT SCHEME (4TH ORDER ONE SIDED)
acoef1 = 4.0_8
bcoef1 = -3.0_8
ccoef1 = 4.0_8/3.0_8
dcoef1 = -1.0_8/4.0_8

acof1x = acoef1
bcof1x = bcoef1
ccof1x = ccoef1
dcof1x = dcoef1

acof1y = acoef1
bcof1y = bcoef1
ccof1y = ccoef1
dcof1y = dcoef1

acof1z = acoef1
bcof1z = bcoef1
ccof1z = ccoef1
dcof1z = dcoef1

!     SECOND POINT SCHEME (4TH ORDER MIXED)
acoef2 = -1.0_8/4.0_8
bcoef2 = 3.0_8/2.0_8
ccoef2 = -1.0_8/2.0_8
dcoef2 = 1.0_8/12.0_8

acof2x = acoef2
bcof2x = bcoef2
ccof2x = ccoef2
dcof2x = dcoef2

acof2y = acoef2
bcof2y = bcoef2
ccof2y = ccoef2
dcof2y = dcoef2

acof2z = acoef2
bcof2z = bcoef2
ccof2z = ccoef2
dcof2z = dcoef2

!     3RD POINT SCHEME (4TH ORDER EXPLICIT CENTRED)
acoef3 = 4.0_8/3.0_8
bcoef3 = -1.0_8/3.0_8

acof3x = acoef3/2.0_8
bcof3x = bcoef3/4.0_8

acof3y = acoef3/2.0_8
bcof3y = bcoef3/4.0_8

acof3z = acoef3/2.0_8
bcof3z = bcoef3/4.0_8

!     4TH POINT SCHEME (6TH ORDER EXPLICIT CENTRED)
acoef4 = 3.0_8/2.0_8
bcoef4 = -3.0_8/5.0_8
ccoef4 = 1.0_8/10.0_8

acof4x = acoef4/2.0_8
bcof4x = bcoef4/4.0_8
ccof4x = ccoef4/6.0_8

acof4y = acoef4/2.0_8
bcof4y = bcoef4/4.0_8
ccof4y = ccoef4/6.0_8

acof4z = acoef4/2.0_8
bcof4z = bcoef4/4.0_8
ccof4z = ccoef4/6.0_8

!     5TH POINT SCHEME (8TH ORDER EXPLICIT CENTRED)
acoef5 = 8.0_8/5.0_8
bcoef5 = -4.0_8/5.0_8
ccoef5 = 8.0_8/35.0_8
dcoef5 = -1.0_8/35.0_8

acof5x = acoef5/2.0_8
bcof5x = bcoef5/4.0_8
ccof5x = ccoef5/6.0_8
dcof5x = dcoef5/8.0_8

acof5y = acoef5/2.0_8
bcof5y = bcoef5/4.0_8
ccof5y = ccoef5/6.0_8
dcof5y = dcoef5/8.0_8

acof5z = acoef5/2.0_8
bcof5z = bcoef5/4.0_8
ccof5z = ccoef5/6.0_8
dcof5z = dcoef5/8.0_8

!     =========================================================================

!     SECOND DERIVATIVES
!     ==================

!     INTERIOR SCHEME
!     ---------------

!     TENTH ORDER EXPLICIT CENTRED DIFFERENCES
acoefs = acoeff
bcoefs = bcoeff
ccoefs = ccoeff
dcoefs = dcoeff
ecoefs = ecoeff

acofsx = acoefs
bcofsx = bcoefs/4.0_8
ccofsx = ccoefs/9.0_8
dcofsx = dcoefs/16.0_8
ecofsx = ecoefs/25.0_8

acofsy = acoefs
bcofsy = bcoefs/4.0_8
ccofsy = ccoefs/9.0_8
dcofsy = dcoefs/16.0_8
ecofsy = ecoefs/25.0_8

acofsz = acoefs
bcofsz = bcoefs/4.0_8
ccofsz = ccoefs/9.0_8
dcofsz = dcoefs/16.0_8
ecofsz = ecoefs/25.0_8


!     BOUNDARY TREATMENT
!     ------------------

!     FIRST POINT SCHEME (4TH ORDER ONE SIDED)
acofs1 = -77.0_8/6.0_8
bcofs1 = 107.0_8/6.0_8
ccofs1 = -13.0_8
dcofs1 = 61.0_8/12.0_8
ecofs1 = -5.0_8/6.0_8

acfs1x = acofs1
bcfs1x = bcofs1
ccfs1x = ccofs1
dcfs1x = dcofs1
ecfs1x = ecofs1

acfs1y = acofs1
bcfs1y = bcofs1
ccfs1y = ccofs1
dcfs1y = dcofs1
ecfs1y = ecofs1

acfs1z = acofs1
bcfs1z = bcofs1
ccfs1z = ccofs1
dcfs1z = dcofs1
ecfs1z = ecofs1

!     SECOND POINT SCHEME (4TH ORDER MIXED)
acofs2 = 5.0_8/6.0_8
bcofs2 = -1.0_8/3.0_8
ccofs2 = 7.0_8/6.0_8
dcofs2 = -1.0_8/2.0_8
ecofs2 = 1.0_8/12.0_8

acfs2x = acofs2
bcfs2x = bcofs2
ccfs2x = ccofs2
dcfs2x = dcofs2
ecfs2x = ecofs2

acfs2y = acofs2
bcfs2y = bcofs2
ccfs2y = ccofs2
dcfs2y = dcofs2
ecfs2y = ecofs2

acfs2z = acofs2
bcfs2z = bcofs2
ccfs2z = ccofs2
dcfs2z = dcofs2
ecfs2z = ecofs2

!     3RD POINT SCHEME (4TH ORDER EXPLICIT CENTRED)
acofs3 = acoef3
bcofs3 = bcoef3

acfs3x = acofs3
bcfs3x = bcofs3/4.0_8

acfs3y = acofs3
bcfs3y = bcofs3/4.0_8

acfs3z = acofs3
bcfs3z = bcofs3/4.0_8

!     4TH POINT SCHEME (6TH ORDER EXPLICIT CENTRED)
acofs4 = acoef4
bcofs4 = bcoef4
ccofs4 = ccoef4

acfs4x = acofs4
bcfs4x = bcofs4/4.0_8
ccfs4x = ccofs4/9.0_8

acfs4y = acofs4
bcfs4y = bcofs4/4.0_8
ccfs4y = ccofs4/9.0_8

acfs4z = acofs4
bcfs4z = bcofs4/4.0_8
ccfs4z = ccofs4/9.0_8

!     5TH POINT SCHEME (8TH ORDER EXPLICIT CENTRED)
acofs5 = acoef5
bcofs5 = bcoef5
ccofs5 = ccoef5
dcofs5 = dcoef5

acfs5x = acofs5
bcfs5x = bcofs5/4.0_8
ccfs5x = ccofs5/9.0_8
dcfs5x = dcofs5/16.0_8

acfs5y = acofs5
bcfs5y = bcofs5/4.0_8
ccfs5y = ccofs5/9.0_8
dcfs5y = dcofs5/16.0_8

acfs5z = acofs5
bcfs5z = bcofs5/4.0_8
ccfs5z = ccofs5/9.0_8
dcfs5z = dcofs5/16.0_8

!     =========================================================================

!     SECOND CROSS-DERIVATIVES
!     ========================

!     TENTH ORDER EXPLICIT CENTRED DIFFERENCES
acoefx = acoeff
bcoefx = bcoeff
ccoefx = ccoeff
dcoefx = dcoeff
ecoefx = ecoeff

acofxy = acoefx/4.0_8
bcofxy = bcoefx/4.0_8/4.0_8
ccofxy = ccoefx/4.0_8/9.0_8
dcofxy = dcoefx/4.0_8/16.0_8
ecofxy = ecoefx/4.0_8/25.0_8

acofxz = acoefx/4.0_8
bcofxz = bcoefx/4.0_8/4.0_8
ccofxz = ccoefx/4.0_8/9.0_8
dcofxz = dcoefx/4.0_8/16.0_8
ecofxz = ecoefx/4.0_8/25.0_8

acofyz = acoefx/4.0_8
bcofyz = bcoefx/4.0_8/4.0_8
ccofyz = ccoefx/4.0_8/9.0_8
dcofyz = dcoefx/4.0_8/16.0_8
ecofyz = ecoefx/4.0_8/25.0_8


!     BOUNDARY TREATMENT
!     ------------------
!     FIRST/SECOND POINT SCHEME (4ND ORDER CENTRED IN TRANSVERSE DIRECTION)
acofx1 = acoef3/2.0_8
bcofx1 = bcoef3/4.0_8

acofy1 = acoef3/2.0_8
bcofy1 = bcoef3/4.0_8

acofz1 = acoef3/2.0_8
bcofz1 = bcoef3/4.0_8

!     FIRST POINT SCHEME (4ND ORDER ONE SIDED/CENTRED)
acf1xy = acoef1
bcf1xy = bcoef1
ccf1xy = ccoef1
dcf1xy = dcoef1

acf1xz = acoef1
bcf1xz = bcoef1
ccf1xz = ccoef1
dcf1xz = dcoef1

acf1yz = acoef1
bcf1yz = bcoef1
ccf1yz = ccoef1
dcf1yz = dcoef1

!     SECOND POINT SCHEME (4TH ORDER MIXED/CENTRED)
acf2xy = acoef2
bcf2xy = bcoef2
ccf2xy = ccoef2
dcf2xy = dcoef2

acf2xz = acoef2
bcf2xz = bcoef2
ccf2xz = ccoef2
dcf2xz = dcoef2

acf2yz = acoef2
bcf2yz = bcoef2
ccf2yz = ccoef2
dcf2yz = dcoef2

!     THIRD POINT SCHEME (4TH ORDER EXPLICIT CENTRED)
acf3xy = acoef3/4.0_8
bcf3xy = bcoef3/4.0_8/4.0_8

acf3xz = acoef3/4.0_8
bcf3xz = bcoef3/4.0_8/4.0_8

acf3yz = acoef3/4.0_8
bcf3yz = bcoef3/4.0_8/4.0_8

!     FOURTH POINT SCHEME (6TH ORDER EXPLICIT CENTRED)
acf4xy = acoef4/4.0_8
bcf4xy = bcoef4/4.0_8/4.0_8
ccf4xy = ccoef4/4.0_8/9.0_8

acf4xz = acoef4/4.0_8
bcf4xz = bcoef4/4.0_8/4.0_8
ccf4xz = ccoef4/4.0_8/9.0_8

acf4yz = acoef4/4.0_8
bcf4yz = bcoef4/4.0_8/4.0_8
ccf4yz = ccoef4/4.0_8/9.0_8

!     FIFTH POINT SCHEME (8TH ORDER EXPLICIT CENTRED)
acf5xy = acoef5/4.0_8
bcf5xy = bcoef5/4.0_8/4.0_8
ccf5xy = ccoef5/4.0_8/9.0_8
dcf5xy = dcoef5/4.0_8/16.0_8

acf5xz = acoef5/4.0_8
bcf5xz = bcoef5/4.0_8/4.0_8
ccf5xz = ccoef5/4.0_8/9.0_8
dcf5xz = dcoef5/4.0_8/16.0_8

acf5yz = acoef5/4.0_8
bcf5yz = bcoef5/4.0_8/4.0_8
ccf5yz = ccoef5/4.0_8/9.0_8
dcf5yz = dcoef5/4.0_8/16.0_8

!     CORNER POINT SCHEME (4TH ORDER ONE SIDED/ONE SIDED)
acofc1 = acoef1
bcofc1 = bcoef1*2.0_8
ccofc1 = ccoef1*3.0_8
dcofc1 = dcoef1*4.0_8

acc1xy = acofc1
bcc1xy = bcofc1/4.0_8
ccc1xy = ccofc1/9.0_8
dcc1xy = dcofc1/16.0_8

acc1xz = acofc1
bcc1xz = bcofc1/4.0_8
ccc1xz = ccofc1/9.0_8
dcc1xz = dcofc1/16.0_8

acc1yz = acofc1
bcc1yz = bcofc1/4.0_8
ccc1yz = ccofc1/9.0_8
dcc1yz = dcofc1/16.0_8

!     SECOND CORNER POINT SCHEME (4TH ORDER MIXED)
acofc2 = -acoef2
bcofc2 =  bcoef2
ccofc2 =  ccoef2*2.0_8
dcofc2 =  dcoef2*3.0_8

acc2xy = acofc2
bcc2xy = bcofc2
ccc2xy = ccofc2/4.0_8
dcc2xy = dcofc2/9.0_8

acc2xz = acofc2
bcc2xz = bcofc2
ccc2xz = ccofc2/4.0_8
dcc2xz = dcofc2/9.0_8

acc2yz = acofc2
bcc2yz = bcofc2
ccc2yz = ccofc2/4.0_8
dcc2yz = dcofc2/9.0_8

!     SECOND EDGE POINT SCHEME (4TH ORDER MIXED/CENTRED)
!     USES FIRST AND SECOND POINT COEFFS

!     =========================================================================

!     RSC 10-NOV-2013 INITIALISE MESH STRETCHING
CALL dfmstr

!     =========================================================================

!     DIFFERENCE COEFFICIENTS FOR WALL BCS
acbcxl(1) = acoef1*ovdelx
acbcxl(2) = bcoef1*ovdelx
acbcxl(3) = ccoef1*ovdelx
acbcxl(4) = dcoef1*ovdelx

acbcxr(1) = acoef1*ovdelx
acbcxr(2) = bcoef1*ovdelx
acbcxr(3) = ccoef1*ovdelx
acbcxr(4) = dcoef1*ovdelx

acbcyl(1) = acoef1*ovdely
acbcyl(2) = bcoef1*ovdely
acbcyl(3) = ccoef1*ovdely
acbcyl(4) = dcoef1*ovdely

acbcyr(1) = acoef1*ovdely
acbcyr(2) = bcoef1*ovdely
acbcyr(3) = ccoef1*ovdely
acbcyr(4) = dcoef1*ovdely

acbczl(1) = acoef1*ovdelz
acbczl(2) = bcoef1*ovdelz
acbczl(3) = ccoef1*ovdelz
acbczl(4) = dcoef1*ovdelz

acbczr(1) = acoef1*ovdelz
acbczr(2) = bcoef1*ovdelz
acbczr(3) = ccoef1*ovdelz
acbczr(4) = dcoef1*ovdelz

!C     ONE-SIDED MESH STRETCHING IN Y DIRECTION
!      ACBCYL(1) = ACOEF1*OVDELY*DGDHAT(1)
!      ACBCYL(2) = BCOEF1*OVDELY*DGDHAT(1)
!      ACBCYL(3) = CCOEF1*OVDELY*DGDHAT(1)
!      ACBCYL(4) = DCOEF1*OVDELY*DGDHAT(1)

!     =========================================================================


RETURN
END SUBROUTINE dfinit

SUBROUTINE dfinit

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga
 
!   *************************************************************************

!   DFINIT
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   28-MAR-2003:  RSC MODIFIED FOR SENGA2
!   10-NOV-2013:  RSC MODIFIED FOR MESH STRETCHING

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   INITIALISES SPATIAL DIFFERENTIATORS
!   10TH ORDER EXPLICIT DIFFERENCING
!   8TH,6TH,4TH,4TH ORDER EXPLICIT BOUNDARY SCHEMES

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   BEGIN
!   =====

!   =========================================================================

!   SPATIAL STEP SIZES
!   ==================
    deltax = xgdlen/DBLE(nxglbl-1)
    ovdelx = one/deltax
    ovdlx2 = ovdelx*ovdelx
    deltay = ygdlen/DBLE(nyglbl-1)
    ovdely = one/deltay
    ovdly2 = ovdely*ovdely
    deltaz = zgdlen/DBLE(nzglbl-1)
    ovdelz = one/deltaz
    ovdlz2 = ovdelz*ovdelz

!   =========================================================================

!   FIRST DERIVATIVES
!   =================

!   INTERIOR SCHEME
!   ---------------

!   TENTH ORDER EXPLICIT CENTRED DIFFERENCES
    acoeff = 5.0_dp/3.0_dp
    bcoeff = -20.0_dp/21.0_dp
    ccoeff = 5.0_dp/14.0_dp
    dcoeff = -5.0_dp/63.0_dp
    ecoeff = 1.0_dp/126.0_dp

    acoffx = acoeff/2.0_dp
    bcoffx = bcoeff/4.0_dp
    ccoffx = ccoeff/6.0_dp
    dcoffx = dcoeff/8.0_dp
    ecoffx = ecoeff/10.0_dp

    acoffy = acoeff/2.0_dp
    bcoffy = bcoeff/4.0_dp
    ccoffy = ccoeff/6.0_dp
    dcoffy = dcoeff/8.0_dp
    ecoffy = ecoeff/10.0_dp

    acoffz = acoeff/2.0_dp
    bcoffz = bcoeff/4.0_dp
    ccoffz = ccoeff/6.0_dp
    dcoffz = dcoeff/8.0_dp
    ecoffz = ecoeff/10.0_dp

!   BOUNDARY TREATMENT
!   ------------------

!   FIRST POINT SCHEME (4TH ORDER ONE SIDED)
    acoef1 = 4.0_dp
    bcoef1 = -3.0_dp
    ccoef1 = 4.0_dp/3.0_dp
    dcoef1 = -1.0_dp/4.0_dp

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

!   SECOND POINT SCHEME (4TH ORDER MIXED)
    acoef2 = -1.0_dp/4.0_dp
    bcoef2 = 3.0_dp/2.0_dp
    ccoef2 = -1.0_dp/2.0_dp
    dcoef2 = 1.0_dp/12.0_dp

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

!   3RD POINT SCHEME (4TH ORDER EXPLICIT CENTRED)
    acoef3 = 4.0_dp/3.0_dp
    bcoef3 = -1.0_dp/3.0_dp

    acof3x = acoef3/2.0_dp
    bcof3x = bcoef3/4.0_dp

    acof3y = acoef3/2.0_dp
    bcof3y = bcoef3/4.0_dp

    acof3z = acoef3/2.0_dp
    bcof3z = bcoef3/4.0_dp

!   4TH POINT SCHEME (6TH ORDER EXPLICIT CENTRED)
    acoef4 = 3.0_dp/2.0_dp
    bcoef4 = -3.0_dp/5.0_dp
    ccoef4 = 1.0_dp/10.0_dp

    acof4x = acoef4/2.0_dp
    bcof4x = bcoef4/4.0_dp
    ccof4x = ccoef4/6.0_dp

    acof4y = acoef4/2.0_dp
    bcof4y = bcoef4/4.0_dp
    ccof4y = ccoef4/6.0_dp

    acof4z = acoef4/2.0_dp
    bcof4z = bcoef4/4.0_dp
    ccof4z = ccoef4/6.0_dp

!   5TH POINT SCHEME (8TH ORDER EXPLICIT CENTRED)
    acoef5 = 8.0_dp/5.0_dp
    bcoef5 = -4.0_dp/5.0_dp
    ccoef5 = 8.0_dp/35.0_dp
    dcoef5 = -1.0_dp/35.0_dp

    acof5x = acoef5/2.0_dp
    bcof5x = bcoef5/4.0_dp
    ccof5x = ccoef5/6.0_dp
    dcof5x = dcoef5/8.0_dp

    acof5y = acoef5/2.0_dp
    bcof5y = bcoef5/4.0_dp
    ccof5y = ccoef5/6.0_dp
    dcof5y = dcoef5/8.0_dp

    acof5z = acoef5/2.0_dp
    bcof5z = bcoef5/4.0_dp
    ccof5z = ccoef5/6.0_dp
    dcof5z = dcoef5/8.0_dp

!   =========================================================================

!   SECOND DERIVATIVES
!   ==================

!   INTERIOR SCHEME
!   ---------------

!   TENTH ORDER EXPLICIT CENTRED DIFFERENCES
    acoefs = acoeff
    bcoefs = bcoeff
    ccoefs = ccoeff
    dcoefs = dcoeff
    ecoefs = ecoeff

    acofsx = acoefs
    bcofsx = bcoefs/4.0_dp
    ccofsx = ccoefs/9.0_dp
    dcofsx = dcoefs/16.0_dp
    ecofsx = ecoefs/25.0_dp

    acofsy = acoefs
    bcofsy = bcoefs/4.0_dp
    ccofsy = ccoefs/9.0_dp
    dcofsy = dcoefs/16.0_dp
    ecofsy = ecoefs/25.0_dp

    acofsz = acoefs
    bcofsz = bcoefs/4.0_dp
    ccofsz = ccoefs/9.0_dp
    dcofsz = dcoefs/16.0_dp
    ecofsz = ecoefs/25.0_dp

!   BOUNDARY TREATMENT
!   ------------------

!   FIRST POINT SCHEME (4TH ORDER ONE SIDED)
    acofs1 = -77.0_dp/6.0_dp
    bcofs1 = 107.0_dp/6.0_dp
    ccofs1 = -13.0_dp
    dcofs1 = 61.0_dp/12.0_dp
    ecofs1 = -5.0_dp/6.0_dp

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

!   SECOND POINT SCHEME (4TH ORDER MIXED)
    acofs2 = 5.0_dp/6.0_dp
    bcofs2 = -1.0_dp/3.0_dp
    ccofs2 = 7.0_dp/6.0_dp
    dcofs2 = -1.0_dp/2.0_dp
    ecofs2 = 1.0_dp/12.0_dp

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

!   3RD POINT SCHEME (4TH ORDER EXPLICIT CENTRED)
    acofs3 = acoef3
    bcofs3 = bcoef3

    acfs3x = acofs3
    bcfs3x = bcofs3/4.0_dp

    acfs3y = acofs3
    bcfs3y = bcofs3/4.0_dp

    acfs3z = acofs3
    bcfs3z = bcofs3/4.0_dp

!   4TH POINT SCHEME (6TH ORDER EXPLICIT CENTRED)
    acofs4 = acoef4
    bcofs4 = bcoef4
    ccofs4 = ccoef4

    acfs4x = acofs4
    bcfs4x = bcofs4/4.0_dp
    ccfs4x = ccofs4/9.0_dp

    acfs4y = acofs4
    bcfs4y = bcofs4/4.0_dp
    ccfs4y = ccofs4/9.0_dp

    acfs4z = acofs4
    bcfs4z = bcofs4/4.0_dp
    ccfs4z = ccofs4/9.0_dp

!   5TH POINT SCHEME (8TH ORDER EXPLICIT CENTRED)
    acofs5 = acoef5
    bcofs5 = bcoef5
    ccofs5 = ccoef5
    dcofs5 = dcoef5

    acfs5x = acofs5
    bcfs5x = bcofs5/4.0_dp
    ccfs5x = ccofs5/9.0_dp
    dcfs5x = dcofs5/16.0_dp

    acfs5y = acofs5
    bcfs5y = bcofs5/4.0_dp
    ccfs5y = ccofs5/9.0_dp
    dcfs5y = dcofs5/16.0_dp

    acfs5z = acofs5
    bcfs5z = bcofs5/4.0_dp
    ccfs5z = ccofs5/9.0_dp
    dcfs5z = dcofs5/16.0_dp

!   =========================================================================

!   SECOND CROSS-DERIVATIVES
!   ========================

!   TENTH ORDER EXPLICIT CENTRED DIFFERENCES
    acoefx = acoeff
    bcoefx = bcoeff
    ccoefx = ccoeff
    dcoefx = dcoeff
    ecoefx = ecoeff

    acofxy = acoefx/4.0_dp
    bcofxy = bcoefx/4.0_dp/4.0_dp
    ccofxy = ccoefx/4.0_dp/9.0_dp
    dcofxy = dcoefx/4.0_dp/16.0_dp
    ecofxy = ecoefx/4.0_dp/25.0_dp

    acofxz = acoefx/4.0_dp
    bcofxz = bcoefx/4.0_dp/4.0_dp
    ccofxz = ccoefx/4.0_dp/9.0_dp
    dcofxz = dcoefx/4.0_dp/16.0_dp
    ecofxz = ecoefx/4.0_dp/25.0_dp

    acofyz = acoefx/4.0_dp
    bcofyz = bcoefx/4.0_dp/4.0_dp
    ccofyz = ccoefx/4.0_dp/9.0_dp
    dcofyz = dcoefx/4.0_dp/16.0_dp
    ecofyz = ecoefx/4.0_dp/25.0_dp

!   BOUNDARY TREATMENT
!   ------------------
!   FIRST/SECOND POINT SCHEME (4ND ORDER CENTRED IN TRANSVERSE DIRECTION)
    acofx1 = acoef3/2.0_dp
    bcofx1 = bcoef3/4.0_dp

    acofy1 = acoef3/2.0_dp
    bcofy1 = bcoef3/4.0_dp

    acofz1 = acoef3/2.0_dp
    bcofz1 = bcoef3/4.0_dp

!   FIRST POINT SCHEME (4ND ORDER ONE SIDED/CENTRED)
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

!   SECOND POINT SCHEME (4TH ORDER MIXED/CENTRED)
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

!   THIRD POINT SCHEME (4TH ORDER EXPLICIT CENTRED)
    acf3xy = acoef3/4.0_dp
    bcf3xy = bcoef3/4.0_dp/4.0_dp

    acf3xz = acoef3/4.0_dp
    bcf3xz = bcoef3/4.0_dp/4.0_dp

    acf3yz = acoef3/4.0_dp
    bcf3yz = bcoef3/4.0_dp/4.0_dp

!   FOURTH POINT SCHEME (6TH ORDER EXPLICIT CENTRED)
    acf4xy = acoef4/4.0_dp
    bcf4xy = bcoef4/4.0_dp/4.0_dp
    ccf4xy = ccoef4/4.0_dp/9.0_dp

    acf4xz = acoef4/4.0_dp
    bcf4xz = bcoef4/4.0_dp/4.0_dp
    ccf4xz = ccoef4/4.0_dp/9.0_dp

    acf4yz = acoef4/4.0_dp
    bcf4yz = bcoef4/4.0_dp/4.0_dp
    ccf4yz = ccoef4/4.0_dp/9.0_dp

!   FIFTH POINT SCHEME (8TH ORDER EXPLICIT CENTRED)
    acf5xy = acoef5/4.0_dp
    bcf5xy = bcoef5/4.0_dp/4.0_dp
    ccf5xy = ccoef5/4.0_dp/9.0_dp
    dcf5xy = dcoef5/4.0_dp/16.0_dp

    acf5xz = acoef5/4.0_dp
    bcf5xz = bcoef5/4.0_dp/4.0_dp
    ccf5xz = ccoef5/4.0_dp/9.0_dp
    dcf5xz = dcoef5/4.0_dp/16.0_dp

    acf5yz = acoef5/4.0_dp
    bcf5yz = bcoef5/4.0_dp/4.0_dp
    ccf5yz = ccoef5/4.0_dp/9.0_dp
    dcf5yz = dcoef5/4.0_dp/16.0_dp

!   CORNER POINT SCHEME (4TH ORDER ONE SIDED/ONE SIDED)
    acofc1 = acoef1
    bcofc1 = bcoef1*2.0_dp
    ccofc1 = ccoef1*3.0_dp
    dcofc1 = dcoef1*4.0_dp

    acc1xy = acofc1
    bcc1xy = bcofc1/4.0_dp
    ccc1xy = ccofc1/9.0_dp
    dcc1xy = dcofc1/16.0_dp

    acc1xz = acofc1
    bcc1xz = bcofc1/4.0_dp
    ccc1xz = ccofc1/9.0_dp
    dcc1xz = dcofc1/16.0_dp

    acc1yz = acofc1
    bcc1yz = bcofc1/4.0_dp
    ccc1yz = ccofc1/9.0_dp
    dcc1yz = dcofc1/16.0_dp

!   SECOND CORNER POINT SCHEME (4TH ORDER MIXED)
    acofc2 = -acoef2
    bcofc2 =  bcoef2
    ccofc2 =  ccoef2*2.0_dp
    dcofc2 =  dcoef2*3.0_dp

    acc2xy = acofc2
    bcc2xy = bcofc2
    ccc2xy = ccofc2/4.0_dp
    dcc2xy = dcofc2/9.0_dp

    acc2xz = acofc2
    bcc2xz = bcofc2
    ccc2xz = ccofc2/4.0_dp
    dcc2xz = dcofc2/9.0_dp

    acc2yz = acofc2
    bcc2yz = bcofc2
    ccc2yz = ccofc2/4.0_dp
    dcc2yz = dcofc2/9.0_dp

!   SECOND EDGE POINT SCHEME (4TH ORDER MIXED/CENTRED)
!   USES FIRST AND SECOND POINT COEFFS

!   =========================================================================

!   RSC 10-NOV-2013 INITIALISE MESH STRETCHING
    CALL dfmstr

!   =========================================================================

!   DIFFERENCE COEFFICIENTS FOR WALL BCS
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

!   ONE-SIDED MESH STRETCHING IN Y DIRECTION
!   ACBCYL(1) = ACOEF1*OVDELY*DGDHAT(1)
!   ACBCYL(2) = BCOEF1*OVDELY*DGDHAT(1)
!   ACBCYL(3) = CCOEF1*OVDELY*DGDHAT(1)
!   ACBCYL(4) = DCOEF1*OVDELY*DGDHAT(1)

!   =========================================================================

#ifdef OPS_WITH_CUDAFOR
    acoffx_opsconstant = acoffx
    bcoffx_opsconstant = bcoffx
    ccoffx_opsconstant = ccoffx
    dcoffx_opsconstant = dcoffx
    ecoffx_opsconstant = ecoffx

    acof1x_opsconstant = acof1x
    bcof1x_opsconstant = bcof1x
    ccof1x_opsconstant = ccof1x
    dcof1x_opsconstant = dcof1x

    acof2x_opsconstant = acof2x
    bcof2x_opsconstant = bcof2x
    ccof2x_opsconstant = ccof2x
    dcof2x_opsconstant = dcof2x

    acof3x_opsconstant = acof3x
    bcof3x_opsconstant = bcof3x

    acof4x_opsconstant = acof4x
    bcof4x_opsconstant = bcof4x
    ccof4x_opsconstant = ccof4x

    acof5x_opsconstant = acof5x
    bcof5x_opsconstant = bcof5x
    ccof5x_opsconstant = ccof5x
    dcof5x_opsconstant = dcof5x

    ovdelx_opsconstant = ovdelx

    acofsx_opsconstant = acofsx
    bcofsx_opsconstant = bcofsx
    ccofsx_opsconstant = ccofsx
    dcofsx_opsconstant = dcofsx
    ecofsx_opsconstant = ecofsx

    acfs1x_opsconstant = acfs1x
    bcfs1x_opsconstant = bcfs1x
    ccfs1x_opsconstant = ccfs1x
    dcfs1x_opsconstant = dcfs1x
    ecfs1x_opsconstant = ecfs1x

    acfs2x_opsconstant = acfs2x
    bcfs2x_opsconstant = bcfs2x
    ccfs2x_opsconstant = ccfs2x
    dcfs2x_opsconstant = dcfs2x
    ecfs2x_opsconstant = ecfs2x

    acfs3x_opsconstant = acfs3x
    bcfs3x_opsconstant = bcfs3x

    acfs4x_opsconstant = acfs4x
    bcfs4x_opsconstant = bcfs4x
    ccfs4x_opsconstant = ccfs4x

    acfs5x_opsconstant = acfs5x
    bcfs5x_opsconstant = bcfs5x
    ccfs5x_opsconstant = ccfs5x
    dcfs5x_opsconstant = dcfs5x

    ovdlx2_opsconstant = ovdlx2

    acoffy_opsconstant = acoffy
    bcoffy_opsconstant = bcoffy
    ccoffy_opsconstant = ccoffy
    dcoffy_opsconstant = dcoffy
    ecoffy_opsconstant = ecoffy

    acof1y_opsconstant = acof1y
    bcof1y_opsconstant = bcof1y
    ccof1y_opsconstant = ccof1y
    dcof1y_opsconstant = dcof1y

    acof2y_opsconstant = acof2y
    bcof2y_opsconstant = bcof2y
    ccof2y_opsconstant = ccof2y
    dcof2y_opsconstant = dcof2y

    acof3y_opsconstant = acof3y
    bcof3y_opsconstant = bcof3y

    acof4y_opsconstant = acof4y
    bcof4y_opsconstant = bcof4y
    ccof4y_opsconstant = ccof4y

    acof5y_opsconstant = acof5y
    bcof5y_opsconstant = bcof5y
    ccof5y_opsconstant = ccof5y
    dcof5y_opsconstant = dcof5y

    ovdely_opsconstant = ovdely

    acofsy_opsconstant = acofsy
    bcofsy_opsconstant = bcofsy
    ccofsy_opsconstant = ccofsy
    dcofsy_opsconstant = dcofsy
    ecofsy_opsconstant = ecofsy

    acfs1y_opsconstant = acfs1y
    bcfs1y_opsconstant = bcfs1y
    ccfs1y_opsconstant = ccfs1y
    dcfs1y_opsconstant = dcfs1y
    ecfs1y_opsconstant = ecfs1y

    acfs2y_opsconstant = acfs2y
    bcfs2y_opsconstant = bcfs2y
    ccfs2y_opsconstant = ccfs2y
    dcfs2y_opsconstant = dcfs2y
    ecfs2y_opsconstant = ecfs2y

    acfs3y_opsconstant = acfs3y
    bcfs3y_opsconstant = bcfs3y

    acfs4y_opsconstant = acfs4y
    bcfs4y_opsconstant = bcfs4y
    ccfs4y_opsconstant = ccfs4y

    acfs5y_opsconstant = acfs5y
    bcfs5y_opsconstant = bcfs5y
    ccfs5y_opsconstant = ccfs5y
    dcfs5y_opsconstant = dcfs5y

    ovdly2_opsconstant = ovdly2
    
    acoffz_opsconstant = acoffz
    bcoffz_opsconstant = bcoffz
    ccoffz_opsconstant = ccoffz
    dcoffz_opsconstant = dcoffz
    ecoffz_opsconstant = ecoffz

    acof1z_opsconstant = acof1z
    bcof1z_opsconstant = bcof1z
    ccof1z_opsconstant = ccof1z
    dcof1z_opsconstant = dcof1z

    acof2z_opsconstant = acof2z
    bcof2z_opsconstant = bcof2z
    ccof2z_opsconstant = ccof2z
    dcof2z_opsconstant = dcof2z

    acof3z_opsconstant = acof3z
    bcof3z_opsconstant = bcof3z

    acof4z_opsconstant = acof4z
    bcof4z_opsconstant = bcof4z
    ccof4z_opsconstant = ccof4z

    acof5z_opsconstant = acof5z
    bcof5z_opsconstant = bcof5z
    ccof5z_opsconstant = ccof5z
    dcof5z_opsconstant = dcof5z

    ovdelz_opsconstant = ovdelz

    acofsz_opsconstant = acofsz
    bcofsz_opsconstant = bcofsz
    ccofsz_opsconstant = ccofsz
    dcofsz_opsconstant = dcofsz
    ecofsz_opsconstant = ecofsz

    acfs1z_opsconstant = acfs1z
    bcfs1z_opsconstant = bcfs1z
    ccfs1z_opsconstant = ccfs1z
    dcfs1z_opsconstant = dcfs1z
    ecfs1z_opsconstant = ecfs1z

    acfs2z_opsconstant = acfs2z
    bcfs2z_opsconstant = bcfs2z
    ccfs2z_opsconstant = ccfs2z
    dcfs2z_opsconstant = dcfs2z
    ecfs2z_opsconstant = ecfs2z

    acfs3z_opsconstant = acfs3z
    bcfs3z_opsconstant = bcfs3z

    acfs4z_opsconstant = acfs4z
    bcfs4z_opsconstant = bcfs4z
    ccfs4z_opsconstant = ccfs4z

    acfs5z_opsconstant = acfs5z
    bcfs5z_opsconstant = bcfs5z
    ccfs5z_opsconstant = ccfs5z
    dcfs5z_opsconstant = dcfs5z

    ovdlz2_opsconstant = ovdlz2

    acofx1_opsconstant = acofx1
    bcofx1_opsconstant = bcofx1
    acofy1_opsconstant = acofy1
    bcofy1_opsconstant = bcofy1
    acofz1_opsconstant = acofz1
    bcofz1_opsconstant = bcofz1

    acofxy_opsconstant = acofxy
    bcofxy_opsconstant = bcofxy
    ccofxy_opsconstant = ccofxy
    dcofxy_opsconstant = dcofxy
    ecofxy_opsconstant = ecofxy

    acf1xy_opsconstant = acf1xy
    bcf1xy_opsconstant = bcf1xy
    ccf1xy_opsconstant = ccf1xy
    dcf1xy_opsconstant = dcf1xy

    acf2xy_opsconstant = acf2xy
    bcf2xy_opsconstant = bcf2xy
    ccf2xy_opsconstant = ccf2xy
    dcf2xy_opsconstant = dcf2xy

    acf3xy_opsconstant = acf3xy
    bcf3xy_opsconstant = bcf3xy

    acf4xy_opsconstant = acf4xy
    bcf4xy_opsconstant = bcf4xy
    ccf4xy_opsconstant = ccf4xy

    acf5xy_opsconstant = acf5xy
    bcf5xy_opsconstant = bcf5xy
    ccf5xy_opsconstant = ccf5xy
    dcf5xy_opsconstant = dcf5xy

    acc1xy_opsconstant = acc1xy
    bcc1xy_opsconstant = bcc1xy
    ccc1xy_opsconstant = ccc1xy
    dcc1xy_opsconstant = dcc1xy

    acc2xy_opsconstant = acc2xy
    bcc2xy_opsconstant = bcc2xy
    ccc2xy_opsconstant = ccc2xy
    dcc2xy_opsconstant = dcc2xy

#endif

    call ops_decl_const("acoffx", 1, "double", acoffx)
    call ops_decl_const("bcoffx", 1, "double", bcoffx)
    call ops_decl_const("ccoffx", 1, "double", ccoffx)
    call ops_decl_const("dcoffx", 1, "double", dcoffx)
    call ops_decl_const("ecoffx", 1, "double", ecoffx)

    call ops_decl_const("acof1x", 1, "double", acof1x)
    call ops_decl_const("bcof1x", 1, "double", bcof1x)
    call ops_decl_const("ccof1x", 1, "double", ccof1x)
    call ops_decl_const("dcof1x", 1, "double", dcof1x)

    call ops_decl_const("acof2x", 1, "double", acof2x)
    call ops_decl_const("bcof2x", 1, "double", bcof2x)
    call ops_decl_const("ccof2x", 1, "double", ccof2x)
    call ops_decl_const("dcof2x", 1, "double", dcof2x)

    call ops_decl_const("acof3x", 1, "double", acof3x)
    call ops_decl_const("bcof3x", 1, "double", bcof3x)

    call ops_decl_const("acof4x", 1, "double", acof4x)
    call ops_decl_const("bcof4x", 1, "double", bcof4x)
    call ops_decl_const("ccof4x", 1, "double", ccof4x)

    call ops_decl_const("acof5x", 1, "double", acof5x)
    call ops_decl_const("bcof5x", 1, "double", bcof5x)
    call ops_decl_const("ccof5x", 1, "double", ccof5x)
    call ops_decl_const("dcof5x", 1, "double", dcof5x)

    call ops_decl_const("ovdelx", 1, "double", ovdelx)

    call ops_decl_const("acofsx", 1, "double", acofsx)
    call ops_decl_const("bcofsx", 1, "double", bcofsx)
    call ops_decl_const("ccofsx", 1, "double", ccofsx)
    call ops_decl_const("dcofsx", 1, "double", dcofsx)
    call ops_decl_const("ecofsx", 1, "double", ecofsx)

    call ops_decl_const("acfs1x", 1, "double", acfs1x)
    call ops_decl_const("bcfs1x", 1, "double", bcfs1x)
    call ops_decl_const("ccfs1x", 1, "double", ccfs1x)
    call ops_decl_const("dcfs1x", 1, "double", dcfs1x)
    call ops_decl_const("ecfs1x", 1, "double", ecfs1x)

    call ops_decl_const("acfs2x", 1, "double", acfs2x)
    call ops_decl_const("bcfs2x", 1, "double", bcfs2x)
    call ops_decl_const("ccfs2x", 1, "double", ccfs2x)
    call ops_decl_const("dcfs2x", 1, "double", dcfs2x)
    call ops_decl_const("ecfs2x", 1, "double", ecfs2x)

    call ops_decl_const("acfs3x", 1, "double", acfs3x)
    call ops_decl_const("bcfs3x", 1, "double", bcfs3x)

    call ops_decl_const("acfs4x", 1, "double", acfs4x)
    call ops_decl_const("bcfs4x", 1, "double", bcfs4x)
    call ops_decl_const("ccfs4x", 1, "double", ccfs4x)

    call ops_decl_const("acfs5x", 1, "double", acfs5x)
    call ops_decl_const("bcfs5x", 1, "double", bcfs5x)
    call ops_decl_const("ccfs5x", 1, "double", ccfs5x)
    call ops_decl_const("dcfs5x", 1, "double", dcfs5x)

    call ops_decl_const("ovdlx2", 1, "double", ovdlx2)

    call ops_decl_const("acoffy", 1, "double", acoffy)
    call ops_decl_const("bcoffy", 1, "double", bcoffy)
    call ops_decl_const("ccoffy", 1, "double", ccoffy)
    call ops_decl_const("dcoffy", 1, "double", dcoffy)
    call ops_decl_const("ecoffy", 1, "double", ecoffy)

    call ops_decl_const("acof1y", 1, "double", acof1y)
    call ops_decl_const("bcof1y", 1, "double", bcof1y)
    call ops_decl_const("ccof1y", 1, "double", ccof1y)
    call ops_decl_const("dcof1y", 1, "double", dcof1y)

    call ops_decl_const("acof2y", 1, "double", acof2y)
    call ops_decl_const("bcof2y", 1, "double", bcof2y)
    call ops_decl_const("ccof2y", 1, "double", ccof2y)
    call ops_decl_const("dcof2y", 1, "double", dcof2y)

    call ops_decl_const("acof3y", 1, "double", acof3y)
    call ops_decl_const("bcof3y", 1, "double", bcof3y)

    call ops_decl_const("acof4y", 1, "double", acof4y)
    call ops_decl_const("bcof4y", 1, "double", bcof4y)
    call ops_decl_const("ccof4y", 1, "double", ccof4y)

    call ops_decl_const("acof5y", 1, "double", acof5y)
    call ops_decl_const("bcof5y", 1, "double", bcof5y)
    call ops_decl_const("ccof5y", 1, "double", ccof5y)
    call ops_decl_const("dcof5y", 1, "double", dcof5y)

    call ops_decl_const("ovdely", 1, "double", ovdely)

    call ops_decl_const("acofsy", 1, "double", acofsy)
    call ops_decl_const("bcofsy", 1, "double", bcofsy)
    call ops_decl_const("ccofsy", 1, "double", ccofsy)
    call ops_decl_const("dcofsy", 1, "double", dcofsy)
    call ops_decl_const("ecofsy", 1, "double", ecofsy)

    call ops_decl_const("acfs1y", 1, "double", acfs1y)
    call ops_decl_const("bcfs1y", 1, "double", bcfs1y)
    call ops_decl_const("ccfs1y", 1, "double", ccfs1y)
    call ops_decl_const("dcfs1y", 1, "double", dcfs1y)
    call ops_decl_const("ecfs1y", 1, "double", ecfs1y)

    call ops_decl_const("acfs2y", 1, "double", acfs2y)
    call ops_decl_const("bcfs2y", 1, "double", bcfs2y)
    call ops_decl_const("ccfs2y", 1, "double", ccfs2y)
    call ops_decl_const("dcfs2y", 1, "double", dcfs2y)
    call ops_decl_const("ecfs2y", 1, "double", ecfs2y)

    call ops_decl_const("acfs3y", 1, "double", acfs3y)
    call ops_decl_const("bcfs3y", 1, "double", bcfs3y)

    call ops_decl_const("acfs4y", 1, "double", acfs4y)
    call ops_decl_const("bcfs4y", 1, "double", bcfs4y)
    call ops_decl_const("ccfs4y", 1, "double", ccfs4y)

    call ops_decl_const("acfs5y", 1, "double", acfs5y)
    call ops_decl_const("bcfs5y", 1, "double", bcfs5y)
    call ops_decl_const("ccfs5y", 1, "double", ccfs5y)
    call ops_decl_const("dcfs5y", 1, "double", dcfs5y)

    call ops_decl_const("ovdly2", 1, "double", ovdly2)

    call ops_decl_const("acoffz", 1, "double", acoffz)
    call ops_decl_const("bcoffz", 1, "double", bcoffz)
    call ops_decl_const("ccoffz", 1, "double", ccoffz)
    call ops_decl_const("dcoffz", 1, "double", dcoffz)
    call ops_decl_const("ecoffz", 1, "double", ecoffz)

    call ops_decl_const("acof1z", 1, "double", acof1z)
    call ops_decl_const("bcof1z", 1, "double", bcof1z)
    call ops_decl_const("ccof1z", 1, "double", ccof1z)
    call ops_decl_const("dcof1z", 1, "double", dcof1z)

    call ops_decl_const("acof2z", 1, "double", acof2z)
    call ops_decl_const("bcof2z", 1, "double", bcof2z)
    call ops_decl_const("ccof2z", 1, "double", ccof2z)
    call ops_decl_const("dcof2z", 1, "double", dcof2z)

    call ops_decl_const("acof3z", 1, "double", acof3z)
    call ops_decl_const("bcof3z", 1, "double", bcof3z)

    call ops_decl_const("acof4z", 1, "double", acof4z)
    call ops_decl_const("bcof4z", 1, "double", bcof4z)
    call ops_decl_const("ccof4z", 1, "double", ccof4z)

    call ops_decl_const("acof5z", 1, "double", acof5z)
    call ops_decl_const("bcof5z", 1, "double", bcof5z)
    call ops_decl_const("ccof5z", 1, "double", ccof5z)
    call ops_decl_const("dcof5z", 1, "double", dcof5z)

    call ops_decl_const("ovdelz", 1, "double", ovdelz)

    call ops_decl_const("acofsz", 1, "double", acofsz)
    call ops_decl_const("bcofsz", 1, "double", bcofsz)
    call ops_decl_const("ccofsz", 1, "double", ccofsz)
    call ops_decl_const("dcofsz", 1, "double", dcofsz)
    call ops_decl_const("ecofsz", 1, "double", ecofsz)

    call ops_decl_const("acfs1z", 1, "double", acfs1z)
    call ops_decl_const("bcfs1z", 1, "double", bcfs1z)
    call ops_decl_const("ccfs1z", 1, "double", ccfs1z)
    call ops_decl_const("dcfs1z", 1, "double", dcfs1z)
    call ops_decl_const("ecfs1z", 1, "double", ecfs1z)

    call ops_decl_const("acfs2z", 1, "double", acfs2z)
    call ops_decl_const("bcfs2z", 1, "double", bcfs2z)
    call ops_decl_const("ccfs2z", 1, "double", ccfs2z)
    call ops_decl_const("dcfs2z", 1, "double", dcfs2z)
    call ops_decl_const("ecfs2z", 1, "double", ecfs2z)

    call ops_decl_const("acfs3z", 1, "double", acfs3z)
    call ops_decl_const("bcfs3z", 1, "double", bcfs3z)

    call ops_decl_const("acfs4z", 1, "double", acfs4z)
    call ops_decl_const("bcfs4z", 1, "double", bcfs4z)
    call ops_decl_const("ccfs4z", 1, "double", ccfs4z)

    call ops_decl_const("acfs5z", 1, "double", acfs5z)
    call ops_decl_const("bcfs5z", 1, "double", bcfs5z)
    call ops_decl_const("ccfs5z", 1, "double", ccfs5z)
    call ops_decl_const("dcfs5z", 1, "double", dcfs5z)

    call ops_decl_const("ovdlz2", 1, "double", ovdlz2)

    call ops_decl_const("acofx1", 1, "double", acofx1)
    call ops_decl_const("bcofx1", 1, "double", bcofx1)
    call ops_decl_const("acofy1", 1, "double", acofy1)
    call ops_decl_const("bcofy1", 1, "double", bcofy1)
    call ops_decl_const("acofz1", 1, "double", acofz1)
    call ops_decl_const("bcofz1", 1, "double", bcofz1)

    call ops_decl_const("acofxy", 1, "double", acofxy)
    call ops_decl_const("bcofxy", 1, "double", bcofxy)
    call ops_decl_const("ccofxy", 1, "double", ccofxy)
    call ops_decl_const("dcofxy", 1, "double", dcofxy)
    call ops_decl_const("ecofxy", 1, "double", ecofxy)

    call ops_decl_const("acf1xy", 1, "double", acf1xy)
    call ops_decl_const("bcf1xy", 1, "double", bcf1xy)
    call ops_decl_const("ccf1xy", 1, "double", ccf1xy)
    call ops_decl_const("dcf1xy", 1, "double", dcf1xy)

    call ops_decl_const("acf2xy", 1, "double", acf2xy)
    call ops_decl_const("bcf2xy", 1, "double", bcf2xy)
    call ops_decl_const("ccf2xy", 1, "double", ccf2xy)
    call ops_decl_const("dcf2xy", 1, "double", dcf2xy)

    call ops_decl_const("acf3xy", 1, "double", acf3xy)
    call ops_decl_const("bcf3xy", 1, "double", bcf3xy)

    call ops_decl_const("acf4xy", 1, "double", acf4xy)
    call ops_decl_const("bcf4xy", 1, "double", bcf4xy)
    call ops_decl_const("ccf4xy", 1, "double", ccf4xy)

    call ops_decl_const("acf5xy", 1, "double", acf5xy)
    call ops_decl_const("bcf5xy", 1, "double", bcf5xy)
    call ops_decl_const("ccf5xy", 1, "double", ccf5xy)
    call ops_decl_const("dcf5xy", 1, "double", dcf5xy)

    call ops_decl_const("acc1xy", 1, "double", acc1xy)
    call ops_decl_const("bcc1xy", 1, "double", bcc1xy)
    call ops_decl_const("ccc1xy", 1, "double", ccc1xy)
    call ops_decl_const("dcc1xy", 1, "double", dcc1xy)

    call ops_decl_const("acc2xy", 1, "double", acc2xy)
    call ops_decl_const("bcc2xy", 1, "double", bcc2xy)
    call ops_decl_const("ccc2xy", 1, "double", ccc2xy)
    call ops_decl_const("dcc2xy", 1, "double", dcc2xy)

END SUBROUTINE dfinit

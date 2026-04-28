
! Auto-generated at 2026-04-28 18:43:09.036281 by ops-translator
SUBROUTINE dfinit

  USE ops_fortran_declarations
  USE ops_fortran_rt_support

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

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
  deltax = xgdlen / ((nxglbl - 1) * one)
  ovdelx = one / deltax
  ovdlx2 = ovdelx * ovdelx
  deltay = ygdlen / ((nyglbl - 1) * one)
  ovdely = one / deltay
  ovdly2 = ovdely * ovdely
  deltaz = zgdlen / ((nzglbl - 1) * one)
  ovdelz = one / deltaz
  ovdlz2 = ovdelz * ovdelz

  !   =========================================================================

  !   FIRST DERIVATIVES
  !   =================

  !   INTERIOR SCHEME
  !   ---------------

  !   TENTH ORDER EXPLICIT CENTRED DIFFERENCES
  acoeff = 5.0_8 / 3.0_8
  bcoeff = - 20.0_8 / 21.0_8
  ccoeff = 5.0_8 / 14.0_8
  dcoeff = - 5.0_8 / 63.0_8
  ecoeff = 1.0_8 / 126.0_8

  acoffx = acoeff / 2.0_8
  bcoffx = bcoeff / 4.0_8
  ccoffx = ccoeff / 6.0_8
  dcoffx = dcoeff / 8.0_8
  ecoffx = ecoeff / 10.0_8

  acoffy = acoeff / 2.0_8
  bcoffy = bcoeff / 4.0_8
  ccoffy = ccoeff / 6.0_8
  dcoffy = dcoeff / 8.0_8
  ecoffy = ecoeff / 10.0_8

  acoffz = acoeff / 2.0_8
  bcoffz = bcoeff / 4.0_8
  ccoffz = ccoeff / 6.0_8
  dcoffz = dcoeff / 8.0_8
  ecoffz = ecoeff / 10.0_8

  !   BOUNDARY TREATMENT
  !   ------------------

  !   FIRST POINT SCHEME (4TH ORDER ONE SIDED)
  acoef1 = 4.0_8
  bcoef1 = - 3.0_8
  ccoef1 = 4.0_8 / 3.0_8
  dcoef1 = - 1.0_8 / 4.0_8

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
  acoef2 = - 1.0_8 / 4.0_8
  bcoef2 = 3.0_8 / 2.0_8
  ccoef2 = - 1.0_8 / 2.0_8
  dcoef2 = 1.0_8 / 12.0_8

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
  acoef3 = 4.0_8 / 3.0_8
  bcoef3 = - 1.0_8 / 3.0_8

  acof3x = acoef3 / 2.0_8
  bcof3x = bcoef3 / 4.0_8

  acof3y = acoef3 / 2.0_8
  bcof3y = bcoef3 / 4.0_8

  acof3z = acoef3 / 2.0_8
  bcof3z = bcoef3 / 4.0_8

  !   4TH POINT SCHEME (6TH ORDER EXPLICIT CENTRED)
  acoef4 = 3.0_8 / 2.0_8
  bcoef4 = - 3.0_8 / 5.0_8
  ccoef4 = 1.0_8 / 10.0_8

  acof4x = acoef4 / 2.0_8
  bcof4x = bcoef4 / 4.0_8
  ccof4x = ccoef4 / 6.0_8

  acof4y = acoef4 / 2.0_8
  bcof4y = bcoef4 / 4.0_8
  ccof4y = ccoef4 / 6.0_8

  acof4z = acoef4 / 2.0_8
  bcof4z = bcoef4 / 4.0_8
  ccof4z = ccoef4 / 6.0_8

  !   5TH POINT SCHEME (8TH ORDER EXPLICIT CENTRED)
  acoef5 = 8.0_8 / 5.0_8
  bcoef5 = - 4.0_8 / 5.0_8
  ccoef5 = 8.0_8 / 35.0_8
  dcoef5 = - 1.0_8 / 35.0_8

  acof5x = acoef5 / 2.0_8
  bcof5x = bcoef5 / 4.0_8
  ccof5x = ccoef5 / 6.0_8
  dcof5x = dcoef5 / 8.0_8

  acof5y = acoef5 / 2.0_8
  bcof5y = bcoef5 / 4.0_8
  ccof5y = ccoef5 / 6.0_8
  dcof5y = dcoef5 / 8.0_8

  acof5z = acoef5 / 2.0_8
  bcof5z = bcoef5 / 4.0_8
  ccof5z = ccoef5 / 6.0_8
  dcof5z = dcoef5 / 8.0_8

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
  bcofsx = bcoefs / 4.0_8
  ccofsx = ccoefs / 9.0_8
  dcofsx = dcoefs / 16.0_8
  ecofsx = ecoefs / 25.0_8

  acofsy = acoefs
  bcofsy = bcoefs / 4.0_8
  ccofsy = ccoefs / 9.0_8
  dcofsy = dcoefs / 16.0_8
  ecofsy = ecoefs / 25.0_8

  acofsz = acoefs
  bcofsz = bcoefs / 4.0_8
  ccofsz = ccoefs / 9.0_8
  dcofsz = dcoefs / 16.0_8
  ecofsz = ecoefs / 25.0_8

  !   BOUNDARY TREATMENT
  !   ------------------

  !   FIRST POINT SCHEME (4TH ORDER ONE SIDED)
  acofs1 = - 77.0_8 / 6.0_8
  bcofs1 = 107.0_8 / 6.0_8
  ccofs1 = - 13.0_8
  dcofs1 = 61.0_8 / 12.0_8
  ecofs1 = - 5.0_8 / 6.0_8

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
  acofs2 = 5.0_8 / 6.0_8
  bcofs2 = - 1.0_8 / 3.0_8
  ccofs2 = 7.0_8 / 6.0_8
  dcofs2 = - 1.0_8 / 2.0_8
  ecofs2 = 1.0_8 / 12.0_8

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
  bcfs3x = bcofs3 / 4.0_8

  acfs3y = acofs3
  bcfs3y = bcofs3 / 4.0_8

  acfs3z = acofs3
  bcfs3z = bcofs3 / 4.0_8

  !   4TH POINT SCHEME (6TH ORDER EXPLICIT CENTRED)
  acofs4 = acoef4
  bcofs4 = bcoef4
  ccofs4 = ccoef4

  acfs4x = acofs4
  bcfs4x = bcofs4 / 4.0_8
  ccfs4x = ccofs4 / 9.0_8

  acfs4y = acofs4
  bcfs4y = bcofs4 / 4.0_8
  ccfs4y = ccofs4 / 9.0_8

  acfs4z = acofs4
  bcfs4z = bcofs4 / 4.0_8
  ccfs4z = ccofs4 / 9.0_8

  !   5TH POINT SCHEME (8TH ORDER EXPLICIT CENTRED)
  acofs5 = acoef5
  bcofs5 = bcoef5
  ccofs5 = ccoef5
  dcofs5 = dcoef5

  acfs5x = acofs5
  bcfs5x = bcofs5 / 4.0_8
  ccfs5x = ccofs5 / 9.0_8
  dcfs5x = dcofs5 / 16.0_8

  acfs5y = acofs5
  bcfs5y = bcofs5 / 4.0_8
  ccfs5y = ccofs5 / 9.0_8
  dcfs5y = dcofs5 / 16.0_8

  acfs5z = acofs5
  bcfs5z = bcofs5 / 4.0_8
  ccfs5z = ccofs5 / 9.0_8
  dcfs5z = dcofs5 / 16.0_8

  !   =========================================================================

  !   SECOND CROSS-DERIVATIVES
  !   ========================

  !   TENTH ORDER EXPLICIT CENTRED DIFFERENCES
  acoefx = acoeff
  bcoefx = bcoeff
  ccoefx = ccoeff
  dcoefx = dcoeff
  ecoefx = ecoeff

  acofxy = acoefx / 4.0_8
  bcofxy = bcoefx / 4.0_8 / 4.0_8
  ccofxy = ccoefx / 4.0_8 / 9.0_8
  dcofxy = dcoefx / 4.0_8 / 16.0_8
  ecofxy = ecoefx / 4.0_8 / 25.0_8

  acofxz = acoefx / 4.0_8
  bcofxz = bcoefx / 4.0_8 / 4.0_8
  ccofxz = ccoefx / 4.0_8 / 9.0_8
  dcofxz = dcoefx / 4.0_8 / 16.0_8
  ecofxz = ecoefx / 4.0_8 / 25.0_8

  acofyz = acoefx / 4.0_8
  bcofyz = bcoefx / 4.0_8 / 4.0_8
  ccofyz = ccoefx / 4.0_8 / 9.0_8
  dcofyz = dcoefx / 4.0_8 / 16.0_8
  ecofyz = ecoefx / 4.0_8 / 25.0_8

  !   BOUNDARY TREATMENT
  !   ------------------
  !   FIRST/SECOND POINT SCHEME (4ND ORDER CENTRED IN TRANSVERSE DIRECTION)
  acofx1 = acoef3 / 2.0_8
  bcofx1 = bcoef3 / 4.0_8

  acofy1 = acoef3 / 2.0_8
  bcofy1 = bcoef3 / 4.0_8

  acofz1 = acoef3 / 2.0_8
  bcofz1 = bcoef3 / 4.0_8

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
  acf3xy = acoef3 / 4.0_8
  bcf3xy = bcoef3 / 4.0_8 / 4.0_8

  acf3xz = acoef3 / 4.0_8
  bcf3xz = bcoef3 / 4.0_8 / 4.0_8

  acf3yz = acoef3 / 4.0_8
  bcf3yz = bcoef3 / 4.0_8 / 4.0_8

  !   FOURTH POINT SCHEME (6TH ORDER EXPLICIT CENTRED)
  acf4xy = acoef4 / 4.0_8
  bcf4xy = bcoef4 / 4.0_8 / 4.0_8
  ccf4xy = ccoef4 / 4.0_8 / 9.0_8

  acf4xz = acoef4 / 4.0_8
  bcf4xz = bcoef4 / 4.0_8 / 4.0_8
  ccf4xz = ccoef4 / 4.0_8 / 9.0_8

  acf4yz = acoef4 / 4.0_8
  bcf4yz = bcoef4 / 4.0_8 / 4.0_8
  ccf4yz = ccoef4 / 4.0_8 / 9.0_8

  !   FIFTH POINT SCHEME (8TH ORDER EXPLICIT CENTRED)
  acf5xy = acoef5 / 4.0_8
  bcf5xy = bcoef5 / 4.0_8 / 4.0_8
  ccf5xy = ccoef5 / 4.0_8 / 9.0_8
  dcf5xy = dcoef5 / 4.0_8 / 16.0_8

  acf5xz = acoef5 / 4.0_8
  bcf5xz = bcoef5 / 4.0_8 / 4.0_8
  ccf5xz = ccoef5 / 4.0_8 / 9.0_8
  dcf5xz = dcoef5 / 4.0_8 / 16.0_8

  acf5yz = acoef5 / 4.0_8
  bcf5yz = bcoef5 / 4.0_8 / 4.0_8
  ccf5yz = ccoef5 / 4.0_8 / 9.0_8
  dcf5yz = dcoef5 / 4.0_8 / 16.0_8

  !   CORNER POINT SCHEME (4TH ORDER ONE SIDED/ONE SIDED)
  acofc1 = acoef1
  bcofc1 = bcoef1 * 2.0_8
  ccofc1 = ccoef1 * 3.0_8
  dcofc1 = dcoef1 * 4.0_8

  acc1xy = acofc1
  bcc1xy = bcofc1 / 4.0_8
  ccc1xy = ccofc1 / 9.0_8
  dcc1xy = dcofc1 / 16.0_8

  acc1xz = acofc1
  bcc1xz = bcofc1 / 4.0_8
  ccc1xz = ccofc1 / 9.0_8
  dcc1xz = dcofc1 / 16.0_8

  acc1yz = acofc1
  bcc1yz = bcofc1 / 4.0_8
  ccc1yz = ccofc1 / 9.0_8
  dcc1yz = dcofc1 / 16.0_8

  !   SECOND CORNER POINT SCHEME (4TH ORDER MIXED)
  acofc2 = - acoef2
  bcofc2 = bcoef2
  ccofc2 = ccoef2 * 2.0_8
  dcofc2 = dcoef2 * 3.0_8

  acc2xy = acofc2
  bcc2xy = bcofc2
  ccc2xy = ccofc2 / 4.0_8
  dcc2xy = dcofc2 / 9.0_8

  acc2xz = acofc2
  bcc2xz = bcofc2
  ccc2xz = ccofc2 / 4.0_8
  dcc2xz = dcofc2 / 9.0_8

  acc2yz = acofc2
  bcc2yz = bcofc2
  ccc2yz = ccofc2 / 4.0_8
  dcc2yz = dcofc2 / 9.0_8

  !   SECOND EDGE POINT SCHEME (4TH ORDER MIXED/CENTRED)
  !   USES FIRST AND SECOND POINT COEFFS

  !   =========================================================================

  !   RSC 10-NOV-2013 INITIALISE MESH STRETCHING
  CALL dfmstr

  !   =========================================================================

  !   DIFFERENCE COEFFICIENTS FOR WALL BCS
  acbcxl(1) = acoef1 * ovdelx
  acbcxl(2) = bcoef1 * ovdelx
  acbcxl(3) = ccoef1 * ovdelx
  acbcxl(4) = dcoef1 * ovdelx

  acbcxr(1) = acoef1 * ovdelx
  acbcxr(2) = bcoef1 * ovdelx
  acbcxr(3) = ccoef1 * ovdelx
  acbcxr(4) = dcoef1 * ovdelx

  acbcyl(1) = acoef1 * ovdely
  acbcyl(2) = bcoef1 * ovdely
  acbcyl(3) = ccoef1 * ovdely
  acbcyl(4) = dcoef1 * ovdely

  acbcyr(1) = acoef1 * ovdely
  acbcyr(2) = bcoef1 * ovdely
  acbcyr(3) = ccoef1 * ovdely
  acbcyr(4) = dcoef1 * ovdely

  acbczl(1) = acoef1 * ovdelz
  acbczl(2) = bcoef1 * ovdelz
  acbczl(3) = ccoef1 * ovdelz
  acbczl(4) = dcoef1 * ovdelz

  acbczr(1) = acoef1 * ovdelz
  acbczr(2) = bcoef1 * ovdelz
  acbczr(3) = ccoef1 * ovdelz
  acbczr(4) = dcoef1 * ovdelz

  !   ONE-SIDED MESH STRETCHING IN Y DIRECTION
  !   ACBCYL(1) = ACOEF1*OVDELY*DGDHAT(1)
  !   ACBCYL(2) = BCOEF1*OVDELY*DGDHAT(1)
  !   ACBCYL(3) = CCOEF1*OVDELY*DGDHAT(1)
  !   ACBCYL(4) = DCOEF1*OVDELY*DGDHAT(1)

  !   ==============================================================
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acoffx", 1, "real(kind=8)", acoffx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcoffx", 1, "real(kind=8)", bcoffx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccoffx", 1, "real(kind=8)", ccoffx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcoffx", 1, "real(kind=8)", dcoffx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecoffx", 1, "real(kind=8)", ecoffx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof1x", 1, "real(kind=8)", acof1x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof1x", 1, "real(kind=8)", bcof1x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof1x", 1, "real(kind=8)", ccof1x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcof1x", 1, "real(kind=8)", dcof1x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof2x", 1, "real(kind=8)", acof2x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof2x", 1, "real(kind=8)", bcof2x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof2x", 1, "real(kind=8)", ccof2x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcof2x", 1, "real(kind=8)", dcof2x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof3x", 1, "real(kind=8)", acof3x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof3x", 1, "real(kind=8)", bcof3x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof4x", 1, "real(kind=8)", acof4x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof4x", 1, "real(kind=8)", bcof4x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof4x", 1, "real(kind=8)", ccof4x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof5x", 1, "real(kind=8)", acof5x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof5x", 1, "real(kind=8)", bcof5x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof5x", 1, "real(kind=8)", ccof5x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcof5x", 1, "real(kind=8)", dcof5x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ovdelx", 1, "real(kind=8)", ovdelx)
#endif
!   ==============================================================
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acofsx", 1, "real(kind=8)", acofsx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcofsx", 1, "real(kind=8)", bcofsx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccofsx", 1, "real(kind=8)", ccofsx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcofsx", 1, "real(kind=8)", dcofsx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecofsx", 1, "real(kind=8)", ecofsx)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs1x", 1, "real(kind=8)", acfs1x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs1x", 1, "real(kind=8)", bcfs1x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs1x", 1, "real(kind=8)", ccfs1x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcfs1x", 1, "real(kind=8)", dcfs1x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecfs1x", 1, "real(kind=8)", ecfs1x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs2x", 1, "real(kind=8)", acfs2x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs2x", 1, "real(kind=8)", bcfs2x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs2x", 1, "real(kind=8)", ccfs2x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcfs2x", 1, "real(kind=8)", dcfs2x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecfs2x", 1, "real(kind=8)", ecfs2x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs3x", 1, "real(kind=8)", acfs3x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs3x", 1, "real(kind=8)", bcfs3x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs4x", 1, "real(kind=8)", acfs4x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs4x", 1, "real(kind=8)", bcfs4x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs4x", 1, "real(kind=8)", ccfs4x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs5x", 1, "real(kind=8)", acfs5x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs5x", 1, "real(kind=8)", bcfs5x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs5x", 1, "real(kind=8)", ccfs5x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcfs5x", 1, "real(kind=8)", dcfs5x)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ovdlx2", 1, "real(kind=8)", ovdlx2)
#endif
!   ==============================================================
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acoffy", 1, "real(kind=8)", acoffy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcoffy", 1, "real(kind=8)", bcoffy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccoffy", 1, "real(kind=8)", ccoffy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcoffy", 1, "real(kind=8)", dcoffy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecoffy", 1, "real(kind=8)", ecoffy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof1y", 1, "real(kind=8)", acof1y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof1y", 1, "real(kind=8)", bcof1y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof1y", 1, "real(kind=8)", ccof1y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcof1y", 1, "real(kind=8)", dcof1y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof2y", 1, "real(kind=8)", acof2y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof2y", 1, "real(kind=8)", bcof2y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof2y", 1, "real(kind=8)", ccof2y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcof2y", 1, "real(kind=8)", dcof2y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof3y", 1, "real(kind=8)", acof3y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof3y", 1, "real(kind=8)", bcof3y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof4y", 1, "real(kind=8)", acof4y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof4y", 1, "real(kind=8)", bcof4y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof4y", 1, "real(kind=8)", ccof4y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof5y", 1, "real(kind=8)", acof5y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof5y", 1, "real(kind=8)", bcof5y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof5y", 1, "real(kind=8)", ccof5y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcof5y", 1, "real(kind=8)", dcof5y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ovdely", 1, "real(kind=8)", ovdely)
#endif
!   ==============================================================
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acofsy", 1, "real(kind=8)", acofsy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcofsy", 1, "real(kind=8)", bcofsy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccofsy", 1, "real(kind=8)", ccofsy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcofsy", 1, "real(kind=8)", dcofsy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecofsy", 1, "real(kind=8)", ecofsy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs1y", 1, "real(kind=8)", acfs1y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs1y", 1, "real(kind=8)", bcfs1y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs1y", 1, "real(kind=8)", ccfs1y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcfs1y", 1, "real(kind=8)", dcfs1y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecfs1y", 1, "real(kind=8)", ecfs1y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs2y", 1, "real(kind=8)", acfs2y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs2y", 1, "real(kind=8)", bcfs2y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs2y", 1, "real(kind=8)", ccfs2y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcfs2y", 1, "real(kind=8)", dcfs2y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecfs2y", 1, "real(kind=8)", ecfs2y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs3y", 1, "real(kind=8)", acfs3y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs3y", 1, "real(kind=8)", bcfs3y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs4y", 1, "real(kind=8)", acfs4y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs4y", 1, "real(kind=8)", bcfs4y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs4y", 1, "real(kind=8)", ccfs4y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs5y", 1, "real(kind=8)", acfs5y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs5y", 1, "real(kind=8)", bcfs5y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs5y", 1, "real(kind=8)", ccfs5y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcfs5y", 1, "real(kind=8)", dcfs5y)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ovdly2", 1, "real(kind=8)", ovdly2)
#endif
!   ==============================================================
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acoffz", 1, "real(kind=8)", acoffz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcoffz", 1, "real(kind=8)", bcoffz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccoffz", 1, "real(kind=8)", ccoffz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcoffz", 1, "real(kind=8)", dcoffz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecoffz", 1, "real(kind=8)", ecoffz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof1z", 1, "real(kind=8)", acof1z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof1z", 1, "real(kind=8)", bcof1z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof1z", 1, "real(kind=8)", ccof1z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcof1z", 1, "real(kind=8)", dcof1z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof2z", 1, "real(kind=8)", acof2z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof2z", 1, "real(kind=8)", bcof2z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof2z", 1, "real(kind=8)", ccof2z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcof2z", 1, "real(kind=8)", dcof2z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof3z", 1, "real(kind=8)", acof3z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof3z", 1, "real(kind=8)", bcof3z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof4z", 1, "real(kind=8)", acof4z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof4z", 1, "real(kind=8)", bcof4z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof4z", 1, "real(kind=8)", ccof4z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acof5z", 1, "real(kind=8)", acof5z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcof5z", 1, "real(kind=8)", bcof5z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccof5z", 1, "real(kind=8)", ccof5z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcof5z", 1, "real(kind=8)", dcof5z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ovdelz", 1, "real(kind=8)", ovdelz)
#endif
!   ==============================================================
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acofsz", 1, "real(kind=8)", acofsz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcofsz", 1, "real(kind=8)", bcofsz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccofsz", 1, "real(kind=8)", ccofsz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcofsz", 1, "real(kind=8)", dcofsz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecofsz", 1, "real(kind=8)", ecofsz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs1z", 1, "real(kind=8)", acfs1z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs1z", 1, "real(kind=8)", bcfs1z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs1z", 1, "real(kind=8)", ccfs1z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcfs1z", 1, "real(kind=8)", dcfs1z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecfs1z", 1, "real(kind=8)", ecfs1z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs2z", 1, "real(kind=8)", acfs2z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs2z", 1, "real(kind=8)", bcfs2z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs2z", 1, "real(kind=8)", ccfs2z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcfs2z", 1, "real(kind=8)", dcfs2z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecfs2z", 1, "real(kind=8)", ecfs2z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs3z", 1, "real(kind=8)", acfs3z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs3z", 1, "real(kind=8)", bcfs3z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs4z", 1, "real(kind=8)", acfs4z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs4z", 1, "real(kind=8)", bcfs4z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs4z", 1, "real(kind=8)", ccfs4z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acfs5z", 1, "real(kind=8)", acfs5z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcfs5z", 1, "real(kind=8)", bcfs5z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccfs5z", 1, "real(kind=8)", ccfs5z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcfs5z", 1, "real(kind=8)", dcfs5z)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ovdlz2", 1, "real(kind=8)", ovdlz2)
#endif
!   ==============================================================
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acofx1", 1, "real(kind=8)", acofx1)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcofx1", 1, "real(kind=8)", bcofx1)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acofy1", 1, "real(kind=8)", acofy1)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcofy1", 1, "real(kind=8)", bcofy1)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acofz1", 1, "real(kind=8)", acofz1)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcofz1", 1, "real(kind=8)", bcofz1)
#endif
!   ==============================================================
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acofxy", 1, "real(kind=8)", acofxy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcofxy", 1, "real(kind=8)", bcofxy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccofxy", 1, "real(kind=8)", ccofxy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcofxy", 1, "real(kind=8)", dcofxy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecofxy", 1, "real(kind=8)", ecofxy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf1xy", 1, "real(kind=8)", acf1xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf1xy", 1, "real(kind=8)", bcf1xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf1xy", 1, "real(kind=8)", ccf1xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcf1xy", 1, "real(kind=8)", dcf1xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf2xy", 1, "real(kind=8)", acf2xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf2xy", 1, "real(kind=8)", bcf2xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf2xy", 1, "real(kind=8)", ccf2xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcf2xy", 1, "real(kind=8)", dcf2xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf3xy", 1, "real(kind=8)", acf3xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf3xy", 1, "real(kind=8)", bcf3xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf4xy", 1, "real(kind=8)", acf4xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf4xy", 1, "real(kind=8)", bcf4xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf4xy", 1, "real(kind=8)", ccf4xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf5xy", 1, "real(kind=8)", acf5xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf5xy", 1, "real(kind=8)", bcf5xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf5xy", 1, "real(kind=8)", ccf5xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcf5xy", 1, "real(kind=8)", dcf5xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acc1xy", 1, "real(kind=8)", acc1xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcc1xy", 1, "real(kind=8)", bcc1xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccc1xy", 1, "real(kind=8)", ccc1xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcc1xy", 1, "real(kind=8)", dcc1xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acc2xy", 1, "real(kind=8)", acc2xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcc2xy", 1, "real(kind=8)", bcc2xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccc2xy", 1, "real(kind=8)", ccc2xy)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcc2xy", 1, "real(kind=8)", dcc2xy)
#endif
!   ==============================================================
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acofxz", 1, "real(kind=8)", acofxz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcofxz", 1, "real(kind=8)", bcofxz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccofxz", 1, "real(kind=8)", ccofxz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcofxz", 1, "real(kind=8)", dcofxz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecofxz", 1, "real(kind=8)", ecofxz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf1xz", 1, "real(kind=8)", acf1xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf1xz", 1, "real(kind=8)", bcf1xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf1xz", 1, "real(kind=8)", ccf1xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcf1xz", 1, "real(kind=8)", dcf1xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf2xz", 1, "real(kind=8)", acf2xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf2xz", 1, "real(kind=8)", bcf2xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf2xz", 1, "real(kind=8)", ccf2xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcf2xz", 1, "real(kind=8)", dcf2xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf3xz", 1, "real(kind=8)", acf3xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf3xz", 1, "real(kind=8)", bcf3xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf4xz", 1, "real(kind=8)", acf4xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf4xz", 1, "real(kind=8)", bcf4xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf4xz", 1, "real(kind=8)", ccf4xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf5xz", 1, "real(kind=8)", acf5xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf5xz", 1, "real(kind=8)", bcf5xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf5xz", 1, "real(kind=8)", ccf5xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcf5xz", 1, "real(kind=8)", dcf5xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acc1xz", 1, "real(kind=8)", acc1xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcc1xz", 1, "real(kind=8)", bcc1xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccc1xz", 1, "real(kind=8)", ccc1xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcc1xz", 1, "real(kind=8)", dcc1xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acc2xz", 1, "real(kind=8)", acc2xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcc2xz", 1, "real(kind=8)", bcc2xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccc2xz", 1, "real(kind=8)", ccc2xz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcc2xz", 1, "real(kind=8)", dcc2xz)
#endif
!   ==============================================================
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acofyz", 1, "real(kind=8)", acofyz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcofyz", 1, "real(kind=8)", bcofyz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccofyz", 1, "real(kind=8)", ccofyz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcofyz", 1, "real(kind=8)", dcofyz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ecofyz", 1, "real(kind=8)", ecofyz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf1yz", 1, "real(kind=8)", acf1yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf1yz", 1, "real(kind=8)", bcf1yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf1yz", 1, "real(kind=8)", ccf1yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcf1yz", 1, "real(kind=8)", dcf1yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf2yz", 1, "real(kind=8)", acf2yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf2yz", 1, "real(kind=8)", bcf2yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf2yz", 1, "real(kind=8)", ccf2yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcf2yz", 1, "real(kind=8)", dcf2yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf3yz", 1, "real(kind=8)", acf3yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf3yz", 1, "real(kind=8)", bcf3yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf4yz", 1, "real(kind=8)", acf4yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf4yz", 1, "real(kind=8)", bcf4yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf4yz", 1, "real(kind=8)", ccf4yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acf5yz", 1, "real(kind=8)", acf5yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcf5yz", 1, "real(kind=8)", bcf5yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccf5yz", 1, "real(kind=8)", ccf5yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcf5yz", 1, "real(kind=8)", dcf5yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acc1yz", 1, "real(kind=8)", acc1yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcc1yz", 1, "real(kind=8)", bcc1yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccc1yz", 1, "real(kind=8)", ccc1yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcc1yz", 1, "real(kind=8)", dcc1yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("acc2yz", 1, "real(kind=8)", acc2yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("bcc2yz", 1, "real(kind=8)", bcc2yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("ccc2yz", 1, "real(kind=8)", ccc2yz)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("dcc2yz", 1, "real(kind=8)", dcc2yz)
#endif
#ifdef OPS_WITH_OMPOFFLOADFOR
!$OMP TARGET UPDATE TO(acoffx)
!$OMP TARGET UPDATE TO(bcoffx)
!$OMP TARGET UPDATE TO(ccoffx)
!$OMP TARGET UPDATE TO(dcoffx)
!$OMP TARGET UPDATE TO(ecoffx)
!$OMP TARGET UPDATE TO(acof1x)
!$OMP TARGET UPDATE TO(bcof1x)
!$OMP TARGET UPDATE TO(ccof1x)
!$OMP TARGET UPDATE TO(dcof1x)
!$OMP TARGET UPDATE TO(acof2x)
!$OMP TARGET UPDATE TO(bcof2x)
!$OMP TARGET UPDATE TO(ccof2x)
!$OMP TARGET UPDATE TO(dcof2x)
!$OMP TARGET UPDATE TO(acof3x)
!$OMP TARGET UPDATE TO(bcof3x)
!$OMP TARGET UPDATE TO(acof4x)
!$OMP TARGET UPDATE TO(bcof4x)
!$OMP TARGET UPDATE TO(ccof4x)
!$OMP TARGET UPDATE TO(acof5x)
!$OMP TARGET UPDATE TO(bcof5x)
!$OMP TARGET UPDATE TO(ccof5x)
!$OMP TARGET UPDATE TO(dcof5x)
!$OMP TARGET UPDATE TO(ovdelx)
!$OMP TARGET UPDATE TO(acofsx)
!$OMP TARGET UPDATE TO(bcofsx)
!$OMP TARGET UPDATE TO(ccofsx)
!$OMP TARGET UPDATE TO(dcofsx)
!$OMP TARGET UPDATE TO(ecofsx)
!$OMP TARGET UPDATE TO(acfs1x)
!$OMP TARGET UPDATE TO(bcfs1x)
!$OMP TARGET UPDATE TO(ccfs1x)
!$OMP TARGET UPDATE TO(dcfs1x)
!$OMP TARGET UPDATE TO(ecfs1x)
!$OMP TARGET UPDATE TO(acfs2x)
!$OMP TARGET UPDATE TO(bcfs2x)
!$OMP TARGET UPDATE TO(ccfs2x)
!$OMP TARGET UPDATE TO(dcfs2x)
!$OMP TARGET UPDATE TO(ecfs2x)
!$OMP TARGET UPDATE TO(acfs3x)
!$OMP TARGET UPDATE TO(bcfs3x)
!$OMP TARGET UPDATE TO(acfs4x)
!$OMP TARGET UPDATE TO(bcfs4x)
!$OMP TARGET UPDATE TO(ccfs4x)
!$OMP TARGET UPDATE TO(acfs5x)
!$OMP TARGET UPDATE TO(bcfs5x)
!$OMP TARGET UPDATE TO(ccfs5x)
!$OMP TARGET UPDATE TO(dcfs5x)
!$OMP TARGET UPDATE TO(ovdlx2)
!$OMP TARGET UPDATE TO(acoffy)
!$OMP TARGET UPDATE TO(bcoffy)
!$OMP TARGET UPDATE TO(ccoffy)
!$OMP TARGET UPDATE TO(dcoffy)
!$OMP TARGET UPDATE TO(ecoffy)
!$OMP TARGET UPDATE TO(acof1y)
!$OMP TARGET UPDATE TO(bcof1y)
!$OMP TARGET UPDATE TO(ccof1y)
!$OMP TARGET UPDATE TO(dcof1y)
!$OMP TARGET UPDATE TO(acof2y)
!$OMP TARGET UPDATE TO(bcof2y)
!$OMP TARGET UPDATE TO(ccof2y)
!$OMP TARGET UPDATE TO(dcof2y)
!$OMP TARGET UPDATE TO(acof3y)
!$OMP TARGET UPDATE TO(bcof3y)
!$OMP TARGET UPDATE TO(acof4y)
!$OMP TARGET UPDATE TO(bcof4y)
!$OMP TARGET UPDATE TO(ccof4y)
!$OMP TARGET UPDATE TO(acof5y)
!$OMP TARGET UPDATE TO(bcof5y)
!$OMP TARGET UPDATE TO(ccof5y)
!$OMP TARGET UPDATE TO(dcof5y)
!$OMP TARGET UPDATE TO(ovdely)
!$OMP TARGET UPDATE TO(acofsy)
!$OMP TARGET UPDATE TO(bcofsy)
!$OMP TARGET UPDATE TO(ccofsy)
!$OMP TARGET UPDATE TO(dcofsy)
!$OMP TARGET UPDATE TO(ecofsy)
!$OMP TARGET UPDATE TO(acfs1y)
!$OMP TARGET UPDATE TO(bcfs1y)
!$OMP TARGET UPDATE TO(ccfs1y)
!$OMP TARGET UPDATE TO(dcfs1y)
!$OMP TARGET UPDATE TO(ecfs1y)
!$OMP TARGET UPDATE TO(acfs2y)
!$OMP TARGET UPDATE TO(bcfs2y)
!$OMP TARGET UPDATE TO(ccfs2y)
!$OMP TARGET UPDATE TO(dcfs2y)
!$OMP TARGET UPDATE TO(ecfs2y)
!$OMP TARGET UPDATE TO(acfs3y)
!$OMP TARGET UPDATE TO(bcfs3y)
!$OMP TARGET UPDATE TO(acfs4y)
!$OMP TARGET UPDATE TO(bcfs4y)
!$OMP TARGET UPDATE TO(ccfs4y)
!$OMP TARGET UPDATE TO(acfs5y)
!$OMP TARGET UPDATE TO(bcfs5y)
!$OMP TARGET UPDATE TO(ccfs5y)
!$OMP TARGET UPDATE TO(dcfs5y)
!$OMP TARGET UPDATE TO(ovdly2)
!$OMP TARGET UPDATE TO(acoffz)
!$OMP TARGET UPDATE TO(bcoffz)
!$OMP TARGET UPDATE TO(ccoffz)
!$OMP TARGET UPDATE TO(dcoffz)
!$OMP TARGET UPDATE TO(ecoffz)
!$OMP TARGET UPDATE TO(acof1z)
!$OMP TARGET UPDATE TO(bcof1z)
!$OMP TARGET UPDATE TO(ccof1z)
!$OMP TARGET UPDATE TO(dcof1z)
!$OMP TARGET UPDATE TO(acof2z)
!$OMP TARGET UPDATE TO(bcof2z)
!$OMP TARGET UPDATE TO(ccof2z)
!$OMP TARGET UPDATE TO(dcof2z)
!$OMP TARGET UPDATE TO(acof3z)
!$OMP TARGET UPDATE TO(bcof3z)
!$OMP TARGET UPDATE TO(acof4z)
!$OMP TARGET UPDATE TO(bcof4z)
!$OMP TARGET UPDATE TO(ccof4z)
!$OMP TARGET UPDATE TO(acof5z)
!$OMP TARGET UPDATE TO(bcof5z)
!$OMP TARGET UPDATE TO(ccof5z)
!$OMP TARGET UPDATE TO(dcof5z)
!$OMP TARGET UPDATE TO(ovdelz)
!$OMP TARGET UPDATE TO(acofsz)
!$OMP TARGET UPDATE TO(bcofsz)
!$OMP TARGET UPDATE TO(ccofsz)
!$OMP TARGET UPDATE TO(dcofsz)
!$OMP TARGET UPDATE TO(ecofsz)
!$OMP TARGET UPDATE TO(acfs1z)
!$OMP TARGET UPDATE TO(bcfs1z)
!$OMP TARGET UPDATE TO(ccfs1z)
!$OMP TARGET UPDATE TO(dcfs1z)
!$OMP TARGET UPDATE TO(ecfs1z)
!$OMP TARGET UPDATE TO(acfs2z)
!$OMP TARGET UPDATE TO(bcfs2z)
!$OMP TARGET UPDATE TO(ccfs2z)
!$OMP TARGET UPDATE TO(dcfs2z)
!$OMP TARGET UPDATE TO(ecfs2z)
!$OMP TARGET UPDATE TO(acfs3z)
!$OMP TARGET UPDATE TO(bcfs3z)
!$OMP TARGET UPDATE TO(acfs4z)
!$OMP TARGET UPDATE TO(bcfs4z)
!$OMP TARGET UPDATE TO(ccfs4z)
!$OMP TARGET UPDATE TO(acfs5z)
!$OMP TARGET UPDATE TO(bcfs5z)
!$OMP TARGET UPDATE TO(ccfs5z)
!$OMP TARGET UPDATE TO(dcfs5z)
!$OMP TARGET UPDATE TO(ovdlz2)
!$OMP TARGET UPDATE TO(acofx1)
!$OMP TARGET UPDATE TO(bcofx1)
!$OMP TARGET UPDATE TO(acofy1)
!$OMP TARGET UPDATE TO(bcofy1)
!$OMP TARGET UPDATE TO(acofz1)
!$OMP TARGET UPDATE TO(bcofz1)
!$OMP TARGET UPDATE TO(acofxy)
!$OMP TARGET UPDATE TO(bcofxy)
!$OMP TARGET UPDATE TO(ccofxy)
!$OMP TARGET UPDATE TO(dcofxy)
!$OMP TARGET UPDATE TO(ecofxy)
!$OMP TARGET UPDATE TO(acf1xy)
!$OMP TARGET UPDATE TO(bcf1xy)
!$OMP TARGET UPDATE TO(ccf1xy)
!$OMP TARGET UPDATE TO(dcf1xy)
!$OMP TARGET UPDATE TO(acf2xy)
!$OMP TARGET UPDATE TO(bcf2xy)
!$OMP TARGET UPDATE TO(ccf2xy)
!$OMP TARGET UPDATE TO(dcf2xy)
!$OMP TARGET UPDATE TO(acf3xy)
!$OMP TARGET UPDATE TO(bcf3xy)
!$OMP TARGET UPDATE TO(acf4xy)
!$OMP TARGET UPDATE TO(bcf4xy)
!$OMP TARGET UPDATE TO(ccf4xy)
!$OMP TARGET UPDATE TO(acf5xy)
!$OMP TARGET UPDATE TO(bcf5xy)
!$OMP TARGET UPDATE TO(ccf5xy)
!$OMP TARGET UPDATE TO(dcf5xy)
!$OMP TARGET UPDATE TO(acc1xy)
!$OMP TARGET UPDATE TO(bcc1xy)
!$OMP TARGET UPDATE TO(ccc1xy)
!$OMP TARGET UPDATE TO(dcc1xy)
!$OMP TARGET UPDATE TO(acc2xy)
!$OMP TARGET UPDATE TO(bcc2xy)
!$OMP TARGET UPDATE TO(ccc2xy)
!$OMP TARGET UPDATE TO(dcc2xy)
!$OMP TARGET UPDATE TO(acofxz)
!$OMP TARGET UPDATE TO(bcofxz)
!$OMP TARGET UPDATE TO(ccofxz)
!$OMP TARGET UPDATE TO(dcofxz)
!$OMP TARGET UPDATE TO(ecofxz)
!$OMP TARGET UPDATE TO(acf1xz)
!$OMP TARGET UPDATE TO(bcf1xz)
!$OMP TARGET UPDATE TO(ccf1xz)
!$OMP TARGET UPDATE TO(dcf1xz)
!$OMP TARGET UPDATE TO(acf2xz)
!$OMP TARGET UPDATE TO(bcf2xz)
!$OMP TARGET UPDATE TO(ccf2xz)
!$OMP TARGET UPDATE TO(dcf2xz)
!$OMP TARGET UPDATE TO(acf3xz)
!$OMP TARGET UPDATE TO(bcf3xz)
!$OMP TARGET UPDATE TO(acf4xz)
!$OMP TARGET UPDATE TO(bcf4xz)
!$OMP TARGET UPDATE TO(ccf4xz)
!$OMP TARGET UPDATE TO(acf5xz)
!$OMP TARGET UPDATE TO(bcf5xz)
!$OMP TARGET UPDATE TO(ccf5xz)
!$OMP TARGET UPDATE TO(dcf5xz)
!$OMP TARGET UPDATE TO(acc1xz)
!$OMP TARGET UPDATE TO(bcc1xz)
!$OMP TARGET UPDATE TO(ccc1xz)
!$OMP TARGET UPDATE TO(dcc1xz)
!$OMP TARGET UPDATE TO(acc2xz)
!$OMP TARGET UPDATE TO(bcc2xz)
!$OMP TARGET UPDATE TO(ccc2xz)
!$OMP TARGET UPDATE TO(dcc2xz)
!$OMP TARGET UPDATE TO(acofyz)
!$OMP TARGET UPDATE TO(bcofyz)
!$OMP TARGET UPDATE TO(ccofyz)
!$OMP TARGET UPDATE TO(dcofyz)
!$OMP TARGET UPDATE TO(ecofyz)
!$OMP TARGET UPDATE TO(acf1yz)
!$OMP TARGET UPDATE TO(bcf1yz)
!$OMP TARGET UPDATE TO(ccf1yz)
!$OMP TARGET UPDATE TO(dcf1yz)
!$OMP TARGET UPDATE TO(acf2yz)
!$OMP TARGET UPDATE TO(bcf2yz)
!$OMP TARGET UPDATE TO(ccf2yz)
!$OMP TARGET UPDATE TO(dcf2yz)
!$OMP TARGET UPDATE TO(acf3yz)
!$OMP TARGET UPDATE TO(bcf3yz)
!$OMP TARGET UPDATE TO(acf4yz)
!$OMP TARGET UPDATE TO(bcf4yz)
!$OMP TARGET UPDATE TO(ccf4yz)
!$OMP TARGET UPDATE TO(acf5yz)
!$OMP TARGET UPDATE TO(bcf5yz)
!$OMP TARGET UPDATE TO(ccf5yz)
!$OMP TARGET UPDATE TO(dcf5yz)
!$OMP TARGET UPDATE TO(acc1yz)
!$OMP TARGET UPDATE TO(bcc1yz)
!$OMP TARGET UPDATE TO(ccc1yz)
!$OMP TARGET UPDATE TO(dcc1yz)
!$OMP TARGET UPDATE TO(acc2yz)
!$OMP TARGET UPDATE TO(bcc2yz)
!$OMP TARGET UPDATE TO(ccc2yz)
!$OMP TARGET UPDATE TO(dcc2yz)
#endif

!   ==============================================================

END SUBROUTINE dfinit
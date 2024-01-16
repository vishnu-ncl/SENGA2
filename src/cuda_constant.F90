SUBROUTINE cuda_const_init
    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

#ifdef OPS_WITH_CUDAFOR

!   ======================================

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

!   ======================================

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

!   ======================================

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

!   ======================================

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

!   ======================================

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

!   ======================================

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

!   ======================================

    acofx1_opsconstant = acofx1
    bcofx1_opsconstant = bcofx1
    acofy1_opsconstant = acofy1
    bcofy1_opsconstant = bcofy1
    acofz1_opsconstant = acofz1
    bcofz1_opsconstant = bcofz1

!   ======================================

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

!   ======================================

    acofxz_opsconstant = acofxz
    bcofxz_opsconstant = bcofxz
    ccofxz_opsconstant = ccofxz
    dcofxz_opsconstant = dcofxz
    ecofxz_opsconstant = ecofxz

    acf1xz_opsconstant = acf1xz
    bcf1xz_opsconstant = bcf1xz
    ccf1xz_opsconstant = ccf1xz
    dcf1xz_opsconstant = dcf1xz

    acf2xz_opsconstant = acf2xz
    bcf2xz_opsconstant = bcf2xz
    ccf2xz_opsconstant = ccf2xz
    dcf2xz_opsconstant = dcf2xz

    acf3xz_opsconstant = acf3xz
    bcf3xz_opsconstant = bcf3xz

    acf4xz_opsconstant = acf4xz
    bcf4xz_opsconstant = bcf4xz
    ccf4xz_opsconstant = ccf4xz

    acf5xz_opsconstant = acf5xz
    bcf5xz_opsconstant = bcf5xz
    ccf5xz_opsconstant = ccf5xz
    dcf5xz_opsconstant = dcf5xz

    acc1xz_opsconstant = acc1xz
    bcc1xz_opsconstant = bcc1xz
    ccc1xz_opsconstant = ccc1xz
    dcc1xz_opsconstant = dcc1xz

    acc2xz_opsconstant = acc2xz
    bcc2xz_opsconstant = bcc2xz
    ccc2xz_opsconstant = ccc2xz
    dcc2xz_opsconstant = dcc2xz

!   ======================================

    acofyz_opsconstant = acofyz
    bcofyz_opsconstant = bcofyz
    ccofyz_opsconstant = ccofyz
    dcofyz_opsconstant = dcofyz
    ecofyz_opsconstant = ecofyz

    acf1yz_opsconstant = acf1yz
    bcf1yz_opsconstant = bcf1yz
    ccf1yz_opsconstant = ccf1yz
    dcf1yz_opsconstant = dcf1yz

    acf2yz_opsconstant = acf2yz
    bcf2yz_opsconstant = bcf2yz
    ccf2yz_opsconstant = ccf2yz
    dcf2yz_opsconstant = dcf2yz

    acf3yz_opsconstant = acf3yz
    bcf3yz_opsconstant = bcf3yz

    acf4yz_opsconstant = acf4yz
    bcf4yz_opsconstant = bcf4yz
    ccf4yz_opsconstant = ccf4yz

    acf5yz_opsconstant = acf5yz
    bcf5yz_opsconstant = bcf5yz
    ccf5yz_opsconstant = ccf5yz
    dcf5yz_opsconstant = dcf5yz

    acc1yz_opsconstant = acc1yz
    bcc1yz_opsconstant = bcc1yz
    ccc1yz_opsconstant = ccc1yz
    dcc1yz_opsconstant = dcc1yz

    acc2yz_opsconstant = acc2yz
    bcc2yz_opsconstant = bcc2yz
    ccc2yz_opsconstant = ccc2yz
    dcc2yz_opsconstant = dcc2yz

!   ======================================

    foursb_opsconstant = foursb
    trfrth_opsconstant = trfrth

!   ======================================

    rlamda_opsconstant = rlamda
    alamda_opsconstant = alamda

!   ======================================

#endif

END SUBROUTINE cuda_const_init

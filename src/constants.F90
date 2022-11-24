MODULE OPS_CONSTANTS
use data_types
#ifdef OPS_WITH_CUDAFOR
    use cudafor

    real(kind=dp), constant :: acoffx_opsconstant, bcoffx_opsconstant, ccoffx_opsconstant, dcoffx_opsconstant, ecoffx_opsconstant
    real(kind=dp), constant :: acof1x_opsconstant, bcof1x_opsconstant, ccof1x_opsconstant, dcof1x_opsconstant
    real(kind=dp), constant :: acof2x_opsconstant, bcof2x_opsconstant, ccof2x_opsconstant, dcof2x_opsconstant
    real(kind=dp), constant :: acof3x_opsconstant, bcof3x_opsconstant
    real(kind=dp), constant :: acof4x_opsconstant, bcof4x_opsconstant, ccof4x_opsconstant
    real(kind=dp), constant :: acof5x_opsconstant, bcof5x_opsconstant, ccof5x_opsconstant, dcof5x_opsconstant

    real(kind=dp), constant :: ovdelx_opsconstant
    
    real(kind=dp) :: acoffx, bcoffx, ccoffx, dcoffx, ecoffx
    real(kind=dp) :: acof1x, bcof1x, ccof1x, dcof1x
    real(kind=dp) :: acof2x, bcof2x, ccof2x, dcof2x
    real(kind=dp) :: acof3x, bcof3x
    real(kind=dp) :: acof4x, bcof4x, ccof4x
    real(kind=dp) :: acof5x, bcof5x, ccof5x, dcof5x

    real(kind=dp) :: ovdelx

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=dp), constant :: acofsx_opsconstant,bcofsx_opsconstant,ccofsx_opsconstant,dcofsx_opsconstant,ecofsx_opsconstant
    real(kind=dp), constant :: acfs1x_opsconstant,bcfs1x_opsconstant,ccfs1x_opsconstant,dcfs1x_opsconstant,ecfs1x_opsconstant
    real(kind=dp), constant :: acfs2x_opsconstant,bcfs2x_opsconstant,ccfs2x_opsconstant,dcfs2x_opsconstant,ecfs2x_opsconstant
    real(kind=dp), constant :: acfs3x_opsconstant,bcfs3x_opsconstant
    real(kind=dp), constant :: acfs4x_opsconstant,bcfs4x_opsconstant,ccfs4x_opsconstant
    real(kind=dp), constant :: acfs5x_opsconstant,bcfs5x_opsconstant,ccfs5x_opsconstant,dcfs5x_opsconstant

    real(kind=dp), constant :: ovdlx2_opsconstant

    real(kind=dp) :: acofsx,bcofsx,ccofsx,dcofsx,ecofsx
    real(kind=dp) :: acfs1x,bcfs1x,ccfs1x,dcfs1x,ecfs1x
    real(kind=dp) :: acfs2x,bcfs2x,ccfs2x,dcfs2x,ecfs2x
    real(kind=dp) :: acfs3x,bcfs3x
    real(kind=dp) :: acfs4x,bcfs4x,ccfs4x
    real(kind=dp) :: acfs5x,bcfs5x,ccfs5x,dcfs5x

    real(kind=dp) :: ovdlx2

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=dp), constant :: acoffy_opsconstant, bcoffy_opsconstant, ccoffy_opsconstant, dcoffy_opsconstant, ecoffy_opsconstant
    real(kind=dp), constant :: acof1y_opsconstant, bcof1y_opsconstant, ccof1y_opsconstant, dcof1y_opsconstant
    real(kind=dp), constant :: acof2y_opsconstant, bcof2y_opsconstant, ccof2y_opsconstant, dcof2y_opsconstant
    real(kind=dp), constant :: acof3y_opsconstant, bcof3y_opsconstant
    real(kind=dp), constant :: acof4y_opsconstant, bcof4y_opsconstant, ccof4y_opsconstant
    real(kind=dp), constant :: acof5y_opsconstant, bcof5y_opsconstant, ccof5y_opsconstant, dcof5y_opsconstant

    real(kind=dp), constant :: ovdely_opsconstant

    real(kind=dp) :: acoffy, bcoffy, ccoffy, dcoffy, ecoffy
    real(kind=dp) :: acof1y, bcof1y, ccof1y, dcof1y
    real(kind=dp) :: acof2y, bcof2y, ccof2y, dcof2y
    real(kind=dp) :: acof3y, bcof3y
    real(kind=dp) :: acof4y, bcof4y, ccof4y
    real(kind=dp) :: acof5y, bcof5y, ccof5y, dcof5y

    real(kind=dp) :: ovdely

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=dp), constant :: acofsy_opsconstant,bcofsy_opsconstant,ccofsy_opsconstant,dcofsy_opsconstant,ecofsy_opsconstant
    real(kind=dp), constant :: acfs1y_opsconstant,bcfs1y_opsconstant,ccfs1y_opsconstant,dcfs1y_opsconstant,ecfs1y_opsconstant
    real(kind=dp), constant :: acfs2y_opsconstant,bcfs2y_opsconstant,ccfs2y_opsconstant,dcfs2y_opsconstant,ecfs2y_opsconstant
    real(kind=dp), constant :: acfs3y_opsconstant,bcfs3y_opsconstant
    real(kind=dp), constant :: acfs4y_opsconstant,bcfs4y_opsconstant,ccfs4y_opsconstant
    real(kind=dp), constant :: acfs5y_opsconstant,bcfs5y_opsconstant,ccfs5y_opsconstant,dcfs5y_opsconstant

    real(kind=dp), constant :: ovdly2_opsconstant

    real(kind=dp) :: acofsy,bcofsy,ccofsy,dcofsy,ecofsy
    real(kind=dp) :: acfs1y,bcfs1y,ccfs1y,dcfs1y,ecfs1y
    real(kind=dp) :: acfs2y,bcfs2y,ccfs2y,dcfs2y,ecfs2y
    real(kind=dp) :: acfs3y,bcfs3y
    real(kind=dp) :: acfs4y,bcfs4y,ccfs4y
    real(kind=dp) :: acfs5y,bcfs5y,ccfs5y,dcfs5y

    real(kind=dp) :: ovdly2

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=dp), constant :: acoffz_opsconstant, bcoffz_opsconstant, ccoffz_opsconstant, dcoffz_opsconstant, ecoffz_opsconstant
    real(kind=dp), constant :: acof1z_opsconstant, bcof1z_opsconstant, ccof1z_opsconstant, dcof1z_opsconstant
    real(kind=dp), constant :: acof2z_opsconstant, bcof2z_opsconstant, ccof2z_opsconstant, dcof2z_opsconstant
    real(kind=dp), constant :: acof3z_opsconstant, bcof3z_opsconstant
    real(kind=dp), constant :: acof4z_opsconstant, bcof4z_opsconstant, ccof4z_opsconstant
    real(kind=dp), constant :: acof5z_opsconstant, bcof5z_opsconstant, ccof5z_opsconstant, dcof5z_opsconstant

    real(kind=dp), constant :: ovdelz_opsconstant

    real(kind=dp) :: acoffz, bcoffz, ccoffz, dcoffz, ecoffz
    real(kind=dp) :: acof1z, bcof1z, ccof1z, dcof1z
    real(kind=dp) :: acof2z, bcof2z, ccof2z, dcof2z
    real(kind=dp) :: acof3z, bcof3z
    real(kind=dp) :: acof4z, bcof4z, ccof4z
    real(kind=dp) :: acof5z, bcof5z, ccof5z, dcof5z

    real(kind=dp) :: ovdelz

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=dp), constant :: acofsz_opsconstant,bcofsz_opsconstant,ccofsz_opsconstant,dcofsz_opsconstant,ecofsz_opsconstant
    real(kind=dp), constant :: acfs1z_opsconstant,bcfs1z_opsconstant,ccfs1z_opsconstant,dcfs1z_opsconstant,ecfs1z_opsconstant
    real(kind=dp), constant :: acfs2z_opsconstant,bcfs2z_opsconstant,ccfs2z_opsconstant,dcfs2z_opsconstant,ecfs2z_opsconstant
    real(kind=dp), constant :: acfs3z_opsconstant,bcfs3z_opsconstant
    real(kind=dp), constant :: acfs4z_opsconstant,bcfs4z_opsconstant,ccfs4z_opsconstant
    real(kind=dp), constant :: acfs5z_opsconstant,bcfs5z_opsconstant,ccfs5z_opsconstant,dcfs5z_opsconstant

    real(kind=dp), constant :: ovdlz2_opsconstant

    real(kind=dp) :: acofsz,bcofsz,ccofsz,dcofsz,ecofsz
    real(kind=dp) :: acfs1z,bcfs1z,ccfs1z,dcfs1z,ecfs1z
    real(kind=dp) :: acfs2z,bcfs2z,ccfs2z,dcfs2z,ecfs2z
    real(kind=dp) :: acfs3z,bcfs3z
    real(kind=dp) :: acfs4z,bcfs4z,ccfs4z
    real(kind=dp) :: acfs5z,bcfs5z,ccfs5z,dcfs5z

    real(kind=dp) :: ovdlz2

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=dp), constant :: acofx1_opsconstant,bcofx1_opsconstant
    real(kind=dp), constant :: acofy1_opsconstant,bcofy1_opsconstant
    real(kind=dp), constant :: acofz1_opsconstant,bcofz1_opsconstant

    real(kind=dp), constant :: acofxy_opsconstant,bcofxy_opsconstant,ccofxy_opsconstant,dcofxy_opsconstant,ecofxy_opsconstant

    real(kind=dp), constant :: acf1xy_opsconstant,bcf1xy_opsconstant,ccf1xy_opsconstant,dcf1xy_opsconstant
    real(kind=dp), constant :: acf2xy_opsconstant,bcf2xy_opsconstant,ccf2xy_opsconstant,dcf2xy_opsconstant
    real(kind=dp), constant :: acf3xy_opsconstant,bcf3xy_opsconstant
    real(kind=dp), constant :: acf4xy_opsconstant,bcf4xy_opsconstant,ccf4xy_opsconstant
    real(kind=dp), constant :: acf5xy_opsconstant,bcf5xy_opsconstant,ccf5xy_opsconstant,dcf5xy_opsconstant

    real(kind=dp), constant :: acc1xy_opsconstant,bcc1xy_opsconstant,ccc1xy_opsconstant,dcc1xy_opsconstant
    real(kind=dp), constant :: acc2xy_opsconstant,bcc2xy_opsconstant,ccc2xy_opsconstant,dcc2xy_opsconstant


    real(kind=dp) :: acofx1,bcofx1,acofy1,bcofy1,acofz1,bcofz1

    real(kind=dp) :: acofxy,bcofxy,ccofxy,dcofxy,ecofxy

    real(kind=dp) :: acf1xy,bcf1xy,ccf1xy,dcf1xy
    real(kind=dp) :: acf2xy,bcf2xy,ccf2xy,dcf2xy
    real(kind=dp) :: acf3xy,bcf3xy
    real(kind=dp) :: acf4xy,bcf4xy,ccf4xy
    real(kind=dp) :: acf5xy,bcf5xy,ccf5xy,dcf5xy

    real(kind=dp) :: acc1xy,bcc1xy,ccc1xy,dcc1xy
    real(kind=dp) :: acc2xy,bcc2xy,ccc2xy,dcc2xy            
    
!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=dp), constant :: foursb_opsconstant,trfrth_opsconstant
    
    real(kind=dp) :: foursb,trfrth
!---------------------------------------------------------------------------------------------------------------------------------
#else
!---------------------------------------------------------------------------------------------------------------------------------
    real(kind=dp) :: acoffx, bcoffx, ccoffx, dcoffx, ecoffx
    real(kind=dp) :: acof1x, bcof1x, ccof1x, dcof1x
    real(kind=dp) :: acof2x, bcof2x, ccof2x, dcof2x
    real(kind=dp) :: acof3x, bcof3x
    real(kind=dp) :: acof4x, bcof4x, ccof4x
    real(kind=dp) :: acof5x, bcof5x, ccof5x, dcof5x

    real(kind=dp) :: ovdelx

!-------------------------------------------------------------------

    real(kind=dp) :: acofsx,bcofsx,ccofsx,dcofsx,ecofsx
    real(kind=dp) :: acfs1x,bcfs1x,ccfs1x,dcfs1x,ecfs1x
    real(kind=dp) :: acfs2x,bcfs2x,ccfs2x,dcfs2x,ecfs2x
    real(kind=dp) :: acfs3x,bcfs3x
    real(kind=dp) :: acfs4x,bcfs4x,ccfs4x
    real(kind=dp) :: acfs5x,bcfs5x,ccfs5x,dcfs5x

    real(kind=dp) :: ovdlx2

!-------------------------------------------------------------------

    real(kind=dp) :: acoffy, bcoffy, ccoffy, dcoffy, ecoffy
    real(kind=dp) :: acof1y, bcof1y, ccof1y, dcof1y
    real(kind=dp) :: acof2y, bcof2y, ccof2y, dcof2y
    real(kind=dp) :: acof3y, bcof3y
    real(kind=dp) :: acof4y, bcof4y, ccof4y
    real(kind=dp) :: acof5y, bcof5y, ccof5y, dcof5y

    real(kind=dp) :: ovdely

!-------------------------------------------------------------------

    real(kind=dp) :: acofsy,bcofsy,ccofsy,dcofsy,ecofsy
    real(kind=dp) :: acfs1y,bcfs1y,ccfs1y,dcfs1y,ecfs1y
    real(kind=dp) :: acfs2y,bcfs2y,ccfs2y,dcfs2y,ecfs2y
    real(kind=dp) :: acfs3y,bcfs3y
    real(kind=dp) :: acfs4y,bcfs4y,ccfs4y
    real(kind=dp) :: acfs5y,bcfs5y,ccfs5y,dcfs5y

    real(kind=dp) :: ovdly2

!-------------------------------------------------------------------

    real(kind=dp) :: acoffz, bcoffz, ccoffz, dcoffz, ecoffz
    real(kind=dp) :: acof1z, bcof1z, ccof1z, dcof1z
    real(kind=dp) :: acof2z, bcof2z, ccof2z, dcof2z
    real(kind=dp) :: acof3z, bcof3z
    real(kind=dp) :: acof4z, bcof4z, ccof4z
    real(kind=dp) :: acof5z, bcof5z, ccof5z, dcof5z

    real(kind=dp) :: ovdelz

!-------------------------------------------------------------------

    real(kind=dp) :: acofsz,bcofsz,ccofsz,dcofsz,ecofsz
    real(kind=dp) :: acfs1z,bcfs1z,ccfs1z,dcfs1z,ecfs1z
    real(kind=dp) :: acfs2z,bcfs2z,ccfs2z,dcfs2z,ecfs2z
    real(kind=dp) :: acfs3z,bcfs3z
    real(kind=dp) :: acfs4z,bcfs4z,ccfs4z
    real(kind=dp) :: acfs5z,bcfs5z,ccfs5z,dcfs5z

    real(kind=dp) :: ovdlz2

!-------------------------------------------------------------------

    real(kind=dp) :: acofx1,bcofx1,acofy1,bcofy1,acofz1,bcofz1
    
    real(kind=dp) :: acofxy,bcofxy,ccofxy,dcofxy,ecofxy

    real(kind=dp) :: acf1xy,bcf1xy,ccf1xy,dcf1xy
    real(kind=dp) :: acf2xy,bcf2xy,ccf2xy,dcf2xy
    real(kind=dp) :: acf3xy,bcf3xy
    real(kind=dp) :: acf4xy,bcf4xy,ccf4xy
    real(kind=dp) :: acf5xy,bcf5xy,ccf5xy,dcf5xy

    real(kind=dp) :: acc1xy,bcc1xy,ccc1xy,dcc1xy
    real(kind=dp) :: acc2xy,bcc2xy,ccc2xy,dcc2xy

!-------------------------------------------------------------------

    real(kind=dp) :: foursb,trfrth

#endif

END MODULE OPS_CONSTANTS

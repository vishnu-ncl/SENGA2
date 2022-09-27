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

    real(kind=dp), constant :: acofsx_opsconstant,bcofsx_opsconstant,ccofsx_opsconstant,dcofsx_opsconstant,ecofsx_opsconstant
    real(kind=dp), constant :: acfs1x_opsconstant,bcfs1x_opsconstant,ccfs1x_opsconstant,dcfs1x_opsconstant,ecfs1x_opsconstant
    real(kind=dp), constant :: acfs2x_opsconstant,bcfs2x_opsconstant,ccfs2x_opsconstant,dcfs2x_opsconstant,ecfs2x_opsconstant
    real(kind=dp), constant :: acfs3x_opsconstant,bcfs3x_opsconstant
    real(kind=dp), constant :: acfs4x_opsconstant,bcfs4x_opsconstant,ccfs4x_opsconstant
    real(kind=dp), constant :: acfs5x_opsconstant,bcfs5x_opsconstant,ccfs5x_opsconstant,dcfs5x_opsconstant

    real(kind=dp), constant :: ovdlx2_opsconstant

    real(kind=dp) :: acoffx, bcoffx, ccoffx, dcoffx, ecoffx
    real(kind=dp) :: acof1x, bcof1x, ccof1x, dcof1x
    real(kind=dp) :: acof2x, bcof2x, ccof2x, dcof2x
    real(kind=dp) :: acof3x, bcof3x
    real(kind=dp) :: acof4x, bcof4x, ccof4x
    real(kind=dp) :: acof5x, bcof5x, ccof5x, dcof5x

    real(kind=dp) :: ovdelx

    real(kind=dp) :: acofsx,bcofsx,ccofsx,dcofsx,ecofsx
    real(kind=dp) :: acfs1x,bcfs1x,ccfs1x,dcfs1x,ecfs1x
    real(kind=dp) :: acfs2x,bcfs2x,ccfs2x,dcfs2x,ecfs2x
    real(kind=dp) :: acfs3x,bcfs3x
    real(kind=dp) :: acfs4x,bcfs4x,ccfs4x
    real(kind=dp) :: acfs5x,bcfs5x,ccfs5x,dcfs5x

    real(kind=dp) :: ovdlx2

#else
    real(kind=dp) :: acoffx, bcoffx, ccoffx, dcoffx, ecoffx
    real(kind=dp) :: acof1x, bcof1x, ccof1x, dcof1x
    real(kind=dp) :: acof2x, bcof2x, ccof2x, dcof2x
    real(kind=dp) :: acof3x, bcof3x
    real(kind=dp) :: acof4x, bcof4x, ccof4x
    real(kind=dp) :: acof5x, bcof5x, ccof5x, dcof5x

    real(kind=dp) :: ovdelx

    real(kind=dp) :: acofsx,bcofsx,ccofsx,dcofsx,ecofsx
    real(kind=dp) :: acfs1x,bcfs1x,ccfs1x,dcfs1x,ecfs1x
    real(kind=dp) :: acfs2x,bcfs2x,ccfs2x,dcfs2x,ecfs2x
    real(kind=dp) :: acfs3x,bcfs3x
    real(kind=dp) :: acfs4x,bcfs4x,ccfs4x
    real(kind=dp) :: acfs5x,bcfs5x,ccfs5x,dcfs5x

    real(kind=dp) :: ovdlx2

#endif

END MODULE OPS_CONSTANTS

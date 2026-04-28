MODULE OPS_CONSTANTS
#ifdef OPS_WITH_CUDAFOR
    use cudafor

    real(kind=8), constant :: acoffx_opsconstant, bcoffx_opsconstant, ccoffx_opsconstant, dcoffx_opsconstant, ecoffx_opsconstant
    real(kind=8), constant :: acof1x_opsconstant, bcof1x_opsconstant, ccof1x_opsconstant, dcof1x_opsconstant
    real(kind=8), constant :: acof2x_opsconstant, bcof2x_opsconstant, ccof2x_opsconstant, dcof2x_opsconstant
    real(kind=8), constant :: acof3x_opsconstant, bcof3x_opsconstant
    real(kind=8), constant :: acof4x_opsconstant, bcof4x_opsconstant, ccof4x_opsconstant
    real(kind=8), constant :: acof5x_opsconstant, bcof5x_opsconstant, ccof5x_opsconstant, dcof5x_opsconstant

    real(kind=8), constant :: ovdelx_opsconstant

    real(kind=8) :: acoffx, bcoffx, ccoffx, dcoffx, ecoffx
    real(kind=8) :: acof1x, bcof1x, ccof1x, dcof1x
    real(kind=8) :: acof2x, bcof2x, ccof2x, dcof2x
    real(kind=8) :: acof3x, bcof3x
    real(kind=8) :: acof4x, bcof4x, ccof4x
    real(kind=8) :: acof5x, bcof5x, ccof5x, dcof5x

    real(kind=8) :: ovdelx

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=8), constant :: acofsx_opsconstant,bcofsx_opsconstant,ccofsx_opsconstant,dcofsx_opsconstant,ecofsx_opsconstant
    real(kind=8), constant :: acfs1x_opsconstant,bcfs1x_opsconstant,ccfs1x_opsconstant,dcfs1x_opsconstant,ecfs1x_opsconstant
    real(kind=8), constant :: acfs2x_opsconstant,bcfs2x_opsconstant,ccfs2x_opsconstant,dcfs2x_opsconstant,ecfs2x_opsconstant
    real(kind=8), constant :: acfs3x_opsconstant,bcfs3x_opsconstant
    real(kind=8), constant :: acfs4x_opsconstant,bcfs4x_opsconstant,ccfs4x_opsconstant
    real(kind=8), constant :: acfs5x_opsconstant,bcfs5x_opsconstant,ccfs5x_opsconstant,dcfs5x_opsconstant

    real(kind=8), constant :: ovdlx2_opsconstant

    real(kind=8) :: acofsx,bcofsx,ccofsx,dcofsx,ecofsx
    real(kind=8) :: acfs1x,bcfs1x,ccfs1x,dcfs1x,ecfs1x
    real(kind=8) :: acfs2x,bcfs2x,ccfs2x,dcfs2x,ecfs2x
    real(kind=8) :: acfs3x,bcfs3x
    real(kind=8) :: acfs4x,bcfs4x,ccfs4x
    real(kind=8) :: acfs5x,bcfs5x,ccfs5x,dcfs5x

    real(kind=8) :: ovdlx2

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=8), constant :: acoffy_opsconstant, bcoffy_opsconstant, ccoffy_opsconstant, dcoffy_opsconstant, ecoffy_opsconstant
    real(kind=8), constant :: acof1y_opsconstant, bcof1y_opsconstant, ccof1y_opsconstant, dcof1y_opsconstant
    real(kind=8), constant :: acof2y_opsconstant, bcof2y_opsconstant, ccof2y_opsconstant, dcof2y_opsconstant
    real(kind=8), constant :: acof3y_opsconstant, bcof3y_opsconstant
    real(kind=8), constant :: acof4y_opsconstant, bcof4y_opsconstant, ccof4y_opsconstant
    real(kind=8), constant :: acof5y_opsconstant, bcof5y_opsconstant, ccof5y_opsconstant, dcof5y_opsconstant

    real(kind=8), constant :: ovdely_opsconstant

    real(kind=8) :: acoffy, bcoffy, ccoffy, dcoffy, ecoffy
    real(kind=8) :: acof1y, bcof1y, ccof1y, dcof1y
    real(kind=8) :: acof2y, bcof2y, ccof2y, dcof2y
    real(kind=8) :: acof3y, bcof3y
    real(kind=8) :: acof4y, bcof4y, ccof4y
    real(kind=8) :: acof5y, bcof5y, ccof5y, dcof5y

    real(kind=8) :: ovdely

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=8), constant :: acofsy_opsconstant,bcofsy_opsconstant,ccofsy_opsconstant,dcofsy_opsconstant,ecofsy_opsconstant
    real(kind=8), constant :: acfs1y_opsconstant,bcfs1y_opsconstant,ccfs1y_opsconstant,dcfs1y_opsconstant,ecfs1y_opsconstant
    real(kind=8), constant :: acfs2y_opsconstant,bcfs2y_opsconstant,ccfs2y_opsconstant,dcfs2y_opsconstant,ecfs2y_opsconstant
    real(kind=8), constant :: acfs3y_opsconstant,bcfs3y_opsconstant
    real(kind=8), constant :: acfs4y_opsconstant,bcfs4y_opsconstant,ccfs4y_opsconstant
    real(kind=8), constant :: acfs5y_opsconstant,bcfs5y_opsconstant,ccfs5y_opsconstant,dcfs5y_opsconstant

    real(kind=8), constant :: ovdly2_opsconstant

    real(kind=8) :: acofsy,bcofsy,ccofsy,dcofsy,ecofsy
    real(kind=8) :: acfs1y,bcfs1y,ccfs1y,dcfs1y,ecfs1y
    real(kind=8) :: acfs2y,bcfs2y,ccfs2y,dcfs2y,ecfs2y
    real(kind=8) :: acfs3y,bcfs3y
    real(kind=8) :: acfs4y,bcfs4y,ccfs4y
    real(kind=8) :: acfs5y,bcfs5y,ccfs5y,dcfs5y

    real(kind=8) :: ovdly2

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=8), constant :: acoffz_opsconstant, bcoffz_opsconstant, ccoffz_opsconstant, dcoffz_opsconstant, ecoffz_opsconstant
    real(kind=8), constant :: acof1z_opsconstant, bcof1z_opsconstant, ccof1z_opsconstant, dcof1z_opsconstant
    real(kind=8), constant :: acof2z_opsconstant, bcof2z_opsconstant, ccof2z_opsconstant, dcof2z_opsconstant
    real(kind=8), constant :: acof3z_opsconstant, bcof3z_opsconstant
    real(kind=8), constant :: acof4z_opsconstant, bcof4z_opsconstant, ccof4z_opsconstant
    real(kind=8), constant :: acof5z_opsconstant, bcof5z_opsconstant, ccof5z_opsconstant, dcof5z_opsconstant

    real(kind=8), constant :: ovdelz_opsconstant

    real(kind=8) :: acoffz, bcoffz, ccoffz, dcoffz, ecoffz
    real(kind=8) :: acof1z, bcof1z, ccof1z, dcof1z
    real(kind=8) :: acof2z, bcof2z, ccof2z, dcof2z
    real(kind=8) :: acof3z, bcof3z
    real(kind=8) :: acof4z, bcof4z, ccof4z
    real(kind=8) :: acof5z, bcof5z, ccof5z, dcof5z

    real(kind=8) :: ovdelz

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=8), constant :: acofsz_opsconstant,bcofsz_opsconstant,ccofsz_opsconstant,dcofsz_opsconstant,ecofsz_opsconstant
    real(kind=8), constant :: acfs1z_opsconstant,bcfs1z_opsconstant,ccfs1z_opsconstant,dcfs1z_opsconstant,ecfs1z_opsconstant
    real(kind=8), constant :: acfs2z_opsconstant,bcfs2z_opsconstant,ccfs2z_opsconstant,dcfs2z_opsconstant,ecfs2z_opsconstant
    real(kind=8), constant :: acfs3z_opsconstant,bcfs3z_opsconstant
    real(kind=8), constant :: acfs4z_opsconstant,bcfs4z_opsconstant,ccfs4z_opsconstant
    real(kind=8), constant :: acfs5z_opsconstant,bcfs5z_opsconstant,ccfs5z_opsconstant,dcfs5z_opsconstant

    real(kind=8), constant :: ovdlz2_opsconstant

    real(kind=8) :: acofsz,bcofsz,ccofsz,dcofsz,ecofsz
    real(kind=8) :: acfs1z,bcfs1z,ccfs1z,dcfs1z,ecfs1z
    real(kind=8) :: acfs2z,bcfs2z,ccfs2z,dcfs2z,ecfs2z
    real(kind=8) :: acfs3z,bcfs3z
    real(kind=8) :: acfs4z,bcfs4z,ccfs4z
    real(kind=8) :: acfs5z,bcfs5z,ccfs5z,dcfs5z

    real(kind=8) :: ovdlz2

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=8), constant :: acofx1_opsconstant,bcofx1_opsconstant
    real(kind=8), constant :: acofy1_opsconstant,bcofy1_opsconstant
    real(kind=8), constant :: acofz1_opsconstant,bcofz1_opsconstant

    real(kind=8) :: acofx1,bcofx1,acofy1,bcofy1,acofz1,bcofz1

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=8), constant :: acofxy_opsconstant,bcofxy_opsconstant,ccofxy_opsconstant,dcofxy_opsconstant,ecofxy_opsconstant

    real(kind=8), constant :: acf1xy_opsconstant,bcf1xy_opsconstant,ccf1xy_opsconstant,dcf1xy_opsconstant
    real(kind=8), constant :: acf2xy_opsconstant,bcf2xy_opsconstant,ccf2xy_opsconstant,dcf2xy_opsconstant
    real(kind=8), constant :: acf3xy_opsconstant,bcf3xy_opsconstant
    real(kind=8), constant :: acf4xy_opsconstant,bcf4xy_opsconstant,ccf4xy_opsconstant
    real(kind=8), constant :: acf5xy_opsconstant,bcf5xy_opsconstant,ccf5xy_opsconstant,dcf5xy_opsconstant

    real(kind=8), constant :: acc1xy_opsconstant,bcc1xy_opsconstant,ccc1xy_opsconstant,dcc1xy_opsconstant
    real(kind=8), constant :: acc2xy_opsconstant,bcc2xy_opsconstant,ccc2xy_opsconstant,dcc2xy_opsconstant

    real(kind=8) :: acofxy,bcofxy,ccofxy,dcofxy,ecofxy

    real(kind=8) :: acf1xy,bcf1xy,ccf1xy,dcf1xy
    real(kind=8) :: acf2xy,bcf2xy,ccf2xy,dcf2xy
    real(kind=8) :: acf3xy,bcf3xy
    real(kind=8) :: acf4xy,bcf4xy,ccf4xy
    real(kind=8) :: acf5xy,bcf5xy,ccf5xy,dcf5xy

    real(kind=8) :: acc1xy,bcc1xy,ccc1xy,dcc1xy
    real(kind=8) :: acc2xy,bcc2xy,ccc2xy,dcc2xy

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=8), constant :: acofxz_opsconstant,bcofxz_opsconstant,ccofxz_opsconstant,dcofxz_opsconstant,ecofxz_opsconstant

    real(kind=8), constant :: acf1xz_opsconstant,bcf1xz_opsconstant,ccf1xz_opsconstant,dcf1xz_opsconstant
    real(kind=8), constant :: acf2xz_opsconstant,bcf2xz_opsconstant,ccf2xz_opsconstant,dcf2xz_opsconstant
    real(kind=8), constant :: acf3xz_opsconstant,bcf3xz_opsconstant
    real(kind=8), constant :: acf4xz_opsconstant,bcf4xz_opsconstant,ccf4xz_opsconstant
    real(kind=8), constant :: acf5xz_opsconstant,bcf5xz_opsconstant,ccf5xz_opsconstant,dcf5xz_opsconstant

    real(kind=8), constant :: acc1xz_opsconstant,bcc1xz_opsconstant,ccc1xz_opsconstant,dcc1xz_opsconstant
    real(kind=8), constant :: acc2xz_opsconstant,bcc2xz_opsconstant,ccc2xz_opsconstant,dcc2xz_opsconstant

    real(kind=8) :: acofxz,bcofxz,ccofxz,dcofxz,ecofxz

    real(kind=8) :: acf1xz,bcf1xz,ccf1xz,dcf1xz
    real(kind=8) :: acf2xz,bcf2xz,ccf2xz,dcf2xz
    real(kind=8) :: acf3xz,bcf3xz
    real(kind=8) :: acf4xz,bcf4xz,ccf4xz
    real(kind=8) :: acf5xz,bcf5xz,ccf5xz,dcf5xz

    real(kind=8) :: acc1xz,bcc1xz,ccc1xz,dcc1xz
    real(kind=8) :: acc2xz,bcc2xz,ccc2xz,dcc2xz

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=8), constant :: acofyz_opsconstant,bcofyz_opsconstant,ccofyz_opsconstant,dcofyz_opsconstant,ecofyz_opsconstant

    real(kind=8), constant :: acf1yz_opsconstant,bcf1yz_opsconstant,ccf1yz_opsconstant,dcf1yz_opsconstant
    real(kind=8), constant :: acf2yz_opsconstant,bcf2yz_opsconstant,ccf2yz_opsconstant,dcf2yz_opsconstant
    real(kind=8), constant :: acf3yz_opsconstant,bcf3yz_opsconstant
    real(kind=8), constant :: acf4yz_opsconstant,bcf4yz_opsconstant,ccf4yz_opsconstant
    real(kind=8), constant :: acf5yz_opsconstant,bcf5yz_opsconstant,ccf5yz_opsconstant,dcf5yz_opsconstant

    real(kind=8), constant :: acc1yz_opsconstant,bcc1yz_opsconstant,ccc1yz_opsconstant,dcc1yz_opsconstant
    real(kind=8), constant :: acc2yz_opsconstant,bcc2yz_opsconstant,ccc2yz_opsconstant,dcc2yz_opsconstant

    real(kind=8) :: acofyz,bcofyz,ccofyz,dcofyz,ecofyz

    real(kind=8) :: acf1yz,bcf1yz,ccf1yz,dcf1yz
    real(kind=8) :: acf2yz,bcf2yz,ccf2yz,dcf2yz
    real(kind=8) :: acf3yz,bcf3yz
    real(kind=8) :: acf4yz,bcf4yz,ccf4yz
    real(kind=8) :: acf5yz,bcf5yz,ccf5yz,dcf5yz

    real(kind=8) :: acc1yz,bcc1yz,ccc1yz,dcc1yz
    real(kind=8) :: acc2yz,bcc2yz,ccc2yz,dcc2yz

!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=8), constant :: foursb_opsconstant,trfrth_opsconstant

    real(kind=8) :: foursb,trfrth
!---------------------------------------------------------------------------------------------------------------------------------

    real(kind=8), constant :: rlamda_opsconstant,alamda_opsconstant

    real(kind=8) :: rlamda,alamda
!---------------------------------------------------------------------------------------------------------------------------------

    integer(kind=4), constant :: nspcmx_opsconstant,nstpmx_opsconstant
    integer(kind=4), constant :: ntinmx_opsconstant,ncofmx_opsconstant
    integer(kind=4), constant :: nctmax_opsconstant,nctmm1_opsconstant
    integer(kind=4), constant :: nspimx_opsconstant,ntbase_opsconstant,nintmx_opsconstant
    integer(kind=4), constant :: nssmax_opsconstant,nrsmax_opsconstant
    integer(kind=4), constant :: ndcfmx_opsconstant,nvcfmx_opsconstant,nccfmx_opsconstant
    integer(kind=4), constant :: ncfrmx_opsconstant
    integer(kind=4), constant :: nrkmax_opsconstant
    integer(kind=4), constant :: nbcpri_opsconstant,nbcprr_opsconstant
    integer(kind=4), constant :: ncbcsz_opsconstant

    integer(kind=4), parameter :: nspcmx=22,  nstpmx=104
    integer(kind=4), parameter :: ntinmx=2,  ncofmx=7
    integer(kind=4), parameter :: nctmax=5,  nctmm1=nctmax-1
    integer(kind=4), parameter :: nspimx=15, ntbase=4, nintmx=3
    integer(kind=4), parameter :: nssmax=10, nrsmax=10
    integer(kind=4), parameter :: ndcfmx=4,  nvcfmx=4, nccfmx = 4
    integer(kind=4), parameter :: ncfrmx=6
    integer(kind=4), parameter :: nrkmax=5
    integer(kind=4), parameter :: nbcpri=4,  nbcprr=4
    integer(kind=4), parameter :: ncbcsz=5

!---------------------------------------------------------------------------------------------------------------------------------
#else
!---------------------------------------------------------------------------------------------------------------------------------
    real(kind=8) :: acoffx, bcoffx, ccoffx, dcoffx, ecoffx
    real(kind=8) :: acof1x, bcof1x, ccof1x, dcof1x
    real(kind=8) :: acof2x, bcof2x, ccof2x, dcof2x
    real(kind=8) :: acof3x, bcof3x
    real(kind=8) :: acof4x, bcof4x, ccof4x
    real(kind=8) :: acof5x, bcof5x, ccof5x, dcof5x

    real(kind=8) :: ovdelx

!-------------------------------------------------------------------

    real(kind=8) :: acofsx,bcofsx,ccofsx,dcofsx,ecofsx
    real(kind=8) :: acfs1x,bcfs1x,ccfs1x,dcfs1x,ecfs1x
    real(kind=8) :: acfs2x,bcfs2x,ccfs2x,dcfs2x,ecfs2x
    real(kind=8) :: acfs3x,bcfs3x
    real(kind=8) :: acfs4x,bcfs4x,ccfs4x
    real(kind=8) :: acfs5x,bcfs5x,ccfs5x,dcfs5x

    real(kind=8) :: ovdlx2

!-------------------------------------------------------------------

    real(kind=8) :: acoffy, bcoffy, ccoffy, dcoffy, ecoffy
    real(kind=8) :: acof1y, bcof1y, ccof1y, dcof1y
    real(kind=8) :: acof2y, bcof2y, ccof2y, dcof2y
    real(kind=8) :: acof3y, bcof3y
    real(kind=8) :: acof4y, bcof4y, ccof4y
    real(kind=8) :: acof5y, bcof5y, ccof5y, dcof5y

    real(kind=8) :: ovdely

!-------------------------------------------------------------------

    real(kind=8) :: acofsy,bcofsy,ccofsy,dcofsy,ecofsy
    real(kind=8) :: acfs1y,bcfs1y,ccfs1y,dcfs1y,ecfs1y
    real(kind=8) :: acfs2y,bcfs2y,ccfs2y,dcfs2y,ecfs2y
    real(kind=8) :: acfs3y,bcfs3y
    real(kind=8) :: acfs4y,bcfs4y,ccfs4y
    real(kind=8) :: acfs5y,bcfs5y,ccfs5y,dcfs5y

    real(kind=8) :: ovdly2

!-------------------------------------------------------------------

    real(kind=8) :: acoffz, bcoffz, ccoffz, dcoffz, ecoffz
    real(kind=8) :: acof1z, bcof1z, ccof1z, dcof1z
    real(kind=8) :: acof2z, bcof2z, ccof2z, dcof2z
    real(kind=8) :: acof3z, bcof3z
    real(kind=8) :: acof4z, bcof4z, ccof4z
    real(kind=8) :: acof5z, bcof5z, ccof5z, dcof5z

    real(kind=8) :: ovdelz

!-------------------------------------------------------------------

    real(kind=8) :: acofsz,bcofsz,ccofsz,dcofsz,ecofsz
    real(kind=8) :: acfs1z,bcfs1z,ccfs1z,dcfs1z,ecfs1z
    real(kind=8) :: acfs2z,bcfs2z,ccfs2z,dcfs2z,ecfs2z
    real(kind=8) :: acfs3z,bcfs3z
    real(kind=8) :: acfs4z,bcfs4z,ccfs4z
    real(kind=8) :: acfs5z,bcfs5z,ccfs5z,dcfs5z

    real(kind=8) :: ovdlz2

!-------------------------------------------------------------------

    real(kind=8) :: acofx1,bcofx1,acofy1,bcofy1,acofz1,bcofz1

!-------------------------------------------------------------------

    real(kind=8) :: acofxy,bcofxy,ccofxy,dcofxy,ecofxy

    real(kind=8) :: acf1xy,bcf1xy,ccf1xy,dcf1xy
    real(kind=8) :: acf2xy,bcf2xy,ccf2xy,dcf2xy
    real(kind=8) :: acf3xy,bcf3xy
    real(kind=8) :: acf4xy,bcf4xy,ccf4xy
    real(kind=8) :: acf5xy,bcf5xy,ccf5xy,dcf5xy

    real(kind=8) :: acc1xy,bcc1xy,ccc1xy,dcc1xy
    real(kind=8) :: acc2xy,bcc2xy,ccc2xy,dcc2xy

!-------------------------------------------------------------------

    real(kind=8) :: acofxz,bcofxz,ccofxz,dcofxz,ecofxz

    real(kind=8) :: acf1xz,bcf1xz,ccf1xz,dcf1xz
    real(kind=8) :: acf2xz,bcf2xz,ccf2xz,dcf2xz
    real(kind=8) :: acf3xz,bcf3xz
    real(kind=8) :: acf4xz,bcf4xz,ccf4xz
    real(kind=8) :: acf5xz,bcf5xz,ccf5xz,dcf5xz

    real(kind=8) :: acc1xz,bcc1xz,ccc1xz,dcc1xz
    real(kind=8) :: acc2xz,bcc2xz,ccc2xz,dcc2xz

!-------------------------------------------------------------------

    real(kind=8) :: acofyz,bcofyz,ccofyz,dcofyz,ecofyz

    real(kind=8) :: acf1yz,bcf1yz,ccf1yz,dcf1yz
    real(kind=8) :: acf2yz,bcf2yz,ccf2yz,dcf2yz
    real(kind=8) :: acf3yz,bcf3yz
    real(kind=8) :: acf4yz,bcf4yz,ccf4yz
    real(kind=8) :: acf5yz,bcf5yz,ccf5yz,dcf5yz

    real(kind=8) :: acc1yz,bcc1yz,ccc1yz,dcc1yz
    real(kind=8) :: acc2yz,bcc2yz,ccc2yz,dcc2yz

!-------------------------------------------------------------------

    real(kind=8) :: foursb,trfrth

!-------------------------------------------------------------------

    real(kind=8) :: rlamda,alamda

!-------------------------------------------------------------------

    integer(kind=4), parameter :: nspcmx=22,  nstpmx=104
    integer(kind=4), parameter :: ntinmx=2,  ncofmx=7
    integer(kind=4), parameter :: nctmax=5,  nctmm1=nctmax-1
    integer(kind=4), parameter :: nspimx=15, ntbase=4, nintmx=3
    integer(kind=4), parameter :: nssmax=10, nrsmax=10
    integer(kind=4), parameter :: ndcfmx=4,  nvcfmx=4, nccfmx = 4
    integer(kind=4), parameter :: ncfrmx=6
    integer(kind=4), parameter :: nrkmax=5
    integer(kind=4), parameter :: nbcpri=4,  nbcprr=4
    integer(kind=4), parameter :: ncbcsz=5

!-------------------------------------------------------------------


#ifdef OPS_WITH_OMPOFFLOADFOR
!$OMP DECLARE TARGET(acoffx)
!$OMP DECLARE TARGET(bcoffx)
!$OMP DECLARE TARGET(ccoffx)
!$OMP DECLARE TARGET(dcoffx)
!$OMP DECLARE TARGET(ecoffx)
!$OMP DECLARE TARGET(acof1x)
!$OMP DECLARE TARGET(bcof1x)
!$OMP DECLARE TARGET(ccof1x)
!$OMP DECLARE TARGET(dcof1x)
!$OMP DECLARE TARGET(acof2x)
!$OMP DECLARE TARGET(bcof2x)
!$OMP DECLARE TARGET(ccof2x)
!$OMP DECLARE TARGET(dcof2x)
!$OMP DECLARE TARGET(acof3x)
!$OMP DECLARE TARGET(bcof3x)
!$OMP DECLARE TARGET(acof4x)
!$OMP DECLARE TARGET(bcof4x)
!$OMP DECLARE TARGET(ccof4x)
!$OMP DECLARE TARGET(acof5x)
!$OMP DECLARE TARGET(bcof5x)
!$OMP DECLARE TARGET(ccof5x)
!$OMP DECLARE TARGET(dcof5x)
!$OMP DECLARE TARGET(ovdelx)
!$OMP DECLARE TARGET(acofsx)
!$OMP DECLARE TARGET(bcofsx)
!$OMP DECLARE TARGET(ccofsx)
!$OMP DECLARE TARGET(dcofsx)
!$OMP DECLARE TARGET(ecofsx)
!$OMP DECLARE TARGET(acfs1x)
!$OMP DECLARE TARGET(bcfs1x)
!$OMP DECLARE TARGET(ccfs1x)
!$OMP DECLARE TARGET(dcfs1x)
!$OMP DECLARE TARGET(ecfs1x)
!$OMP DECLARE TARGET(acfs2x)
!$OMP DECLARE TARGET(bcfs2x)
!$OMP DECLARE TARGET(ccfs2x)
!$OMP DECLARE TARGET(dcfs2x)
!$OMP DECLARE TARGET(ecfs2x)
!$OMP DECLARE TARGET(acfs3x)
!$OMP DECLARE TARGET(bcfs3x)
!$OMP DECLARE TARGET(acfs4x)
!$OMP DECLARE TARGET(bcfs4x)
!$OMP DECLARE TARGET(ccfs4x)
!$OMP DECLARE TARGET(acfs5x)
!$OMP DECLARE TARGET(bcfs5x)
!$OMP DECLARE TARGET(ccfs5x)
!$OMP DECLARE TARGET(dcfs5x)
!$OMP DECLARE TARGET(ovdlx2)
!$OMP DECLARE TARGET(acoffy)
!$OMP DECLARE TARGET(bcoffy)
!$OMP DECLARE TARGET(ccoffy)
!$OMP DECLARE TARGET(dcoffy)
!$OMP DECLARE TARGET(ecoffy)
!$OMP DECLARE TARGET(acof1y)
!$OMP DECLARE TARGET(bcof1y)
!$OMP DECLARE TARGET(ccof1y)
!$OMP DECLARE TARGET(dcof1y)
!$OMP DECLARE TARGET(acof2y)
!$OMP DECLARE TARGET(bcof2y)
!$OMP DECLARE TARGET(ccof2y)
!$OMP DECLARE TARGET(dcof2y)
!$OMP DECLARE TARGET(acof3y)
!$OMP DECLARE TARGET(bcof3y)
!$OMP DECLARE TARGET(acof4y)
!$OMP DECLARE TARGET(bcof4y)
!$OMP DECLARE TARGET(ccof4y)
!$OMP DECLARE TARGET(acof5y)
!$OMP DECLARE TARGET(bcof5y)
!$OMP DECLARE TARGET(ccof5y)
!$OMP DECLARE TARGET(dcof5y)
!$OMP DECLARE TARGET(ovdely)
!$OMP DECLARE TARGET(acofsy)
!$OMP DECLARE TARGET(bcofsy)
!$OMP DECLARE TARGET(ccofsy)
!$OMP DECLARE TARGET(dcofsy)
!$OMP DECLARE TARGET(ecofsy)
!$OMP DECLARE TARGET(acfs1y)
!$OMP DECLARE TARGET(bcfs1y)
!$OMP DECLARE TARGET(ccfs1y)
!$OMP DECLARE TARGET(dcfs1y)
!$OMP DECLARE TARGET(ecfs1y)
!$OMP DECLARE TARGET(acfs2y)
!$OMP DECLARE TARGET(bcfs2y)
!$OMP DECLARE TARGET(ccfs2y)
!$OMP DECLARE TARGET(dcfs2y)
!$OMP DECLARE TARGET(ecfs2y)
!$OMP DECLARE TARGET(acfs3y)
!$OMP DECLARE TARGET(bcfs3y)
!$OMP DECLARE TARGET(acfs4y)
!$OMP DECLARE TARGET(bcfs4y)
!$OMP DECLARE TARGET(ccfs4y)
!$OMP DECLARE TARGET(acfs5y)
!$OMP DECLARE TARGET(bcfs5y)
!$OMP DECLARE TARGET(ccfs5y)
!$OMP DECLARE TARGET(dcfs5y)
!$OMP DECLARE TARGET(ovdly2)
!$OMP DECLARE TARGET(acoffz)
!$OMP DECLARE TARGET(bcoffz)
!$OMP DECLARE TARGET(ccoffz)
!$OMP DECLARE TARGET(dcoffz)
!$OMP DECLARE TARGET(ecoffz)
!$OMP DECLARE TARGET(acof1z)
!$OMP DECLARE TARGET(bcof1z)
!$OMP DECLARE TARGET(ccof1z)
!$OMP DECLARE TARGET(dcof1z)
!$OMP DECLARE TARGET(acof2z)
!$OMP DECLARE TARGET(bcof2z)
!$OMP DECLARE TARGET(ccof2z)
!$OMP DECLARE TARGET(dcof2z)
!$OMP DECLARE TARGET(acof3z)
!$OMP DECLARE TARGET(bcof3z)
!$OMP DECLARE TARGET(acof4z)
!$OMP DECLARE TARGET(bcof4z)
!$OMP DECLARE TARGET(ccof4z)
!$OMP DECLARE TARGET(acof5z)
!$OMP DECLARE TARGET(bcof5z)
!$OMP DECLARE TARGET(ccof5z)
!$OMP DECLARE TARGET(dcof5z)
!$OMP DECLARE TARGET(ovdelz)
!$OMP DECLARE TARGET(acofsz)
!$OMP DECLARE TARGET(bcofsz)
!$OMP DECLARE TARGET(ccofsz)
!$OMP DECLARE TARGET(dcofsz)
!$OMP DECLARE TARGET(ecofsz)
!$OMP DECLARE TARGET(acfs1z)
!$OMP DECLARE TARGET(bcfs1z)
!$OMP DECLARE TARGET(ccfs1z)
!$OMP DECLARE TARGET(dcfs1z)
!$OMP DECLARE TARGET(ecfs1z)
!$OMP DECLARE TARGET(acfs2z)
!$OMP DECLARE TARGET(bcfs2z)
!$OMP DECLARE TARGET(ccfs2z)
!$OMP DECLARE TARGET(dcfs2z)
!$OMP DECLARE TARGET(ecfs2z)
!$OMP DECLARE TARGET(acfs3z)
!$OMP DECLARE TARGET(bcfs3z)
!$OMP DECLARE TARGET(acfs4z)
!$OMP DECLARE TARGET(bcfs4z)
!$OMP DECLARE TARGET(ccfs4z)
!$OMP DECLARE TARGET(acfs5z)
!$OMP DECLARE TARGET(bcfs5z)
!$OMP DECLARE TARGET(ccfs5z)
!$OMP DECLARE TARGET(dcfs5z)
!$OMP DECLARE TARGET(ovdlz2)
!$OMP DECLARE TARGET(acofx1)
!$OMP DECLARE TARGET(bcofx1)
!$OMP DECLARE TARGET(acofy1)
!$OMP DECLARE TARGET(bcofy1)
!$OMP DECLARE TARGET(acofz1)
!$OMP DECLARE TARGET(bcofz1)
!$OMP DECLARE TARGET(acofxy)
!$OMP DECLARE TARGET(bcofxy)
!$OMP DECLARE TARGET(ccofxy)
!$OMP DECLARE TARGET(dcofxy)
!$OMP DECLARE TARGET(ecofxy)
!$OMP DECLARE TARGET(acf1xy)
!$OMP DECLARE TARGET(bcf1xy)
!$OMP DECLARE TARGET(ccf1xy)
!$OMP DECLARE TARGET(dcf1xy)
!$OMP DECLARE TARGET(acf2xy)
!$OMP DECLARE TARGET(bcf2xy)
!$OMP DECLARE TARGET(ccf2xy)
!$OMP DECLARE TARGET(dcf2xy)
!$OMP DECLARE TARGET(acf3xy)
!$OMP DECLARE TARGET(bcf3xy)
!$OMP DECLARE TARGET(acf4xy)
!$OMP DECLARE TARGET(bcf4xy)
!$OMP DECLARE TARGET(ccf4xy)
!$OMP DECLARE TARGET(acf5xy)
!$OMP DECLARE TARGET(bcf5xy)
!$OMP DECLARE TARGET(ccf5xy)
!$OMP DECLARE TARGET(dcf5xy)
!$OMP DECLARE TARGET(acc1xy)
!$OMP DECLARE TARGET(bcc1xy)
!$OMP DECLARE TARGET(ccc1xy)
!$OMP DECLARE TARGET(dcc1xy)
!$OMP DECLARE TARGET(acc2xy)
!$OMP DECLARE TARGET(bcc2xy)
!$OMP DECLARE TARGET(ccc2xy)
!$OMP DECLARE TARGET(dcc2xy)
!$OMP DECLARE TARGET(acofxz)
!$OMP DECLARE TARGET(bcofxz)
!$OMP DECLARE TARGET(ccofxz)
!$OMP DECLARE TARGET(dcofxz)
!$OMP DECLARE TARGET(ecofxz)
!$OMP DECLARE TARGET(acf1xz)
!$OMP DECLARE TARGET(bcf1xz)
!$OMP DECLARE TARGET(ccf1xz)
!$OMP DECLARE TARGET(dcf1xz)
!$OMP DECLARE TARGET(acf2xz)
!$OMP DECLARE TARGET(bcf2xz)
!$OMP DECLARE TARGET(ccf2xz)
!$OMP DECLARE TARGET(dcf2xz)
!$OMP DECLARE TARGET(acf3xz)
!$OMP DECLARE TARGET(bcf3xz)
!$OMP DECLARE TARGET(acf4xz)
!$OMP DECLARE TARGET(bcf4xz)
!$OMP DECLARE TARGET(ccf4xz)
!$OMP DECLARE TARGET(acf5xz)
!$OMP DECLARE TARGET(bcf5xz)
!$OMP DECLARE TARGET(ccf5xz)
!$OMP DECLARE TARGET(dcf5xz)
!$OMP DECLARE TARGET(acc1xz)
!$OMP DECLARE TARGET(bcc1xz)
!$OMP DECLARE TARGET(ccc1xz)
!$OMP DECLARE TARGET(dcc1xz)
!$OMP DECLARE TARGET(acc2xz)
!$OMP DECLARE TARGET(bcc2xz)
!$OMP DECLARE TARGET(ccc2xz)
!$OMP DECLARE TARGET(dcc2xz)
!$OMP DECLARE TARGET(acofyz)
!$OMP DECLARE TARGET(bcofyz)
!$OMP DECLARE TARGET(ccofyz)
!$OMP DECLARE TARGET(dcofyz)
!$OMP DECLARE TARGET(ecofyz)
!$OMP DECLARE TARGET(acf1yz)
!$OMP DECLARE TARGET(bcf1yz)
!$OMP DECLARE TARGET(ccf1yz)
!$OMP DECLARE TARGET(dcf1yz)
!$OMP DECLARE TARGET(acf2yz)
!$OMP DECLARE TARGET(bcf2yz)
!$OMP DECLARE TARGET(ccf2yz)
!$OMP DECLARE TARGET(dcf2yz)
!$OMP DECLARE TARGET(acf3yz)
!$OMP DECLARE TARGET(bcf3yz)
!$OMP DECLARE TARGET(acf4yz)
!$OMP DECLARE TARGET(bcf4yz)
!$OMP DECLARE TARGET(ccf4yz)
!$OMP DECLARE TARGET(acf5yz)
!$OMP DECLARE TARGET(bcf5yz)
!$OMP DECLARE TARGET(ccf5yz)
!$OMP DECLARE TARGET(dcf5yz)
!$OMP DECLARE TARGET(acc1yz)
!$OMP DECLARE TARGET(bcc1yz)
!$OMP DECLARE TARGET(ccc1yz)
!$OMP DECLARE TARGET(dcc1yz)
!$OMP DECLARE TARGET(acc2yz)
!$OMP DECLARE TARGET(bcc2yz)
!$OMP DECLARE TARGET(ccc2yz)
!$OMP DECLARE TARGET(dcc2yz)
!$OMP DECLARE TARGET(foursb)
!$OMP DECLARE TARGET(trfrth)
!$OMP DECLARE TARGET(rlamda)
!$OMP DECLARE TARGET(alamda)
#endif
#endif

END MODULE OPS_CONSTANTS

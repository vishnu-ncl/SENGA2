// Auto-generated at 2026-04-28 18:44:59.830325 by ops-translator

// headers
#define OPS_3D
#define OPS_SOA
#define OPS_API 2
#include "ops_lib_core.h"

#ifdef OPS_MPI
#include "ops_mpi_core.h"
#include <limits>
#endif

#include "ops_f2c_prelude_v2.h"

namespace f2c = ops::f2c;

#include "ops_sycl_rt_support.h"
#include "ops_sycl_reduction.h"

// global constants
cl::sycl::buffer<int,1> *ncofmx_p = nullptr;
extern int ncofmx;
cl::sycl::buffer<int,1> *ntinmx_p = nullptr;
extern int ntinmx;
cl::sycl::buffer<int,1> *nspcmx_p = nullptr;
extern int nspcmx;
cl::sycl::buffer<int,1> *nssmax_p = nullptr;
extern int nssmax;
cl::sycl::buffer<int,1> *nstpmx_p = nullptr;
extern int nstpmx;
cl::sycl::buffer<int,1> *ndcfmx_p = nullptr;
extern int ndcfmx;
cl::sycl::buffer<int,1> *nvcfmx_p = nullptr;
extern int nvcfmx;
cl::sycl::buffer<int,1> *nccfmx_p = nullptr;
extern int nccfmx;
cl::sycl::buffer<int,1> *nrkmax_p = nullptr;
extern int nrkmax;
cl::sycl::buffer<int,1> *ncbcsz_p = nullptr;
extern int ncbcsz;
cl::sycl::buffer<int,1> *nbcprr_p = nullptr;
extern int nbcprr;
cl::sycl::buffer<int,1> *nspimx_p = nullptr;
extern int nspimx;
cl::sycl::buffer<int,1> *ntbase_p = nullptr;
extern int ntbase;
cl::sycl::buffer<int,1> *nintmx_p = nullptr;
extern int nintmx;
cl::sycl::buffer<int,1> *nctmax_p = nullptr;
extern int nctmax;
cl::sycl::buffer<int,1> *nctmm1_p = nullptr;
extern int nctmm1;
cl::sycl::buffer<int,1> *nrsmax_p = nullptr;
extern int nrsmax;
cl::sycl::buffer<int,1> *nbcpri_p = nullptr;
extern int nbcpri;
cl::sycl::buffer<int,1> *ncfrmx_p = nullptr;
extern int ncfrmx;
cl::sycl::buffer<double,1> *acoffx_p = nullptr;
extern double acoffx;
cl::sycl::buffer<double,1> *bcoffx_p = nullptr;
extern double bcoffx;
cl::sycl::buffer<double,1> *ccoffx_p = nullptr;
extern double ccoffx;
cl::sycl::buffer<double,1> *dcoffx_p = nullptr;
extern double dcoffx;
cl::sycl::buffer<double,1> *ecoffx_p = nullptr;
extern double ecoffx;
cl::sycl::buffer<double,1> *acof1x_p = nullptr;
extern double acof1x;
cl::sycl::buffer<double,1> *bcof1x_p = nullptr;
extern double bcof1x;
cl::sycl::buffer<double,1> *ccof1x_p = nullptr;
extern double ccof1x;
cl::sycl::buffer<double,1> *dcof1x_p = nullptr;
extern double dcof1x;
cl::sycl::buffer<double,1> *acof2x_p = nullptr;
extern double acof2x;
cl::sycl::buffer<double,1> *bcof2x_p = nullptr;
extern double bcof2x;
cl::sycl::buffer<double,1> *ccof2x_p = nullptr;
extern double ccof2x;
cl::sycl::buffer<double,1> *dcof2x_p = nullptr;
extern double dcof2x;
cl::sycl::buffer<double,1> *acof3x_p = nullptr;
extern double acof3x;
cl::sycl::buffer<double,1> *bcof3x_p = nullptr;
extern double bcof3x;
cl::sycl::buffer<double,1> *acof4x_p = nullptr;
extern double acof4x;
cl::sycl::buffer<double,1> *bcof4x_p = nullptr;
extern double bcof4x;
cl::sycl::buffer<double,1> *ccof4x_p = nullptr;
extern double ccof4x;
cl::sycl::buffer<double,1> *acof5x_p = nullptr;
extern double acof5x;
cl::sycl::buffer<double,1> *bcof5x_p = nullptr;
extern double bcof5x;
cl::sycl::buffer<double,1> *ccof5x_p = nullptr;
extern double ccof5x;
cl::sycl::buffer<double,1> *dcof5x_p = nullptr;
extern double dcof5x;
cl::sycl::buffer<double,1> *ovdelx_p = nullptr;
extern double ovdelx;
cl::sycl::buffer<double,1> *acofsx_p = nullptr;
extern double acofsx;
cl::sycl::buffer<double,1> *bcofsx_p = nullptr;
extern double bcofsx;
cl::sycl::buffer<double,1> *ccofsx_p = nullptr;
extern double ccofsx;
cl::sycl::buffer<double,1> *dcofsx_p = nullptr;
extern double dcofsx;
cl::sycl::buffer<double,1> *ecofsx_p = nullptr;
extern double ecofsx;
cl::sycl::buffer<double,1> *acfs1x_p = nullptr;
extern double acfs1x;
cl::sycl::buffer<double,1> *bcfs1x_p = nullptr;
extern double bcfs1x;
cl::sycl::buffer<double,1> *ccfs1x_p = nullptr;
extern double ccfs1x;
cl::sycl::buffer<double,1> *dcfs1x_p = nullptr;
extern double dcfs1x;
cl::sycl::buffer<double,1> *ecfs1x_p = nullptr;
extern double ecfs1x;
cl::sycl::buffer<double,1> *acfs2x_p = nullptr;
extern double acfs2x;
cl::sycl::buffer<double,1> *bcfs2x_p = nullptr;
extern double bcfs2x;
cl::sycl::buffer<double,1> *ccfs2x_p = nullptr;
extern double ccfs2x;
cl::sycl::buffer<double,1> *dcfs2x_p = nullptr;
extern double dcfs2x;
cl::sycl::buffer<double,1> *ecfs2x_p = nullptr;
extern double ecfs2x;
cl::sycl::buffer<double,1> *acfs3x_p = nullptr;
extern double acfs3x;
cl::sycl::buffer<double,1> *bcfs3x_p = nullptr;
extern double bcfs3x;
cl::sycl::buffer<double,1> *acfs4x_p = nullptr;
extern double acfs4x;
cl::sycl::buffer<double,1> *bcfs4x_p = nullptr;
extern double bcfs4x;
cl::sycl::buffer<double,1> *ccfs4x_p = nullptr;
extern double ccfs4x;
cl::sycl::buffer<double,1> *acfs5x_p = nullptr;
extern double acfs5x;
cl::sycl::buffer<double,1> *bcfs5x_p = nullptr;
extern double bcfs5x;
cl::sycl::buffer<double,1> *ccfs5x_p = nullptr;
extern double ccfs5x;
cl::sycl::buffer<double,1> *dcfs5x_p = nullptr;
extern double dcfs5x;
cl::sycl::buffer<double,1> *ovdlx2_p = nullptr;
extern double ovdlx2;
cl::sycl::buffer<double,1> *acoffy_p = nullptr;
extern double acoffy;
cl::sycl::buffer<double,1> *bcoffy_p = nullptr;
extern double bcoffy;
cl::sycl::buffer<double,1> *ccoffy_p = nullptr;
extern double ccoffy;
cl::sycl::buffer<double,1> *dcoffy_p = nullptr;
extern double dcoffy;
cl::sycl::buffer<double,1> *ecoffy_p = nullptr;
extern double ecoffy;
cl::sycl::buffer<double,1> *acof1y_p = nullptr;
extern double acof1y;
cl::sycl::buffer<double,1> *bcof1y_p = nullptr;
extern double bcof1y;
cl::sycl::buffer<double,1> *ccof1y_p = nullptr;
extern double ccof1y;
cl::sycl::buffer<double,1> *dcof1y_p = nullptr;
extern double dcof1y;
cl::sycl::buffer<double,1> *acof2y_p = nullptr;
extern double acof2y;
cl::sycl::buffer<double,1> *bcof2y_p = nullptr;
extern double bcof2y;
cl::sycl::buffer<double,1> *ccof2y_p = nullptr;
extern double ccof2y;
cl::sycl::buffer<double,1> *dcof2y_p = nullptr;
extern double dcof2y;
cl::sycl::buffer<double,1> *acof3y_p = nullptr;
extern double acof3y;
cl::sycl::buffer<double,1> *bcof3y_p = nullptr;
extern double bcof3y;
cl::sycl::buffer<double,1> *acof4y_p = nullptr;
extern double acof4y;
cl::sycl::buffer<double,1> *bcof4y_p = nullptr;
extern double bcof4y;
cl::sycl::buffer<double,1> *ccof4y_p = nullptr;
extern double ccof4y;
cl::sycl::buffer<double,1> *acof5y_p = nullptr;
extern double acof5y;
cl::sycl::buffer<double,1> *bcof5y_p = nullptr;
extern double bcof5y;
cl::sycl::buffer<double,1> *ccof5y_p = nullptr;
extern double ccof5y;
cl::sycl::buffer<double,1> *dcof5y_p = nullptr;
extern double dcof5y;
cl::sycl::buffer<double,1> *ovdely_p = nullptr;
extern double ovdely;
cl::sycl::buffer<double,1> *acofsy_p = nullptr;
extern double acofsy;
cl::sycl::buffer<double,1> *bcofsy_p = nullptr;
extern double bcofsy;
cl::sycl::buffer<double,1> *ccofsy_p = nullptr;
extern double ccofsy;
cl::sycl::buffer<double,1> *dcofsy_p = nullptr;
extern double dcofsy;
cl::sycl::buffer<double,1> *ecofsy_p = nullptr;
extern double ecofsy;
cl::sycl::buffer<double,1> *acfs1y_p = nullptr;
extern double acfs1y;
cl::sycl::buffer<double,1> *bcfs1y_p = nullptr;
extern double bcfs1y;
cl::sycl::buffer<double,1> *ccfs1y_p = nullptr;
extern double ccfs1y;
cl::sycl::buffer<double,1> *dcfs1y_p = nullptr;
extern double dcfs1y;
cl::sycl::buffer<double,1> *ecfs1y_p = nullptr;
extern double ecfs1y;
cl::sycl::buffer<double,1> *acfs2y_p = nullptr;
extern double acfs2y;
cl::sycl::buffer<double,1> *bcfs2y_p = nullptr;
extern double bcfs2y;
cl::sycl::buffer<double,1> *ccfs2y_p = nullptr;
extern double ccfs2y;
cl::sycl::buffer<double,1> *dcfs2y_p = nullptr;
extern double dcfs2y;
cl::sycl::buffer<double,1> *ecfs2y_p = nullptr;
extern double ecfs2y;
cl::sycl::buffer<double,1> *acfs3y_p = nullptr;
extern double acfs3y;
cl::sycl::buffer<double,1> *bcfs3y_p = nullptr;
extern double bcfs3y;
cl::sycl::buffer<double,1> *acfs4y_p = nullptr;
extern double acfs4y;
cl::sycl::buffer<double,1> *bcfs4y_p = nullptr;
extern double bcfs4y;
cl::sycl::buffer<double,1> *ccfs4y_p = nullptr;
extern double ccfs4y;
cl::sycl::buffer<double,1> *acfs5y_p = nullptr;
extern double acfs5y;
cl::sycl::buffer<double,1> *bcfs5y_p = nullptr;
extern double bcfs5y;
cl::sycl::buffer<double,1> *ccfs5y_p = nullptr;
extern double ccfs5y;
cl::sycl::buffer<double,1> *dcfs5y_p = nullptr;
extern double dcfs5y;
cl::sycl::buffer<double,1> *ovdly2_p = nullptr;
extern double ovdly2;
cl::sycl::buffer<double,1> *acoffz_p = nullptr;
extern double acoffz;
cl::sycl::buffer<double,1> *bcoffz_p = nullptr;
extern double bcoffz;
cl::sycl::buffer<double,1> *ccoffz_p = nullptr;
extern double ccoffz;
cl::sycl::buffer<double,1> *dcoffz_p = nullptr;
extern double dcoffz;
cl::sycl::buffer<double,1> *ecoffz_p = nullptr;
extern double ecoffz;
cl::sycl::buffer<double,1> *acof1z_p = nullptr;
extern double acof1z;
cl::sycl::buffer<double,1> *bcof1z_p = nullptr;
extern double bcof1z;
cl::sycl::buffer<double,1> *ccof1z_p = nullptr;
extern double ccof1z;
cl::sycl::buffer<double,1> *dcof1z_p = nullptr;
extern double dcof1z;
cl::sycl::buffer<double,1> *acof2z_p = nullptr;
extern double acof2z;
cl::sycl::buffer<double,1> *bcof2z_p = nullptr;
extern double bcof2z;
cl::sycl::buffer<double,1> *ccof2z_p = nullptr;
extern double ccof2z;
cl::sycl::buffer<double,1> *dcof2z_p = nullptr;
extern double dcof2z;
cl::sycl::buffer<double,1> *acof3z_p = nullptr;
extern double acof3z;
cl::sycl::buffer<double,1> *bcof3z_p = nullptr;
extern double bcof3z;
cl::sycl::buffer<double,1> *acof4z_p = nullptr;
extern double acof4z;
cl::sycl::buffer<double,1> *bcof4z_p = nullptr;
extern double bcof4z;
cl::sycl::buffer<double,1> *ccof4z_p = nullptr;
extern double ccof4z;
cl::sycl::buffer<double,1> *acof5z_p = nullptr;
extern double acof5z;
cl::sycl::buffer<double,1> *bcof5z_p = nullptr;
extern double bcof5z;
cl::sycl::buffer<double,1> *ccof5z_p = nullptr;
extern double ccof5z;
cl::sycl::buffer<double,1> *dcof5z_p = nullptr;
extern double dcof5z;
cl::sycl::buffer<double,1> *ovdelz_p = nullptr;
extern double ovdelz;
cl::sycl::buffer<double,1> *acofsz_p = nullptr;
extern double acofsz;
cl::sycl::buffer<double,1> *bcofsz_p = nullptr;
extern double bcofsz;
cl::sycl::buffer<double,1> *ccofsz_p = nullptr;
extern double ccofsz;
cl::sycl::buffer<double,1> *dcofsz_p = nullptr;
extern double dcofsz;
cl::sycl::buffer<double,1> *ecofsz_p = nullptr;
extern double ecofsz;
cl::sycl::buffer<double,1> *acfs1z_p = nullptr;
extern double acfs1z;
cl::sycl::buffer<double,1> *bcfs1z_p = nullptr;
extern double bcfs1z;
cl::sycl::buffer<double,1> *ccfs1z_p = nullptr;
extern double ccfs1z;
cl::sycl::buffer<double,1> *dcfs1z_p = nullptr;
extern double dcfs1z;
cl::sycl::buffer<double,1> *ecfs1z_p = nullptr;
extern double ecfs1z;
cl::sycl::buffer<double,1> *acfs2z_p = nullptr;
extern double acfs2z;
cl::sycl::buffer<double,1> *bcfs2z_p = nullptr;
extern double bcfs2z;
cl::sycl::buffer<double,1> *ccfs2z_p = nullptr;
extern double ccfs2z;
cl::sycl::buffer<double,1> *dcfs2z_p = nullptr;
extern double dcfs2z;
cl::sycl::buffer<double,1> *ecfs2z_p = nullptr;
extern double ecfs2z;
cl::sycl::buffer<double,1> *acfs3z_p = nullptr;
extern double acfs3z;
cl::sycl::buffer<double,1> *bcfs3z_p = nullptr;
extern double bcfs3z;
cl::sycl::buffer<double,1> *acfs4z_p = nullptr;
extern double acfs4z;
cl::sycl::buffer<double,1> *bcfs4z_p = nullptr;
extern double bcfs4z;
cl::sycl::buffer<double,1> *ccfs4z_p = nullptr;
extern double ccfs4z;
cl::sycl::buffer<double,1> *acfs5z_p = nullptr;
extern double acfs5z;
cl::sycl::buffer<double,1> *bcfs5z_p = nullptr;
extern double bcfs5z;
cl::sycl::buffer<double,1> *ccfs5z_p = nullptr;
extern double ccfs5z;
cl::sycl::buffer<double,1> *dcfs5z_p = nullptr;
extern double dcfs5z;
cl::sycl::buffer<double,1> *ovdlz2_p = nullptr;
extern double ovdlz2;
cl::sycl::buffer<double,1> *acofx1_p = nullptr;
extern double acofx1;
cl::sycl::buffer<double,1> *bcofx1_p = nullptr;
extern double bcofx1;
cl::sycl::buffer<double,1> *acofy1_p = nullptr;
extern double acofy1;
cl::sycl::buffer<double,1> *bcofy1_p = nullptr;
extern double bcofy1;
cl::sycl::buffer<double,1> *acofz1_p = nullptr;
extern double acofz1;
cl::sycl::buffer<double,1> *bcofz1_p = nullptr;
extern double bcofz1;
cl::sycl::buffer<double,1> *acofxy_p = nullptr;
extern double acofxy;
cl::sycl::buffer<double,1> *bcofxy_p = nullptr;
extern double bcofxy;
cl::sycl::buffer<double,1> *ccofxy_p = nullptr;
extern double ccofxy;
cl::sycl::buffer<double,1> *dcofxy_p = nullptr;
extern double dcofxy;
cl::sycl::buffer<double,1> *ecofxy_p = nullptr;
extern double ecofxy;
cl::sycl::buffer<double,1> *acf1xy_p = nullptr;
extern double acf1xy;
cl::sycl::buffer<double,1> *bcf1xy_p = nullptr;
extern double bcf1xy;
cl::sycl::buffer<double,1> *ccf1xy_p = nullptr;
extern double ccf1xy;
cl::sycl::buffer<double,1> *dcf1xy_p = nullptr;
extern double dcf1xy;
cl::sycl::buffer<double,1> *acf2xy_p = nullptr;
extern double acf2xy;
cl::sycl::buffer<double,1> *bcf2xy_p = nullptr;
extern double bcf2xy;
cl::sycl::buffer<double,1> *ccf2xy_p = nullptr;
extern double ccf2xy;
cl::sycl::buffer<double,1> *dcf2xy_p = nullptr;
extern double dcf2xy;
cl::sycl::buffer<double,1> *acf3xy_p = nullptr;
extern double acf3xy;
cl::sycl::buffer<double,1> *bcf3xy_p = nullptr;
extern double bcf3xy;
cl::sycl::buffer<double,1> *acf4xy_p = nullptr;
extern double acf4xy;
cl::sycl::buffer<double,1> *bcf4xy_p = nullptr;
extern double bcf4xy;
cl::sycl::buffer<double,1> *ccf4xy_p = nullptr;
extern double ccf4xy;
cl::sycl::buffer<double,1> *acf5xy_p = nullptr;
extern double acf5xy;
cl::sycl::buffer<double,1> *bcf5xy_p = nullptr;
extern double bcf5xy;
cl::sycl::buffer<double,1> *ccf5xy_p = nullptr;
extern double ccf5xy;
cl::sycl::buffer<double,1> *dcf5xy_p = nullptr;
extern double dcf5xy;
cl::sycl::buffer<double,1> *acc1xy_p = nullptr;
extern double acc1xy;
cl::sycl::buffer<double,1> *bcc1xy_p = nullptr;
extern double bcc1xy;
cl::sycl::buffer<double,1> *ccc1xy_p = nullptr;
extern double ccc1xy;
cl::sycl::buffer<double,1> *dcc1xy_p = nullptr;
extern double dcc1xy;
cl::sycl::buffer<double,1> *acc2xy_p = nullptr;
extern double acc2xy;
cl::sycl::buffer<double,1> *bcc2xy_p = nullptr;
extern double bcc2xy;
cl::sycl::buffer<double,1> *ccc2xy_p = nullptr;
extern double ccc2xy;
cl::sycl::buffer<double,1> *dcc2xy_p = nullptr;
extern double dcc2xy;
cl::sycl::buffer<double,1> *acofxz_p = nullptr;
extern double acofxz;
cl::sycl::buffer<double,1> *bcofxz_p = nullptr;
extern double bcofxz;
cl::sycl::buffer<double,1> *ccofxz_p = nullptr;
extern double ccofxz;
cl::sycl::buffer<double,1> *dcofxz_p = nullptr;
extern double dcofxz;
cl::sycl::buffer<double,1> *ecofxz_p = nullptr;
extern double ecofxz;
cl::sycl::buffer<double,1> *acf1xz_p = nullptr;
extern double acf1xz;
cl::sycl::buffer<double,1> *bcf1xz_p = nullptr;
extern double bcf1xz;
cl::sycl::buffer<double,1> *ccf1xz_p = nullptr;
extern double ccf1xz;
cl::sycl::buffer<double,1> *dcf1xz_p = nullptr;
extern double dcf1xz;
cl::sycl::buffer<double,1> *acf2xz_p = nullptr;
extern double acf2xz;
cl::sycl::buffer<double,1> *bcf2xz_p = nullptr;
extern double bcf2xz;
cl::sycl::buffer<double,1> *ccf2xz_p = nullptr;
extern double ccf2xz;
cl::sycl::buffer<double,1> *dcf2xz_p = nullptr;
extern double dcf2xz;
cl::sycl::buffer<double,1> *acf3xz_p = nullptr;
extern double acf3xz;
cl::sycl::buffer<double,1> *bcf3xz_p = nullptr;
extern double bcf3xz;
cl::sycl::buffer<double,1> *acf4xz_p = nullptr;
extern double acf4xz;
cl::sycl::buffer<double,1> *bcf4xz_p = nullptr;
extern double bcf4xz;
cl::sycl::buffer<double,1> *ccf4xz_p = nullptr;
extern double ccf4xz;
cl::sycl::buffer<double,1> *acf5xz_p = nullptr;
extern double acf5xz;
cl::sycl::buffer<double,1> *bcf5xz_p = nullptr;
extern double bcf5xz;
cl::sycl::buffer<double,1> *ccf5xz_p = nullptr;
extern double ccf5xz;
cl::sycl::buffer<double,1> *dcf5xz_p = nullptr;
extern double dcf5xz;
cl::sycl::buffer<double,1> *acc1xz_p = nullptr;
extern double acc1xz;
cl::sycl::buffer<double,1> *bcc1xz_p = nullptr;
extern double bcc1xz;
cl::sycl::buffer<double,1> *ccc1xz_p = nullptr;
extern double ccc1xz;
cl::sycl::buffer<double,1> *dcc1xz_p = nullptr;
extern double dcc1xz;
cl::sycl::buffer<double,1> *acc2xz_p = nullptr;
extern double acc2xz;
cl::sycl::buffer<double,1> *bcc2xz_p = nullptr;
extern double bcc2xz;
cl::sycl::buffer<double,1> *ccc2xz_p = nullptr;
extern double ccc2xz;
cl::sycl::buffer<double,1> *dcc2xz_p = nullptr;
extern double dcc2xz;
cl::sycl::buffer<double,1> *acofyz_p = nullptr;
extern double acofyz;
cl::sycl::buffer<double,1> *bcofyz_p = nullptr;
extern double bcofyz;
cl::sycl::buffer<double,1> *ccofyz_p = nullptr;
extern double ccofyz;
cl::sycl::buffer<double,1> *dcofyz_p = nullptr;
extern double dcofyz;
cl::sycl::buffer<double,1> *ecofyz_p = nullptr;
extern double ecofyz;
cl::sycl::buffer<double,1> *acf1yz_p = nullptr;
extern double acf1yz;
cl::sycl::buffer<double,1> *bcf1yz_p = nullptr;
extern double bcf1yz;
cl::sycl::buffer<double,1> *ccf1yz_p = nullptr;
extern double ccf1yz;
cl::sycl::buffer<double,1> *dcf1yz_p = nullptr;
extern double dcf1yz;
cl::sycl::buffer<double,1> *acf2yz_p = nullptr;
extern double acf2yz;
cl::sycl::buffer<double,1> *bcf2yz_p = nullptr;
extern double bcf2yz;
cl::sycl::buffer<double,1> *ccf2yz_p = nullptr;
extern double ccf2yz;
cl::sycl::buffer<double,1> *dcf2yz_p = nullptr;
extern double dcf2yz;
cl::sycl::buffer<double,1> *acf3yz_p = nullptr;
extern double acf3yz;
cl::sycl::buffer<double,1> *bcf3yz_p = nullptr;
extern double bcf3yz;
cl::sycl::buffer<double,1> *acf4yz_p = nullptr;
extern double acf4yz;
cl::sycl::buffer<double,1> *bcf4yz_p = nullptr;
extern double bcf4yz;
cl::sycl::buffer<double,1> *ccf4yz_p = nullptr;
extern double ccf4yz;
cl::sycl::buffer<double,1> *acf5yz_p = nullptr;
extern double acf5yz;
cl::sycl::buffer<double,1> *bcf5yz_p = nullptr;
extern double bcf5yz;
cl::sycl::buffer<double,1> *ccf5yz_p = nullptr;
extern double ccf5yz;
cl::sycl::buffer<double,1> *dcf5yz_p = nullptr;
extern double dcf5yz;
cl::sycl::buffer<double,1> *acc1yz_p = nullptr;
extern double acc1yz;
cl::sycl::buffer<double,1> *bcc1yz_p = nullptr;
extern double bcc1yz;
cl::sycl::buffer<double,1> *ccc1yz_p = nullptr;
extern double ccc1yz;
cl::sycl::buffer<double,1> *dcc1yz_p = nullptr;
extern double dcc1yz;
cl::sycl::buffer<double,1> *acc2yz_p = nullptr;
extern double acc2yz;
cl::sycl::buffer<double,1> *bcc2yz_p = nullptr;
extern double bcc2yz;
cl::sycl::buffer<double,1> *ccc2yz_p = nullptr;
extern double ccc2yz;
cl::sycl::buffer<double,1> *dcc2yz_p = nullptr;
extern double dcc2yz;
cl::sycl::buffer<double,1> *foursb_p = nullptr;
extern double foursb;
cl::sycl::buffer<double,1> *trfrth_p = nullptr;
extern double trfrth;
cl::sycl::buffer<double,1> *rlamda_p = nullptr;
extern double rlamda;
cl::sycl::buffer<double,1> *alamda_p = nullptr;
extern double alamda;

void ops_init_backend(){}

extern "C"
void ops_decl_const_f2c(char const * name, int dim, int size, char * dat) {

    OPS_instance *instance = OPS_instance::getOPSInstance();
    ops_execute(instance);

    if(!strcmp(name, "ncofmx")) {
        if(ncofmx_p == nullptr) ncofmx_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ncofmx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ntinmx")) {
        if(ntinmx_p == nullptr) ntinmx_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ntinmx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nspcmx")) {
        if(nspcmx_p == nullptr) nspcmx_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nspcmx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nssmax")) {
        if(nssmax_p == nullptr) nssmax_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nssmax_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nstpmx")) {
        if(nstpmx_p == nullptr) nstpmx_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nstpmx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ndcfmx")) {
        if(ndcfmx_p == nullptr) ndcfmx_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ndcfmx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nvcfmx")) {
        if(nvcfmx_p == nullptr) nvcfmx_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nvcfmx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nccfmx")) {
        if(nccfmx_p == nullptr) nccfmx_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nccfmx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nrkmax")) {
        if(nrkmax_p == nullptr) nrkmax_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nrkmax_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ncbcsz")) {
        if(ncbcsz_p == nullptr) ncbcsz_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ncbcsz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nbcprr")) {
        if(nbcprr_p == nullptr) nbcprr_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nbcprr_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nspimx")) {
        if(nspimx_p == nullptr) nspimx_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nspimx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ntbase")) {
        if(ntbase_p == nullptr) ntbase_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ntbase_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nintmx")) {
        if(nintmx_p == nullptr) nintmx_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nintmx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nctmax")) {
        if(nctmax_p == nullptr) nctmax_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nctmax_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nctmm1")) {
        if(nctmm1_p == nullptr) nctmm1_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nctmm1_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nrsmax")) {
        if(nrsmax_p == nullptr) nrsmax_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nrsmax_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "nbcpri")) {
        if(nbcpri_p == nullptr) nbcpri_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*nbcpri_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ncfrmx")) {
        if(ncfrmx_p == nullptr) ncfrmx_p = new cl::sycl::buffer<int,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ncfrmx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((int*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acoffx")) {
        if(acoffx_p == nullptr) acoffx_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acoffx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcoffx")) {
        if(bcoffx_p == nullptr) bcoffx_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcoffx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccoffx")) {
        if(ccoffx_p == nullptr) ccoffx_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccoffx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcoffx")) {
        if(dcoffx_p == nullptr) dcoffx_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcoffx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecoffx")) {
        if(ecoffx_p == nullptr) ecoffx_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecoffx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof1x")) {
        if(acof1x_p == nullptr) acof1x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof1x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof1x")) {
        if(bcof1x_p == nullptr) bcof1x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof1x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof1x")) {
        if(ccof1x_p == nullptr) ccof1x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof1x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcof1x")) {
        if(dcof1x_p == nullptr) dcof1x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcof1x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof2x")) {
        if(acof2x_p == nullptr) acof2x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof2x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof2x")) {
        if(bcof2x_p == nullptr) bcof2x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof2x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof2x")) {
        if(ccof2x_p == nullptr) ccof2x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof2x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcof2x")) {
        if(dcof2x_p == nullptr) dcof2x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcof2x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof3x")) {
        if(acof3x_p == nullptr) acof3x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof3x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof3x")) {
        if(bcof3x_p == nullptr) bcof3x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof3x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof4x")) {
        if(acof4x_p == nullptr) acof4x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof4x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof4x")) {
        if(bcof4x_p == nullptr) bcof4x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof4x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof4x")) {
        if(ccof4x_p == nullptr) ccof4x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof4x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof5x")) {
        if(acof5x_p == nullptr) acof5x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof5x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof5x")) {
        if(bcof5x_p == nullptr) bcof5x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof5x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof5x")) {
        if(ccof5x_p == nullptr) ccof5x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof5x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcof5x")) {
        if(dcof5x_p == nullptr) dcof5x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcof5x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ovdelx")) {
        if(ovdelx_p == nullptr) ovdelx_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ovdelx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acofsx")) {
        if(acofsx_p == nullptr) acofsx_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acofsx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcofsx")) {
        if(bcofsx_p == nullptr) bcofsx_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcofsx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccofsx")) {
        if(ccofsx_p == nullptr) ccofsx_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccofsx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcofsx")) {
        if(dcofsx_p == nullptr) dcofsx_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcofsx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecofsx")) {
        if(ecofsx_p == nullptr) ecofsx_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecofsx_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs1x")) {
        if(acfs1x_p == nullptr) acfs1x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs1x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs1x")) {
        if(bcfs1x_p == nullptr) bcfs1x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs1x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs1x")) {
        if(ccfs1x_p == nullptr) ccfs1x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs1x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcfs1x")) {
        if(dcfs1x_p == nullptr) dcfs1x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcfs1x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecfs1x")) {
        if(ecfs1x_p == nullptr) ecfs1x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecfs1x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs2x")) {
        if(acfs2x_p == nullptr) acfs2x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs2x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs2x")) {
        if(bcfs2x_p == nullptr) bcfs2x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs2x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs2x")) {
        if(ccfs2x_p == nullptr) ccfs2x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs2x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcfs2x")) {
        if(dcfs2x_p == nullptr) dcfs2x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcfs2x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecfs2x")) {
        if(ecfs2x_p == nullptr) ecfs2x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecfs2x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs3x")) {
        if(acfs3x_p == nullptr) acfs3x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs3x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs3x")) {
        if(bcfs3x_p == nullptr) bcfs3x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs3x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs4x")) {
        if(acfs4x_p == nullptr) acfs4x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs4x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs4x")) {
        if(bcfs4x_p == nullptr) bcfs4x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs4x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs4x")) {
        if(ccfs4x_p == nullptr) ccfs4x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs4x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs5x")) {
        if(acfs5x_p == nullptr) acfs5x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs5x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs5x")) {
        if(bcfs5x_p == nullptr) bcfs5x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs5x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs5x")) {
        if(ccfs5x_p == nullptr) ccfs5x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs5x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcfs5x")) {
        if(dcfs5x_p == nullptr) dcfs5x_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcfs5x_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ovdlx2")) {
        if(ovdlx2_p == nullptr) ovdlx2_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ovdlx2_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acoffy")) {
        if(acoffy_p == nullptr) acoffy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acoffy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcoffy")) {
        if(bcoffy_p == nullptr) bcoffy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcoffy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccoffy")) {
        if(ccoffy_p == nullptr) ccoffy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccoffy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcoffy")) {
        if(dcoffy_p == nullptr) dcoffy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcoffy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecoffy")) {
        if(ecoffy_p == nullptr) ecoffy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecoffy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof1y")) {
        if(acof1y_p == nullptr) acof1y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof1y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof1y")) {
        if(bcof1y_p == nullptr) bcof1y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof1y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof1y")) {
        if(ccof1y_p == nullptr) ccof1y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof1y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcof1y")) {
        if(dcof1y_p == nullptr) dcof1y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcof1y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof2y")) {
        if(acof2y_p == nullptr) acof2y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof2y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof2y")) {
        if(bcof2y_p == nullptr) bcof2y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof2y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof2y")) {
        if(ccof2y_p == nullptr) ccof2y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof2y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcof2y")) {
        if(dcof2y_p == nullptr) dcof2y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcof2y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof3y")) {
        if(acof3y_p == nullptr) acof3y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof3y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof3y")) {
        if(bcof3y_p == nullptr) bcof3y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof3y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof4y")) {
        if(acof4y_p == nullptr) acof4y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof4y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof4y")) {
        if(bcof4y_p == nullptr) bcof4y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof4y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof4y")) {
        if(ccof4y_p == nullptr) ccof4y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof4y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof5y")) {
        if(acof5y_p == nullptr) acof5y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof5y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof5y")) {
        if(bcof5y_p == nullptr) bcof5y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof5y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof5y")) {
        if(ccof5y_p == nullptr) ccof5y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof5y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcof5y")) {
        if(dcof5y_p == nullptr) dcof5y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcof5y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ovdely")) {
        if(ovdely_p == nullptr) ovdely_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ovdely_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acofsy")) {
        if(acofsy_p == nullptr) acofsy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acofsy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcofsy")) {
        if(bcofsy_p == nullptr) bcofsy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcofsy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccofsy")) {
        if(ccofsy_p == nullptr) ccofsy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccofsy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcofsy")) {
        if(dcofsy_p == nullptr) dcofsy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcofsy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecofsy")) {
        if(ecofsy_p == nullptr) ecofsy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecofsy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs1y")) {
        if(acfs1y_p == nullptr) acfs1y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs1y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs1y")) {
        if(bcfs1y_p == nullptr) bcfs1y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs1y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs1y")) {
        if(ccfs1y_p == nullptr) ccfs1y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs1y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcfs1y")) {
        if(dcfs1y_p == nullptr) dcfs1y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcfs1y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecfs1y")) {
        if(ecfs1y_p == nullptr) ecfs1y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecfs1y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs2y")) {
        if(acfs2y_p == nullptr) acfs2y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs2y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs2y")) {
        if(bcfs2y_p == nullptr) bcfs2y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs2y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs2y")) {
        if(ccfs2y_p == nullptr) ccfs2y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs2y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcfs2y")) {
        if(dcfs2y_p == nullptr) dcfs2y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcfs2y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecfs2y")) {
        if(ecfs2y_p == nullptr) ecfs2y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecfs2y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs3y")) {
        if(acfs3y_p == nullptr) acfs3y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs3y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs3y")) {
        if(bcfs3y_p == nullptr) bcfs3y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs3y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs4y")) {
        if(acfs4y_p == nullptr) acfs4y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs4y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs4y")) {
        if(bcfs4y_p == nullptr) bcfs4y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs4y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs4y")) {
        if(ccfs4y_p == nullptr) ccfs4y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs4y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs5y")) {
        if(acfs5y_p == nullptr) acfs5y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs5y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs5y")) {
        if(bcfs5y_p == nullptr) bcfs5y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs5y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs5y")) {
        if(ccfs5y_p == nullptr) ccfs5y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs5y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcfs5y")) {
        if(dcfs5y_p == nullptr) dcfs5y_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcfs5y_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ovdly2")) {
        if(ovdly2_p == nullptr) ovdly2_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ovdly2_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acoffz")) {
        if(acoffz_p == nullptr) acoffz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acoffz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcoffz")) {
        if(bcoffz_p == nullptr) bcoffz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcoffz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccoffz")) {
        if(ccoffz_p == nullptr) ccoffz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccoffz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcoffz")) {
        if(dcoffz_p == nullptr) dcoffz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcoffz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecoffz")) {
        if(ecoffz_p == nullptr) ecoffz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecoffz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof1z")) {
        if(acof1z_p == nullptr) acof1z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof1z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof1z")) {
        if(bcof1z_p == nullptr) bcof1z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof1z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof1z")) {
        if(ccof1z_p == nullptr) ccof1z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof1z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcof1z")) {
        if(dcof1z_p == nullptr) dcof1z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcof1z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof2z")) {
        if(acof2z_p == nullptr) acof2z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof2z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof2z")) {
        if(bcof2z_p == nullptr) bcof2z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof2z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof2z")) {
        if(ccof2z_p == nullptr) ccof2z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof2z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcof2z")) {
        if(dcof2z_p == nullptr) dcof2z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcof2z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof3z")) {
        if(acof3z_p == nullptr) acof3z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof3z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof3z")) {
        if(bcof3z_p == nullptr) bcof3z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof3z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof4z")) {
        if(acof4z_p == nullptr) acof4z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof4z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof4z")) {
        if(bcof4z_p == nullptr) bcof4z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof4z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof4z")) {
        if(ccof4z_p == nullptr) ccof4z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof4z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acof5z")) {
        if(acof5z_p == nullptr) acof5z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acof5z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcof5z")) {
        if(bcof5z_p == nullptr) bcof5z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcof5z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccof5z")) {
        if(ccof5z_p == nullptr) ccof5z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccof5z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcof5z")) {
        if(dcof5z_p == nullptr) dcof5z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcof5z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ovdelz")) {
        if(ovdelz_p == nullptr) ovdelz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ovdelz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acofsz")) {
        if(acofsz_p == nullptr) acofsz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acofsz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcofsz")) {
        if(bcofsz_p == nullptr) bcofsz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcofsz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccofsz")) {
        if(ccofsz_p == nullptr) ccofsz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccofsz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcofsz")) {
        if(dcofsz_p == nullptr) dcofsz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcofsz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecofsz")) {
        if(ecofsz_p == nullptr) ecofsz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecofsz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs1z")) {
        if(acfs1z_p == nullptr) acfs1z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs1z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs1z")) {
        if(bcfs1z_p == nullptr) bcfs1z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs1z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs1z")) {
        if(ccfs1z_p == nullptr) ccfs1z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs1z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcfs1z")) {
        if(dcfs1z_p == nullptr) dcfs1z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcfs1z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecfs1z")) {
        if(ecfs1z_p == nullptr) ecfs1z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecfs1z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs2z")) {
        if(acfs2z_p == nullptr) acfs2z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs2z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs2z")) {
        if(bcfs2z_p == nullptr) bcfs2z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs2z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs2z")) {
        if(ccfs2z_p == nullptr) ccfs2z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs2z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcfs2z")) {
        if(dcfs2z_p == nullptr) dcfs2z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcfs2z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecfs2z")) {
        if(ecfs2z_p == nullptr) ecfs2z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecfs2z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs3z")) {
        if(acfs3z_p == nullptr) acfs3z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs3z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs3z")) {
        if(bcfs3z_p == nullptr) bcfs3z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs3z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs4z")) {
        if(acfs4z_p == nullptr) acfs4z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs4z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs4z")) {
        if(bcfs4z_p == nullptr) bcfs4z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs4z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs4z")) {
        if(ccfs4z_p == nullptr) ccfs4z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs4z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acfs5z")) {
        if(acfs5z_p == nullptr) acfs5z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acfs5z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcfs5z")) {
        if(bcfs5z_p == nullptr) bcfs5z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcfs5z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccfs5z")) {
        if(ccfs5z_p == nullptr) ccfs5z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccfs5z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcfs5z")) {
        if(dcfs5z_p == nullptr) dcfs5z_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcfs5z_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ovdlz2")) {
        if(ovdlz2_p == nullptr) ovdlz2_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ovdlz2_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acofx1")) {
        if(acofx1_p == nullptr) acofx1_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acofx1_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcofx1")) {
        if(bcofx1_p == nullptr) bcofx1_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcofx1_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acofy1")) {
        if(acofy1_p == nullptr) acofy1_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acofy1_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcofy1")) {
        if(bcofy1_p == nullptr) bcofy1_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcofy1_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acofz1")) {
        if(acofz1_p == nullptr) acofz1_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acofz1_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcofz1")) {
        if(bcofz1_p == nullptr) bcofz1_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcofz1_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acofxy")) {
        if(acofxy_p == nullptr) acofxy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acofxy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcofxy")) {
        if(bcofxy_p == nullptr) bcofxy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcofxy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccofxy")) {
        if(ccofxy_p == nullptr) ccofxy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccofxy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcofxy")) {
        if(dcofxy_p == nullptr) dcofxy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcofxy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecofxy")) {
        if(ecofxy_p == nullptr) ecofxy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecofxy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf1xy")) {
        if(acf1xy_p == nullptr) acf1xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf1xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf1xy")) {
        if(bcf1xy_p == nullptr) bcf1xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf1xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf1xy")) {
        if(ccf1xy_p == nullptr) ccf1xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf1xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcf1xy")) {
        if(dcf1xy_p == nullptr) dcf1xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcf1xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf2xy")) {
        if(acf2xy_p == nullptr) acf2xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf2xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf2xy")) {
        if(bcf2xy_p == nullptr) bcf2xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf2xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf2xy")) {
        if(ccf2xy_p == nullptr) ccf2xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf2xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcf2xy")) {
        if(dcf2xy_p == nullptr) dcf2xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcf2xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf3xy")) {
        if(acf3xy_p == nullptr) acf3xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf3xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf3xy")) {
        if(bcf3xy_p == nullptr) bcf3xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf3xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf4xy")) {
        if(acf4xy_p == nullptr) acf4xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf4xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf4xy")) {
        if(bcf4xy_p == nullptr) bcf4xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf4xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf4xy")) {
        if(ccf4xy_p == nullptr) ccf4xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf4xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf5xy")) {
        if(acf5xy_p == nullptr) acf5xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf5xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf5xy")) {
        if(bcf5xy_p == nullptr) bcf5xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf5xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf5xy")) {
        if(ccf5xy_p == nullptr) ccf5xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf5xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcf5xy")) {
        if(dcf5xy_p == nullptr) dcf5xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcf5xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acc1xy")) {
        if(acc1xy_p == nullptr) acc1xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acc1xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcc1xy")) {
        if(bcc1xy_p == nullptr) bcc1xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcc1xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccc1xy")) {
        if(ccc1xy_p == nullptr) ccc1xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccc1xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcc1xy")) {
        if(dcc1xy_p == nullptr) dcc1xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcc1xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acc2xy")) {
        if(acc2xy_p == nullptr) acc2xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acc2xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcc2xy")) {
        if(bcc2xy_p == nullptr) bcc2xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcc2xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccc2xy")) {
        if(ccc2xy_p == nullptr) ccc2xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccc2xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcc2xy")) {
        if(dcc2xy_p == nullptr) dcc2xy_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcc2xy_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acofxz")) {
        if(acofxz_p == nullptr) acofxz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acofxz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcofxz")) {
        if(bcofxz_p == nullptr) bcofxz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcofxz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccofxz")) {
        if(ccofxz_p == nullptr) ccofxz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccofxz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcofxz")) {
        if(dcofxz_p == nullptr) dcofxz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcofxz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecofxz")) {
        if(ecofxz_p == nullptr) ecofxz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecofxz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf1xz")) {
        if(acf1xz_p == nullptr) acf1xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf1xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf1xz")) {
        if(bcf1xz_p == nullptr) bcf1xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf1xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf1xz")) {
        if(ccf1xz_p == nullptr) ccf1xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf1xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcf1xz")) {
        if(dcf1xz_p == nullptr) dcf1xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcf1xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf2xz")) {
        if(acf2xz_p == nullptr) acf2xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf2xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf2xz")) {
        if(bcf2xz_p == nullptr) bcf2xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf2xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf2xz")) {
        if(ccf2xz_p == nullptr) ccf2xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf2xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcf2xz")) {
        if(dcf2xz_p == nullptr) dcf2xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcf2xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf3xz")) {
        if(acf3xz_p == nullptr) acf3xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf3xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf3xz")) {
        if(bcf3xz_p == nullptr) bcf3xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf3xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf4xz")) {
        if(acf4xz_p == nullptr) acf4xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf4xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf4xz")) {
        if(bcf4xz_p == nullptr) bcf4xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf4xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf4xz")) {
        if(ccf4xz_p == nullptr) ccf4xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf4xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf5xz")) {
        if(acf5xz_p == nullptr) acf5xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf5xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf5xz")) {
        if(bcf5xz_p == nullptr) bcf5xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf5xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf5xz")) {
        if(ccf5xz_p == nullptr) ccf5xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf5xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcf5xz")) {
        if(dcf5xz_p == nullptr) dcf5xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcf5xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acc1xz")) {
        if(acc1xz_p == nullptr) acc1xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acc1xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcc1xz")) {
        if(bcc1xz_p == nullptr) bcc1xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcc1xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccc1xz")) {
        if(ccc1xz_p == nullptr) ccc1xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccc1xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcc1xz")) {
        if(dcc1xz_p == nullptr) dcc1xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcc1xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acc2xz")) {
        if(acc2xz_p == nullptr) acc2xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acc2xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcc2xz")) {
        if(bcc2xz_p == nullptr) bcc2xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcc2xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccc2xz")) {
        if(ccc2xz_p == nullptr) ccc2xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccc2xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcc2xz")) {
        if(dcc2xz_p == nullptr) dcc2xz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcc2xz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acofyz")) {
        if(acofyz_p == nullptr) acofyz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acofyz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcofyz")) {
        if(bcofyz_p == nullptr) bcofyz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcofyz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccofyz")) {
        if(ccofyz_p == nullptr) ccofyz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccofyz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcofyz")) {
        if(dcofyz_p == nullptr) dcofyz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcofyz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ecofyz")) {
        if(ecofyz_p == nullptr) ecofyz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ecofyz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf1yz")) {
        if(acf1yz_p == nullptr) acf1yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf1yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf1yz")) {
        if(bcf1yz_p == nullptr) bcf1yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf1yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf1yz")) {
        if(ccf1yz_p == nullptr) ccf1yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf1yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcf1yz")) {
        if(dcf1yz_p == nullptr) dcf1yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcf1yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf2yz")) {
        if(acf2yz_p == nullptr) acf2yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf2yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf2yz")) {
        if(bcf2yz_p == nullptr) bcf2yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf2yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf2yz")) {
        if(ccf2yz_p == nullptr) ccf2yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf2yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcf2yz")) {
        if(dcf2yz_p == nullptr) dcf2yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcf2yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf3yz")) {
        if(acf3yz_p == nullptr) acf3yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf3yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf3yz")) {
        if(bcf3yz_p == nullptr) bcf3yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf3yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf4yz")) {
        if(acf4yz_p == nullptr) acf4yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf4yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf4yz")) {
        if(bcf4yz_p == nullptr) bcf4yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf4yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf4yz")) {
        if(ccf4yz_p == nullptr) ccf4yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf4yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acf5yz")) {
        if(acf5yz_p == nullptr) acf5yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acf5yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcf5yz")) {
        if(bcf5yz_p == nullptr) bcf5yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcf5yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccf5yz")) {
        if(ccf5yz_p == nullptr) ccf5yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccf5yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcf5yz")) {
        if(dcf5yz_p == nullptr) dcf5yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcf5yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acc1yz")) {
        if(acc1yz_p == nullptr) acc1yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acc1yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcc1yz")) {
        if(bcc1yz_p == nullptr) bcc1yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcc1yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccc1yz")) {
        if(ccc1yz_p == nullptr) ccc1yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccc1yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcc1yz")) {
        if(dcc1yz_p == nullptr) dcc1yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcc1yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "acc2yz")) {
        if(acc2yz_p == nullptr) acc2yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*acc2yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "bcc2yz")) {
        if(bcc2yz_p == nullptr) bcc2yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*bcc2yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "ccc2yz")) {
        if(ccc2yz_p == nullptr) ccc2yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*ccc2yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "dcc2yz")) {
        if(dcc2yz_p == nullptr) dcc2yz_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*dcc2yz_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "foursb")) {
        if(foursb_p == nullptr) foursb_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*foursb_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "trfrth")) {
        if(trfrth_p == nullptr) trfrth_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*trfrth_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "rlamda")) {
        if(rlamda_p == nullptr) rlamda_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*rlamda_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    if(!strcmp(name, "alamda")) {
        if(alamda_p == nullptr) alamda_p = new cl::sycl::buffer<double,1>(cl::sycl::range<1>(dim));
        auto accessor = (*alamda_p).get_access<cl::sycl::access::mode::write>();
        for(int d = 0; d < dim; d++) {
            accessor[d] = ((double*)dat)[d];
        }
    }
    else
    {
        throw OPSException(OPS_RUNTIME_ERROR, "error: unknown const name");
    }
}

// user kernel files
#include "set_zero_kernel_f2c_kernel.cpp"
#include "set_zero_kernel_xdir_f2c_kernel.cpp"
#include "set_zero_kernel_ydir_f2c_kernel.cpp"
#include "set_zero_kernel_zdir_f2c_kernel.cpp"
#include "set_zero_kernel_int_f2c_kernel.cpp"
#include "dfbydx_kernel_main_f2c_kernel.cpp"
#include "d2fdx2_kernel_main_f2c_kernel.cpp"
#include "dfbydy_kernel_null_f2c_kernel.cpp"
#include "dfbydy_kernel_main_f2c_kernel.cpp"
#include "dfbydz_kernel_null_f2c_kernel.cpp"
#include "dfbydz_kernel_main_f2c_kernel.cpp"
#include "d2fdy2_kernel_null_f2c_kernel.cpp"
#include "d2fdy2_kernel_main_f2c_kernel.cpp"
#include "d2fdz2_kernel_null_f2c_kernel.cpp"
#include "d2fdz2_kernel_main_f2c_kernel.cpp"
#include "d2fdxy_kernel_null_f2c_kernel.cpp"
#include "d2fdxy_kernel_interior_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqA_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqB_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqC_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqD_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqE_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqF_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqG_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqH_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqI_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqJ_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqK_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqL_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqM_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqN_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqO_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqP_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqQ_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqR_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqS_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqT_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqU_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqV_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqW_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqX_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqY_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqZ_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAA_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAB_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAC_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAD_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAE_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAF_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAG_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAH_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAI_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAJ_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAK_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAL_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAM_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAN_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAO_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAP_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAQ_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAR_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAS_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAT_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAU_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAV_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAW_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAX_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAY_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqAZ_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqBA_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqBB_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqBC_f2c_kernel.cpp"
#include "d2fdxy_kernel_eqBD_f2c_kernel.cpp"
#include "d2fdxy_kernel_scaling_f2c_kernel.cpp"
#include "d2fdxz_kernel_null_f2c_kernel.cpp"
#include "d2fdxz_kernel_interior_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqA_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqB_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqC_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqD_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqE_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqF_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqG_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqH_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqI_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqJ_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqK_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqL_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqM_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqN_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqO_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqP_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqQ_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqR_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqS_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqT_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqU_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqV_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqW_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqX_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqY_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqZ_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAA_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAB_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAC_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAD_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAE_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAF_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAG_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAH_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAI_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAJ_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAK_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAL_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAM_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAN_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAO_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAP_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAQ_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAR_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAS_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAT_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAU_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAV_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAW_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAX_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAY_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqAZ_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqBA_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqBB_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqBC_f2c_kernel.cpp"
#include "d2fdxz_kernel_eqBD_f2c_kernel.cpp"
#include "d2fdxz_kernel_scaling_f2c_kernel.cpp"
#include "d2fdyz_kernel_null_f2c_kernel.cpp"
#include "d2fdyz_kernel_interior_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqA_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqB_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqC_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqD_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqE_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqF_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqG_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqH_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqI_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqJ_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqK_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqL_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqM_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqN_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqO_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqP_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqQ_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqR_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqS_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqT_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqU_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqV_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqW_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqX_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqY_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqZ_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAA_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAB_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAC_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAD_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAE_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAF_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAG_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAH_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAI_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAJ_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAK_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAL_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAM_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAN_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAO_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAP_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAQ_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAR_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAS_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAT_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAU_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAV_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAW_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAX_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAY_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqAZ_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqBA_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqBB_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqBC_f2c_kernel.cpp"
#include "d2fdyz_kernel_eqBD_f2c_kernel.cpp"
#include "d2fdyz_kernel_scaling_f2c_kernel.cpp"
#include "boundary_kernel_CPandGAS_xdir_f2c_kernel.cpp"
#include "boundary_kernel_CPandGAS_ydir_f2c_kernel.cpp"
#include "boundary_kernel_CPandGAS_zdir_f2c_kernel.cpp"
#include "maths_kernel_eqX_f2c_kernel.cpp"
#include "maths_kernel_eqBR_xdir_f2c_kernel.cpp"
#include "maths_kernel_eqBR_ydir_f2c_kernel.cpp"
#include "maths_kernel_eqBR_zdir_f2c_kernel.cpp"
#include "maths_kernel_eqT_f2c_kernel.cpp"
#include "boundary_kernel_internalenergy_xdir_f2c_kernel.cpp"
#include "boundary_kernel_internalenergy_ydir_f2c_kernel.cpp"
#include "boundary_kernel_internalenergy_zdir_f2c_kernel.cpp"
#include "maths_kernel_eqW_f2c_kernel.cpp"
#include "maths_kernel_eqZ_f2c_kernel.cpp"
#include "maths_kernel_eqAS_f2c_kernel.cpp"
#include "boundary_kernel_temperature_xdir_f2c_kernel.cpp"
#include "boundary_kernel_temperature_ydir_f2c_kernel.cpp"
#include "boundary_kernel_temperature_zdir_f2c_kernel.cpp"
#include "maths_kernel_eqAV_f2c_kernel.cpp"
#include "maths_kernel_eqBC_f2c_kernel.cpp"
#include "maths_kernel_eqBD_f2c_kernel.cpp"
#include "maths_kernel_eqBE_f2c_kernel.cpp"
#include "maths_kernel_eqAH_f2c_kernel.cpp"
#include "hf_kernel_eqA_f2c_kernel.cpp"
#include "hf_kernel_eqB_f2c_kernel.cpp"
#include "hf_kernel_eqC_f2c_kernel.cpp"
#include "hf_kernel_eqD_f2c_kernel.cpp"
#include "hf_kernel_eqE_f2c_kernel.cpp"
#include "hf_kernel_eqF_f2c_kernel.cpp"
#include "maths_kernel_eqAF_f2c_kernel.cpp"
#include "copy_kernel_f2c_kernel.cpp"
#include "copy_kernel_xdir_f2c_kernel.cpp"
#include "copy_kernel_ydir_f2c_kernel.cpp"
#include "copy_kernel_zdir_f2c_kernel.cpp"
#include "maths_kernel_eqA_f2c_kernel.cpp"
#include "maths_kernel_eqAP_f2c_kernel.cpp"
#include "maths_kernel_eqAQ_f2c_kernel.cpp"
#include "boundary_kernel_mass_xdir_f2c_kernel.cpp"
#include "boundary_kernel_mass_ydir_f2c_kernel.cpp"
#include "boundary_kernel_mass_zdir_f2c_kernel.cpp"
#include "maths_kernel_eqAT_f2c_kernel.cpp"
#include "maths_kernel_eqH_f2c_kernel.cpp"
#include "copy_kernel_sdim_to_mdim_f2c_kernel.cpp"
#include "maths_kernel_eqBFG_f2c_kernel.cpp"
#include "maths_kernel_eqBH_f2c_kernel.cpp"
#include "maths_kernel_eqBS_f2c_kernel.cpp"
#include "maths_kernel_eqAA_fused_f2c_kernel.cpp"
#include "maths_kernel_eqBL_f2c_kernel.cpp"
#include "boundary_kernel_speciesH_xdir_f2c_kernel.cpp"
#include "boundary_kernel_speciesH_ydir_f2c_kernel.cpp"
#include "boundary_kernel_speciesH_zdir_f2c_kernel.cpp"
#include "maths_kernel_eqL_f2c_kernel.cpp"
#include "maths_kernel_eqBA_f2c_kernel.cpp"
#include "maths_kernel_eqAI_f2c_kernel.cpp"
#include "hf_kernel_eqS_f2c_kernel.cpp"
#include "hf_kernel_eqT_f2c_kernel.cpp"
#include "hf_kernel_eqU_f2c_kernel.cpp"
#include "hf_kernel_eqV_f2c_kernel.cpp"
#include "hf_kernel_eqW_f2c_kernel.cpp"
#include "hf_kernel_eqX_f2c_kernel.cpp"
#include "maths_kernel_eqBB_f2c_kernel.cpp"
#include "maths_kernel_eqAG_f2c_kernel.cpp"
#include "maths_kernel_eqO_f2c_kernel.cpp"
#include "maths_kernel_eqY_f2c_kernel.cpp"
#include "maths_kernel_eqM_f2c_kernel.cpp"
#include "hf_kernel_eqM_f2c_kernel.cpp"
#include "hf_kernel_eqN_f2c_kernel.cpp"
#include "hf_kernel_eqO_f2c_kernel.cpp"
#include "hf_kernel_eqP_f2c_kernel.cpp"
#include "hf_kernel_eqQ_f2c_kernel.cpp"
#include "hf_kernel_eqR_f2c_kernel.cpp"
#include "hf_kernel_eqG_f2c_kernel.cpp"
#include "hf_kernel_eqH_f2c_kernel.cpp"
#include "hf_kernel_eqI_f2c_kernel.cpp"
#include "hf_kernel_eqJ_f2c_kernel.cpp"
#include "hf_kernel_eqK_f2c_kernel.cpp"
#include "hf_kernel_eqL_f2c_kernel.cpp"
#include "maths_kernel_eqAM_f2c_kernel.cpp"
#include "maths_kernel_eqBIJK_f2c_kernel.cpp"
#include "maths_kernel_eqAN_f2c_kernel.cpp"
#include "boundary_kernel_density_xdir_f2c_kernel.cpp"
#include "boundary_kernel_density_ydir_f2c_kernel.cpp"
#include "boundary_kernel_density_zdir_f2c_kernel.cpp"
#include "maths_kernel_eqU_fused_f2c_kernel.cpp"
#include "boundary_kernel_velcomp_xdir_f2c_kernel.cpp"
#include "boundary_kernel_velcomp_ydir_f2c_kernel.cpp"
#include "boundary_kernel_velcomp_zdir_f2c_kernel.cpp"
#include "boundary_kernel_eqA_xdir_f2c_kernel.cpp"
#include "boundary_kernel_eqA_ydir_f2c_kernel.cpp"
#include "boundary_kernel_eqA_zdir_f2c_kernel.cpp"
#include "maths_kernel_eqQ_f2c_kernel.cpp"
#include "boundary_kernel_velderiv_xdir_f2c_kernel.cpp"
#include "boundary_kernel_velderiv_ydir_f2c_kernel.cpp"
#include "boundary_kernel_velderiv_zdir_f2c_kernel.cpp"
#include "maths_kernel_eqAB_f2c_kernel.cpp"
#include "maths_kernel_eqAO_f2c_kernel.cpp"
#include "boundary_kernel_pressure_xdir_f2c_kernel.cpp"
#include "boundary_kernel_pressure_ydir_f2c_kernel.cpp"
#include "boundary_kernel_pressure_zdir_f2c_kernel.cpp"
#include "maths_kernel_eqS_f2c_kernel.cpp"
#include "maths_kernel_eqAL_f2c_kernel.cpp"
#include "maths_kernel_eqAK_f2c_kernel.cpp"
#include "maths_kernel_eqG_f2c_kernel.cpp"
#include "maths_kernel_eqR_f2c_kernel.cpp"
#include "maths_kernel_eqN_f2c_kernel.cpp"
#include "maths_kernel_eqAD_f2c_kernel.cpp"
#include "maths_kernel_eqV_f2c_kernel.cpp"
#include "maths_kernel_eqAA_f2c_kernel.cpp"
#include "maths_kernel_eqAR_f2c_kernel.cpp"
#include "maths_kernel_eqI_f2c_kernel.cpp"
#include "maths_kernel_eqAJ_f2c_kernel.cpp"
#include "maths_kernel_eqAE_f2c_kernel.cpp"
#include "maths_kernel_eqAC_f2c_kernel.cpp"
#include "maths_kernel_eqtau_f2c_kernel.cpp"
#include "maths_kernel_eqC_f2c_kernel.cpp"
#include "lincom_kernel_main_f2c_kernel.cpp"
#include "lincom_kernel_eqA_f2c_kernel.cpp"
#include "lincom_kernel_eqB_f2c_kernel.cpp"
#include "lincom_kernel_eqC_f2c_kernel.cpp"
#include "lincom_kernel_eqD_f2c_kernel.cpp"
#include "lincom_kernel_eqE_f2c_kernel.cpp"
#include "lincom_kernel_eqF_f2c_kernel.cpp"
#include "fincom_kernel_main_f2c_kernel.cpp"
#include "adaptt_kernel_err_eval_f2c_kernel.cpp"
#include "adaptt_kernel_err_eval_MD_f2c_kernel.cpp"
#include "bcdt_kernel_xdir_f2c_kernel.cpp"
#include "bcdt_kernel_ydir_f2c_kernel.cpp"
#include "bcdt_kernel_zdir_f2c_kernel.cpp"
#include "bcdt_kernel_xdir_eqA_f2c_kernel.cpp"
#include "bcdt_kernel_ydir_eqA_f2c_kernel.cpp"
#include "bcdt_kernel_zdir_eqA_f2c_kernel.cpp"
#include "bcut_kernel_xdir_eqF_f2c_kernel.cpp"
#include "bcut_kernel_xdir_eqG_f2c_kernel.cpp"
#include "bcut_kernel_xdir_eqH_f2c_kernel.cpp"
#include "bcut_kernel_xdir_eqI_f2c_kernel.cpp"
#include "bcut_kernel_ydir_f2c_kernel.cpp"
#include "bcut_kernel_zdir_f2c_kernel.cpp"
#include "bcyt_kernel_xdir_eqA_f2c_kernel.cpp"
#include "bcyt_kernel_xdir_eqB_f2c_kernel.cpp"
#include "bcyt_kernel_xdir_eqC_f2c_kernel.cpp"
#include "bcyt_kernel_xdir_eqD_f2c_kernel.cpp"
#include "bcyt_kernel_ydir_f2c_kernel.cpp"
#include "bcyt_kernel_zdir_f2c_kernel.cpp"
#include "bounds_kernel_eqA_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqB_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqC_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqAA_xdir_f2c_kernel.cpp"
#include "copy_kernel_xxdir_f2c_kernel.cpp"
#include "bounds_kernel_eqD_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqF_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqH_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqP_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqQ_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqI_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqJ_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqR_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqS_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqE_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqG_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqK_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqT_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqL_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqU_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqV_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqAB_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqAC_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqAD_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqAE_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqAF_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqAG_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqAH_xdir_f2c_kernel.cpp"
#include "bounds_kernel_eqM_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqW_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqX_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqN_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqO_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqY_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqZ_xl_f2c_kernel.cpp"
#include "bounds_kernel_eqH_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqP_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqQ_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqI_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqJ_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqR_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqS_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqK_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqT_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqL_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqU_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqV_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqAF_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqM_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqW_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqX_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqN_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqO_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqY_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqZ_xr_f2c_kernel.cpp"
#include "bounds_kernel_eqA_ydir_f2c_kernel.cpp"
#include "bounds_kernel_eqB_ydir_f2c_kernel.cpp"
#include "bounds_kernel_eqC_ydir_f2c_kernel.cpp"
#include "bounds_kernel_eqAA_ydir_f2c_kernel.cpp"
#include "copy_kernel_yydir_f2c_kernel.cpp"
#include "bounds_kernel_eqD_ydir_f2c_kernel.cpp"
#include "bounds_kernel_eqF_ydir_f2c_kernel.cpp"
#include "bounds_kernel_eqH_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqP_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqQ_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqI_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqJ_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqR_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqS_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqE_ydir_f2c_kernel.cpp"
#include "bounds_kernel_eqG_ydir_f2c_kernel.cpp"
#include "bounds_kernel_eqK_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqT_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqL_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqU_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqV_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqAF_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqAG_ydir_f2c_kernel.cpp"
#include "bounds_kernel_eqAH_ydir_f2c_kernel.cpp"
#include "bounds_kernel_eqM_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqW_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqX_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqN_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqO_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqY_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqZ_yl_f2c_kernel.cpp"
#include "bounds_kernel_eqH_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqP_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqQ_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqI_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqJ_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqR_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqS_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqK_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqT_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqL_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqU_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqV_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqAF_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqM_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqW_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqX_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqN_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqO_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqY_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqZ_yr_f2c_kernel.cpp"
#include "bounds_kernel_eqA_zdir_f2c_kernel.cpp"
#include "bounds_kernel_eqB_zdir_f2c_kernel.cpp"
#include "bounds_kernel_eqC_zdir_f2c_kernel.cpp"
#include "bounds_kernel_eqAA_zdir_f2c_kernel.cpp"
#include "copy_kernel_zzdir_f2c_kernel.cpp"
#include "bounds_kernel_eqD_zdir_f2c_kernel.cpp"
#include "bounds_kernel_eqF_zdir_f2c_kernel.cpp"
#include "bounds_kernel_eqH_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqP_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqQ_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqI_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqJ_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqR_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqS_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqE_zdir_f2c_kernel.cpp"
#include "bounds_kernel_eqG_zdir_f2c_kernel.cpp"
#include "bounds_kernel_eqK_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqT_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqL_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqU_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqV_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqAF_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqAG_zdir_f2c_kernel.cpp"
#include "bounds_kernel_eqAH_zdir_f2c_kernel.cpp"
#include "bounds_kernel_eqM_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqW_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqX_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqN_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqO_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqY_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqZ_zl_f2c_kernel.cpp"
#include "bounds_kernel_eqH_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqP_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqQ_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqI_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqJ_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqR_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqS_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqK_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqT_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqL_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqU_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqV_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqAF_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqM_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqW_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqX_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqN_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqO_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqY_zr_f2c_kernel.cpp"
#include "bounds_kernel_eqZ_zr_f2c_kernel.cpp"
#include "boundt_kernel_eqE_xdir_f2c_kernel.cpp"
#include "bountt_kernel_eqA_xdir_f2c_kernel.cpp"
#include "bountt_kernel_eqB_xdir_f2c_kernel.cpp"
#include "bountt_kernel_eqF_xdir_f2c_kernel.cpp"
#include "bountt_kernel_eqD_f2c_kernel.cpp"
#include "bountt_kernel_eqC_xdir_f2c_kernel.cpp"
#include "bountt_kernel_eqE_xdir_f2c_kernel.cpp"
#include "bountt_kernel_eqG_xyz_f2c_kernel.cpp"
#include "boundt_kernel_eqG_xdir_f2c_kernel.cpp"
#include "bountt_kernel_eqH_xdir_f2c_kernel.cpp"
#include "boundt_kernel_eqE_ydir_f2c_kernel.cpp"
#include "bountt_kernel_eqA_ydir_f2c_kernel.cpp"
#include "bountt_kernel_eqB_ydir_f2c_kernel.cpp"
#include "bountt_kernel_eqF_ydir_f2c_kernel.cpp"
#include "bountt_kernel_eqC_ydir_f2c_kernel.cpp"
#include "bountt_kernel_eqE_ydir_f2c_kernel.cpp"
#include "boundt_kernel_eqG_ydir_f2c_kernel.cpp"
#include "bountt_kernel_eqH_ydir_f2c_kernel.cpp"
#include "boundt_kernel_eqE_zdir_f2c_kernel.cpp"
#include "bountt_kernel_eqA_zdir_f2c_kernel.cpp"
#include "bountt_kernel_eqB_zdir_f2c_kernel.cpp"
#include "bountt_kernel_eqF_zdir_f2c_kernel.cpp"
#include "bountt_kernel_eqC_zdir_f2c_kernel.cpp"
#include "bountt_kernel_eqE_zdir_f2c_kernel.cpp"
#include "boundt_kernel_eqG_zdir_f2c_kernel.cpp"
#include "bountt_kernel_eqH_zdir_f2c_kernel.cpp"
#include "boundt_kernel_eqA_xdir_f2c_kernel.cpp"
#include "boundt_kernel_eqB_xdir_f2c_kernel.cpp"
#include "boundt_kernel_eqF_xdir_f2c_kernel.cpp"
#include "boundt_kernel_eqC_xdir_f2c_kernel.cpp"
#include "boundt_kernel_eqD_xdir_f2c_kernel.cpp"
#include "boundt_kernel_eqH_xyz_f2c_kernel.cpp"
#include "boundt_kernel_eqA_ydir_f2c_kernel.cpp"
#include "boundt_kernel_eqB_ydir_f2c_kernel.cpp"
#include "boundt_kernel_eqF_ydir_f2c_kernel.cpp"
#include "boundt_kernel_eqC_ydir_f2c_kernel.cpp"
#include "boundt_kernel_eqD_ydir_f2c_kernel.cpp"
#include "boundt_kernel_eqA_zdir_f2c_kernel.cpp"
#include "boundt_kernel_eqB_zdir_f2c_kernel.cpp"
#include "boundt_kernel_eqF_zdir_f2c_kernel.cpp"
#include "boundt_kernel_eqC_zdir_f2c_kernel.cpp"
#include "boundt_kernel_eqD_zdir_f2c_kernel.cpp"
#include "radcal_kernel_meancoef_f2c_kernel.cpp"
#include "radcal_kernel_addspecies_f2c_kernel.cpp"
#include "radcal_kernel_addradiation_f2c_kernel.cpp"
#include "maths_kernel_eqJ_f2c_kernel.cpp"
#include "maths_kernel_eqAW_f2c_kernel.cpp"
#include "maths_kernel_eqAX_f2c_kernel.cpp"
#include "maths_kernel_eqAY_f2c_kernel.cpp"
#include "maths_kernel_eqAZ_f2c_kernel.cpp"
#include "maths_kernel_eqP_f2c_kernel.cpp"
#include "maths_kernel_eqBM_f2c_kernel.cpp"
#include "maths_kernel_eqBN_f2c_kernel.cpp"
#include "maths_kernel_eqB_f2c_kernel.cpp"
#include "maths_kernel_eqF_f2c_kernel.cpp"
#include "maths_kernel_eqE_f2c_kernel.cpp"
#include "maths_kernel_eqD_f2c_kernel.cpp"
#include "maths_kernel_eqK_f2c_kernel.cpp"
#include "maths_kernel_eqU_f2c_kernel.cpp"
#include "maths_kernel_eqBO_f2c_kernel.cpp"
#include "maths_kernel_eqBP_f2c_kernel.cpp"
#include "maths_kernel_eqAU_f2c_kernel.cpp"
#include "maths_kernel_eqV_fused_f2c_kernel.cpp"
#include "turbin_kernel_eqA_f2c_kernel.cpp"
#include "turbin_kernel_eqB_f2c_kernel.cpp"
#include "turbin_kernel_eqC_f2c_kernel.cpp"
#include "temper_kernel_eqA_f2c_kernel.cpp"
#include "set_zero_kernel_MD5_f2c_kernel.cpp"
#include "temper_kernel_eqB_f2c_kernel.cpp"
#include "temper_kernel_eqC_f2c_kernel.cpp"
#include "temper_kernel_eqD_f2c_kernel.cpp"
#include "temper_kernel_eqE_f2c_kernel.cpp"
#include "temper_kernel_eqF_f2c_kernel.cpp"
#include "tempin_kernel_main_f2c_kernel.cpp"
#include "copy_kernel_mdim_to_sdim_f2c_kernel.cpp"


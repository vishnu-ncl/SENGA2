// Auto-generated at 2026-04-28 18:44:28.306583 by ops-translator

// headers
#define OPS_3D
#define OPS_SOA
#define OPS_API 2
#include "ops_lib_core.h"

#ifdef OPS_MPI
#include "ops_mpi_core.h"
#include <limits>
#endif

#include "ops_f2c_prelude.h"

namespace f2c = ops::f2c;

#include "ops_hip_rt_support.h"
#include "ops_hip_reduction.h"

//  global constants
__constant__ int ncofmx;
__constant__ int ntinmx;
__constant__ int nspcmx;
__constant__ int nssmax;
__constant__ int nstpmx;
__constant__ int ndcfmx;
__constant__ int nvcfmx;
__constant__ int nccfmx;
__constant__ int nrkmax;
__constant__ int ncbcsz;
__constant__ int nbcprr;
__constant__ int nspimx;
__constant__ int ntbase;
__constant__ int nintmx;
__constant__ int nctmax;
__constant__ int nctmm1;
__constant__ int nrsmax;
__constant__ int nbcpri;
__constant__ int ncfrmx;
__constant__ double acoffx;
__constant__ double bcoffx;
__constant__ double ccoffx;
__constant__ double dcoffx;
__constant__ double ecoffx;
__constant__ double acof1x;
__constant__ double bcof1x;
__constant__ double ccof1x;
__constant__ double dcof1x;
__constant__ double acof2x;
__constant__ double bcof2x;
__constant__ double ccof2x;
__constant__ double dcof2x;
__constant__ double acof3x;
__constant__ double bcof3x;
__constant__ double acof4x;
__constant__ double bcof4x;
__constant__ double ccof4x;
__constant__ double acof5x;
__constant__ double bcof5x;
__constant__ double ccof5x;
__constant__ double dcof5x;
__constant__ double ovdelx;
__constant__ double acofsx;
__constant__ double bcofsx;
__constant__ double ccofsx;
__constant__ double dcofsx;
__constant__ double ecofsx;
__constant__ double acfs1x;
__constant__ double bcfs1x;
__constant__ double ccfs1x;
__constant__ double dcfs1x;
__constant__ double ecfs1x;
__constant__ double acfs2x;
__constant__ double bcfs2x;
__constant__ double ccfs2x;
__constant__ double dcfs2x;
__constant__ double ecfs2x;
__constant__ double acfs3x;
__constant__ double bcfs3x;
__constant__ double acfs4x;
__constant__ double bcfs4x;
__constant__ double ccfs4x;
__constant__ double acfs5x;
__constant__ double bcfs5x;
__constant__ double ccfs5x;
__constant__ double dcfs5x;
__constant__ double ovdlx2;
__constant__ double acoffy;
__constant__ double bcoffy;
__constant__ double ccoffy;
__constant__ double dcoffy;
__constant__ double ecoffy;
__constant__ double acof1y;
__constant__ double bcof1y;
__constant__ double ccof1y;
__constant__ double dcof1y;
__constant__ double acof2y;
__constant__ double bcof2y;
__constant__ double ccof2y;
__constant__ double dcof2y;
__constant__ double acof3y;
__constant__ double bcof3y;
__constant__ double acof4y;
__constant__ double bcof4y;
__constant__ double ccof4y;
__constant__ double acof5y;
__constant__ double bcof5y;
__constant__ double ccof5y;
__constant__ double dcof5y;
__constant__ double ovdely;
__constant__ double acofsy;
__constant__ double bcofsy;
__constant__ double ccofsy;
__constant__ double dcofsy;
__constant__ double ecofsy;
__constant__ double acfs1y;
__constant__ double bcfs1y;
__constant__ double ccfs1y;
__constant__ double dcfs1y;
__constant__ double ecfs1y;
__constant__ double acfs2y;
__constant__ double bcfs2y;
__constant__ double ccfs2y;
__constant__ double dcfs2y;
__constant__ double ecfs2y;
__constant__ double acfs3y;
__constant__ double bcfs3y;
__constant__ double acfs4y;
__constant__ double bcfs4y;
__constant__ double ccfs4y;
__constant__ double acfs5y;
__constant__ double bcfs5y;
__constant__ double ccfs5y;
__constant__ double dcfs5y;
__constant__ double ovdly2;
__constant__ double acoffz;
__constant__ double bcoffz;
__constant__ double ccoffz;
__constant__ double dcoffz;
__constant__ double ecoffz;
__constant__ double acof1z;
__constant__ double bcof1z;
__constant__ double ccof1z;
__constant__ double dcof1z;
__constant__ double acof2z;
__constant__ double bcof2z;
__constant__ double ccof2z;
__constant__ double dcof2z;
__constant__ double acof3z;
__constant__ double bcof3z;
__constant__ double acof4z;
__constant__ double bcof4z;
__constant__ double ccof4z;
__constant__ double acof5z;
__constant__ double bcof5z;
__constant__ double ccof5z;
__constant__ double dcof5z;
__constant__ double ovdelz;
__constant__ double acofsz;
__constant__ double bcofsz;
__constant__ double ccofsz;
__constant__ double dcofsz;
__constant__ double ecofsz;
__constant__ double acfs1z;
__constant__ double bcfs1z;
__constant__ double ccfs1z;
__constant__ double dcfs1z;
__constant__ double ecfs1z;
__constant__ double acfs2z;
__constant__ double bcfs2z;
__constant__ double ccfs2z;
__constant__ double dcfs2z;
__constant__ double ecfs2z;
__constant__ double acfs3z;
__constant__ double bcfs3z;
__constant__ double acfs4z;
__constant__ double bcfs4z;
__constant__ double ccfs4z;
__constant__ double acfs5z;
__constant__ double bcfs5z;
__constant__ double ccfs5z;
__constant__ double dcfs5z;
__constant__ double ovdlz2;
__constant__ double acofx1;
__constant__ double bcofx1;
__constant__ double acofy1;
__constant__ double bcofy1;
__constant__ double acofz1;
__constant__ double bcofz1;
__constant__ double acofxy;
__constant__ double bcofxy;
__constant__ double ccofxy;
__constant__ double dcofxy;
__constant__ double ecofxy;
__constant__ double acf1xy;
__constant__ double bcf1xy;
__constant__ double ccf1xy;
__constant__ double dcf1xy;
__constant__ double acf2xy;
__constant__ double bcf2xy;
__constant__ double ccf2xy;
__constant__ double dcf2xy;
__constant__ double acf3xy;
__constant__ double bcf3xy;
__constant__ double acf4xy;
__constant__ double bcf4xy;
__constant__ double ccf4xy;
__constant__ double acf5xy;
__constant__ double bcf5xy;
__constant__ double ccf5xy;
__constant__ double dcf5xy;
__constant__ double acc1xy;
__constant__ double bcc1xy;
__constant__ double ccc1xy;
__constant__ double dcc1xy;
__constant__ double acc2xy;
__constant__ double bcc2xy;
__constant__ double ccc2xy;
__constant__ double dcc2xy;
__constant__ double acofxz;
__constant__ double bcofxz;
__constant__ double ccofxz;
__constant__ double dcofxz;
__constant__ double ecofxz;
__constant__ double acf1xz;
__constant__ double bcf1xz;
__constant__ double ccf1xz;
__constant__ double dcf1xz;
__constant__ double acf2xz;
__constant__ double bcf2xz;
__constant__ double ccf2xz;
__constant__ double dcf2xz;
__constant__ double acf3xz;
__constant__ double bcf3xz;
__constant__ double acf4xz;
__constant__ double bcf4xz;
__constant__ double ccf4xz;
__constant__ double acf5xz;
__constant__ double bcf5xz;
__constant__ double ccf5xz;
__constant__ double dcf5xz;
__constant__ double acc1xz;
__constant__ double bcc1xz;
__constant__ double ccc1xz;
__constant__ double dcc1xz;
__constant__ double acc2xz;
__constant__ double bcc2xz;
__constant__ double ccc2xz;
__constant__ double dcc2xz;
__constant__ double acofyz;
__constant__ double bcofyz;
__constant__ double ccofyz;
__constant__ double dcofyz;
__constant__ double ecofyz;
__constant__ double acf1yz;
__constant__ double bcf1yz;
__constant__ double ccf1yz;
__constant__ double dcf1yz;
__constant__ double acf2yz;
__constant__ double bcf2yz;
__constant__ double ccf2yz;
__constant__ double dcf2yz;
__constant__ double acf3yz;
__constant__ double bcf3yz;
__constant__ double acf4yz;
__constant__ double bcf4yz;
__constant__ double ccf4yz;
__constant__ double acf5yz;
__constant__ double bcf5yz;
__constant__ double ccf5yz;
__constant__ double dcf5yz;
__constant__ double acc1yz;
__constant__ double bcc1yz;
__constant__ double ccc1yz;
__constant__ double dcc1yz;
__constant__ double acc2yz;
__constant__ double bcc2yz;
__constant__ double ccc2yz;
__constant__ double dcc2yz;
__constant__ double foursb;
__constant__ double trfrth;
__constant__ double rlamda;
__constant__ double alamda;

void ops_init_backend(){}

extern "C"
void ops_decl_const_f2c(char const *name, int dim, int size, char *dat) {

    OPS_instance *instance = OPS_instance::getOPSInstance();
    ops_execute(instance);

    if(!strcmp(name, "ncofmx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ncofmx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ntinmx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ntinmx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nspcmx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nspcmx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nssmax")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nssmax), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nstpmx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nstpmx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ndcfmx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ndcfmx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nvcfmx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nvcfmx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nccfmx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nccfmx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nrkmax")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nrkmax), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ncbcsz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ncbcsz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nbcprr")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nbcprr), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nspimx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nspimx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ntbase")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ntbase), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nintmx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nintmx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nctmax")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nctmax), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nctmm1")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nctmm1), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nrsmax")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nrsmax), dat, dim*size));
    } 
    else
    if(!strcmp(name, "nbcpri")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(nbcpri), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ncfrmx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ncfrmx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acoffx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acoffx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcoffx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcoffx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccoffx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccoffx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcoffx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcoffx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecoffx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecoffx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof1x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof1x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof1x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof1x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof1x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof1x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof1x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcof1x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof2x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof2x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof2x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof2x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof2x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof2x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof2x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcof2x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof3x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof3x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof3x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof3x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof4x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof4x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof4x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof4x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof4x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof4x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof5x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof5x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof5x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof5x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof5x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof5x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof5x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcof5x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdelx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ovdelx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofsx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acofsx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofsx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcofsx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofsx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccofsx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofsx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcofsx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofsx")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecofsx), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs1x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs1x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs1x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs1x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs1x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs1x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs1x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcfs1x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs1x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecfs1x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs2x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs2x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs2x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs2x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs2x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs2x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs2x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcfs2x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs2x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecfs2x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs3x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs3x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs3x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs3x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs4x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs4x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs4x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs4x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs4x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs4x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs5x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs5x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs5x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs5x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs5x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs5x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs5x")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcfs5x), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdlx2")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ovdlx2), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acoffy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acoffy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcoffy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcoffy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccoffy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccoffy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcoffy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcoffy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecoffy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecoffy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof1y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof1y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof1y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof1y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof1y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof1y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof1y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcof1y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof2y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof2y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof2y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof2y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof2y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof2y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof2y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcof2y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof3y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof3y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof3y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof3y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof4y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof4y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof4y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof4y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof4y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof4y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof5y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof5y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof5y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof5y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof5y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof5y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof5y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcof5y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdely")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ovdely), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofsy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acofsy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofsy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcofsy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofsy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccofsy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofsy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcofsy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofsy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecofsy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs1y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs1y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs1y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs1y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs1y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs1y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs1y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcfs1y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs1y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecfs1y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs2y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs2y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs2y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs2y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs2y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs2y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs2y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcfs2y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs2y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecfs2y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs3y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs3y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs3y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs3y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs4y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs4y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs4y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs4y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs4y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs4y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs5y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs5y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs5y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs5y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs5y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs5y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs5y")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcfs5y), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdly2")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ovdly2), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acoffz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acoffz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcoffz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcoffz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccoffz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccoffz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcoffz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcoffz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecoffz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecoffz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof1z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof1z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof1z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof1z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof1z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof1z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof1z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcof1z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof2z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof2z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof2z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof2z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof2z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof2z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof2z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcof2z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof3z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof3z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof3z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof3z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof4z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof4z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof4z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof4z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof4z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof4z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof5z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acof5z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof5z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcof5z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof5z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccof5z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof5z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcof5z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdelz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ovdelz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofsz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acofsz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofsz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcofsz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofsz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccofsz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofsz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcofsz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofsz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecofsz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs1z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs1z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs1z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs1z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs1z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs1z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs1z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcfs1z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs1z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecfs1z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs2z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs2z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs2z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs2z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs2z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs2z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs2z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcfs2z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs2z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecfs2z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs3z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs3z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs3z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs3z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs4z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs4z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs4z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs4z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs4z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs4z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs5z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acfs5z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs5z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcfs5z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs5z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccfs5z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs5z")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcfs5z), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdlz2")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ovdlz2), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofx1")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acofx1), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofx1")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcofx1), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofy1")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acofy1), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofy1")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcofy1), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofz1")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acofz1), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofz1")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcofz1), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofxy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acofxy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofxy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcofxy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofxy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccofxy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofxy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcofxy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofxy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecofxy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf1xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf1xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf1xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf1xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf1xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf1xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf1xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcf1xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf2xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf2xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf2xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf2xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf2xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf2xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf2xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcf2xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf3xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf3xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf3xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf3xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf4xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf4xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf4xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf4xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf4xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf4xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf5xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf5xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf5xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf5xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf5xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf5xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf5xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcf5xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc1xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acc1xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc1xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcc1xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc1xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccc1xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc1xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcc1xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc2xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acc2xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc2xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcc2xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc2xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccc2xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc2xy")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcc2xy), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofxz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acofxz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofxz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcofxz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofxz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccofxz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofxz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcofxz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofxz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecofxz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf1xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf1xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf1xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf1xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf1xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf1xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf1xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcf1xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf2xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf2xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf2xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf2xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf2xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf2xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf2xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcf2xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf3xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf3xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf3xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf3xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf4xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf4xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf4xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf4xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf4xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf4xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf5xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf5xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf5xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf5xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf5xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf5xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf5xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcf5xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc1xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acc1xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc1xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcc1xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc1xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccc1xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc1xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcc1xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc2xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acc2xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc2xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcc2xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc2xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccc2xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc2xz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcc2xz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofyz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acofyz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofyz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcofyz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofyz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccofyz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofyz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcofyz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofyz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ecofyz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf1yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf1yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf1yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf1yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf1yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf1yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf1yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcf1yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf2yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf2yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf2yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf2yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf2yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf2yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf2yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcf2yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf3yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf3yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf3yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf3yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf4yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf4yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf4yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf4yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf4yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf4yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf5yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acf5yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf5yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcf5yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf5yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccf5yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf5yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcf5yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc1yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acc1yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc1yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcc1yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc1yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccc1yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc1yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcc1yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc2yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(acc2yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc2yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(bcc2yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc2yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(ccc2yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc2yz")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(dcc2yz), dat, dim*size));
    } 
    else
    if(!strcmp(name, "foursb")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(foursb), dat, dim*size));
    } 
    else
    if(!strcmp(name, "trfrth")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(trfrth), dat, dim*size));
    } 
    else
    if(!strcmp(name, "rlamda")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(rlamda), dat, dim*size));
    } 
    else
    if(!strcmp(name, "alamda")) {
        hipSafeCall(instance->ostream(), hipMemcpyToSymbol(HIP_SYMBOL(alamda), dat, dim*size));
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


// Auto-generated at 2026-04-28 18:44:05.630978 by ops-translator

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

#include "ops_cuda_rt_support.h"
#include "ops_cuda_reduction.h"

#include <cuComplex.h>

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
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ncofmx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ntinmx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ntinmx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nspcmx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nspcmx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nssmax")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nssmax, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nstpmx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nstpmx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ndcfmx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ndcfmx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nvcfmx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nvcfmx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nccfmx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nccfmx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nrkmax")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nrkmax, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ncbcsz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ncbcsz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nbcprr")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nbcprr, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nspimx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nspimx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ntbase")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ntbase, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nintmx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nintmx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nctmax")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nctmax, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nctmm1")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nctmm1, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nrsmax")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nrsmax, dat, dim*size));
    } 
    else
    if(!strcmp(name, "nbcpri")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(nbcpri, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ncfrmx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ncfrmx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acoffx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acoffx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcoffx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcoffx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccoffx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccoffx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcoffx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcoffx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecoffx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecoffx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof1x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof1x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof1x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof1x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof1x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof1x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof1x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcof1x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof2x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof2x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof2x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof2x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof2x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof2x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof2x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcof2x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof3x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof3x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof3x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof3x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof4x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof4x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof4x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof4x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof4x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof4x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof5x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof5x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof5x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof5x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof5x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof5x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof5x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcof5x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdelx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ovdelx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofsx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acofsx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofsx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcofsx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofsx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccofsx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofsx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcofsx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofsx")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecofsx, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs1x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs1x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs1x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs1x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs1x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs1x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs1x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcfs1x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs1x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecfs1x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs2x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs2x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs2x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs2x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs2x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs2x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs2x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcfs2x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs2x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecfs2x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs3x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs3x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs3x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs3x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs4x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs4x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs4x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs4x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs4x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs4x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs5x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs5x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs5x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs5x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs5x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs5x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs5x")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcfs5x, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdlx2")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ovdlx2, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acoffy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acoffy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcoffy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcoffy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccoffy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccoffy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcoffy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcoffy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecoffy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecoffy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof1y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof1y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof1y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof1y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof1y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof1y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof1y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcof1y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof2y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof2y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof2y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof2y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof2y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof2y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof2y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcof2y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof3y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof3y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof3y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof3y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof4y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof4y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof4y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof4y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof4y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof4y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof5y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof5y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof5y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof5y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof5y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof5y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof5y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcof5y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdely")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ovdely, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofsy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acofsy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofsy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcofsy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofsy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccofsy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofsy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcofsy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofsy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecofsy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs1y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs1y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs1y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs1y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs1y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs1y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs1y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcfs1y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs1y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecfs1y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs2y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs2y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs2y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs2y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs2y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs2y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs2y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcfs2y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs2y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecfs2y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs3y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs3y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs3y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs3y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs4y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs4y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs4y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs4y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs4y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs4y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs5y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs5y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs5y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs5y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs5y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs5y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs5y")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcfs5y, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdly2")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ovdly2, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acoffz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acoffz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcoffz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcoffz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccoffz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccoffz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcoffz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcoffz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecoffz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecoffz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof1z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof1z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof1z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof1z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof1z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof1z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof1z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcof1z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof2z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof2z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof2z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof2z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof2z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof2z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof2z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcof2z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof3z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof3z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof3z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof3z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof4z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof4z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof4z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof4z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof4z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof4z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acof5z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acof5z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcof5z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcof5z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccof5z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccof5z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcof5z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcof5z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdelz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ovdelz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofsz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acofsz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofsz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcofsz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofsz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccofsz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofsz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcofsz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofsz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecofsz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs1z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs1z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs1z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs1z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs1z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs1z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs1z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcfs1z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs1z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecfs1z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs2z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs2z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs2z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs2z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs2z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs2z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs2z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcfs2z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecfs2z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecfs2z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs3z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs3z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs3z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs3z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs4z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs4z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs4z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs4z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs4z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs4z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acfs5z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acfs5z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcfs5z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcfs5z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccfs5z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccfs5z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcfs5z")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcfs5z, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ovdlz2")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ovdlz2, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofx1")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acofx1, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofx1")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcofx1, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofy1")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acofy1, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofy1")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcofy1, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofz1")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acofz1, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofz1")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcofz1, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofxy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acofxy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofxy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcofxy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofxy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccofxy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofxy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcofxy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofxy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecofxy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf1xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf1xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf1xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf1xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf1xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf1xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf1xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcf1xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf2xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf2xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf2xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf2xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf2xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf2xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf2xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcf2xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf3xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf3xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf3xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf3xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf4xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf4xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf4xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf4xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf4xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf4xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf5xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf5xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf5xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf5xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf5xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf5xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf5xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcf5xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc1xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acc1xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc1xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcc1xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc1xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccc1xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc1xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcc1xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc2xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acc2xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc2xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcc2xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc2xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccc2xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc2xy")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcc2xy, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofxz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acofxz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofxz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcofxz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofxz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccofxz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofxz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcofxz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofxz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecofxz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf1xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf1xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf1xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf1xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf1xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf1xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf1xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcf1xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf2xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf2xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf2xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf2xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf2xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf2xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf2xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcf2xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf3xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf3xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf3xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf3xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf4xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf4xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf4xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf4xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf4xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf4xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf5xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf5xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf5xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf5xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf5xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf5xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf5xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcf5xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc1xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acc1xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc1xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcc1xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc1xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccc1xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc1xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcc1xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc2xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acc2xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc2xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcc2xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc2xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccc2xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc2xz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcc2xz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acofyz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acofyz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcofyz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcofyz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccofyz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccofyz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcofyz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcofyz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ecofyz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ecofyz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf1yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf1yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf1yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf1yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf1yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf1yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf1yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcf1yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf2yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf2yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf2yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf2yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf2yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf2yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf2yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcf2yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf3yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf3yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf3yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf3yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf4yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf4yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf4yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf4yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf4yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf4yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acf5yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acf5yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcf5yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcf5yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccf5yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccf5yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcf5yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcf5yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc1yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acc1yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc1yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcc1yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc1yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccc1yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc1yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcc1yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "acc2yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(acc2yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "bcc2yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(bcc2yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "ccc2yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(ccc2yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "dcc2yz")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(dcc2yz, dat, dim*size));
    } 
    else
    if(!strcmp(name, "foursb")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(foursb, dat, dim*size));
    } 
    else
    if(!strcmp(name, "trfrth")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(trfrth, dat, dim*size));
    } 
    else
    if(!strcmp(name, "rlamda")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(rlamda, dat, dim*size));
    } 
    else
    if(!strcmp(name, "alamda")) {
        cutilSafeCall(instance->ostream(), cudaMemcpyToSymbol(alamda, dat, dim*size));
    } 
    else
    {
        throw OPSException(OPS_RUNTIME_ERROR, "error: unknown const name");
    }
}

// user kernel files
#include "set_zero_kernel_f2c_kernel.cu"
#include "set_zero_kernel_xdir_f2c_kernel.cu"
#include "set_zero_kernel_ydir_f2c_kernel.cu"
#include "set_zero_kernel_zdir_f2c_kernel.cu"
#include "set_zero_kernel_int_f2c_kernel.cu"
#include "dfbydx_kernel_main_f2c_kernel.cu"
#include "d2fdx2_kernel_main_f2c_kernel.cu"
#include "dfbydy_kernel_null_f2c_kernel.cu"
#include "dfbydy_kernel_main_f2c_kernel.cu"
#include "dfbydz_kernel_null_f2c_kernel.cu"
#include "dfbydz_kernel_main_f2c_kernel.cu"
#include "d2fdy2_kernel_null_f2c_kernel.cu"
#include "d2fdy2_kernel_main_f2c_kernel.cu"
#include "d2fdz2_kernel_null_f2c_kernel.cu"
#include "d2fdz2_kernel_main_f2c_kernel.cu"
#include "d2fdxy_kernel_null_f2c_kernel.cu"
#include "d2fdxy_kernel_interior_f2c_kernel.cu"
#include "d2fdxy_kernel_eqA_f2c_kernel.cu"
#include "d2fdxy_kernel_eqB_f2c_kernel.cu"
#include "d2fdxy_kernel_eqC_f2c_kernel.cu"
#include "d2fdxy_kernel_eqD_f2c_kernel.cu"
#include "d2fdxy_kernel_eqE_f2c_kernel.cu"
#include "d2fdxy_kernel_eqF_f2c_kernel.cu"
#include "d2fdxy_kernel_eqG_f2c_kernel.cu"
#include "d2fdxy_kernel_eqH_f2c_kernel.cu"
#include "d2fdxy_kernel_eqI_f2c_kernel.cu"
#include "d2fdxy_kernel_eqJ_f2c_kernel.cu"
#include "d2fdxy_kernel_eqK_f2c_kernel.cu"
#include "d2fdxy_kernel_eqL_f2c_kernel.cu"
#include "d2fdxy_kernel_eqM_f2c_kernel.cu"
#include "d2fdxy_kernel_eqN_f2c_kernel.cu"
#include "d2fdxy_kernel_eqO_f2c_kernel.cu"
#include "d2fdxy_kernel_eqP_f2c_kernel.cu"
#include "d2fdxy_kernel_eqQ_f2c_kernel.cu"
#include "d2fdxy_kernel_eqR_f2c_kernel.cu"
#include "d2fdxy_kernel_eqS_f2c_kernel.cu"
#include "d2fdxy_kernel_eqT_f2c_kernel.cu"
#include "d2fdxy_kernel_eqU_f2c_kernel.cu"
#include "d2fdxy_kernel_eqV_f2c_kernel.cu"
#include "d2fdxy_kernel_eqW_f2c_kernel.cu"
#include "d2fdxy_kernel_eqX_f2c_kernel.cu"
#include "d2fdxy_kernel_eqY_f2c_kernel.cu"
#include "d2fdxy_kernel_eqZ_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAA_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAB_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAC_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAD_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAE_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAF_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAG_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAH_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAI_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAJ_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAK_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAL_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAM_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAN_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAO_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAP_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAQ_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAR_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAS_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAT_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAU_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAV_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAW_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAX_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAY_f2c_kernel.cu"
#include "d2fdxy_kernel_eqAZ_f2c_kernel.cu"
#include "d2fdxy_kernel_eqBA_f2c_kernel.cu"
#include "d2fdxy_kernel_eqBB_f2c_kernel.cu"
#include "d2fdxy_kernel_eqBC_f2c_kernel.cu"
#include "d2fdxy_kernel_eqBD_f2c_kernel.cu"
#include "d2fdxy_kernel_scaling_f2c_kernel.cu"
#include "d2fdxz_kernel_null_f2c_kernel.cu"
#include "d2fdxz_kernel_interior_f2c_kernel.cu"
#include "d2fdxz_kernel_eqA_f2c_kernel.cu"
#include "d2fdxz_kernel_eqB_f2c_kernel.cu"
#include "d2fdxz_kernel_eqC_f2c_kernel.cu"
#include "d2fdxz_kernel_eqD_f2c_kernel.cu"
#include "d2fdxz_kernel_eqE_f2c_kernel.cu"
#include "d2fdxz_kernel_eqF_f2c_kernel.cu"
#include "d2fdxz_kernel_eqG_f2c_kernel.cu"
#include "d2fdxz_kernel_eqH_f2c_kernel.cu"
#include "d2fdxz_kernel_eqI_f2c_kernel.cu"
#include "d2fdxz_kernel_eqJ_f2c_kernel.cu"
#include "d2fdxz_kernel_eqK_f2c_kernel.cu"
#include "d2fdxz_kernel_eqL_f2c_kernel.cu"
#include "d2fdxz_kernel_eqM_f2c_kernel.cu"
#include "d2fdxz_kernel_eqN_f2c_kernel.cu"
#include "d2fdxz_kernel_eqO_f2c_kernel.cu"
#include "d2fdxz_kernel_eqP_f2c_kernel.cu"
#include "d2fdxz_kernel_eqQ_f2c_kernel.cu"
#include "d2fdxz_kernel_eqR_f2c_kernel.cu"
#include "d2fdxz_kernel_eqS_f2c_kernel.cu"
#include "d2fdxz_kernel_eqT_f2c_kernel.cu"
#include "d2fdxz_kernel_eqU_f2c_kernel.cu"
#include "d2fdxz_kernel_eqV_f2c_kernel.cu"
#include "d2fdxz_kernel_eqW_f2c_kernel.cu"
#include "d2fdxz_kernel_eqX_f2c_kernel.cu"
#include "d2fdxz_kernel_eqY_f2c_kernel.cu"
#include "d2fdxz_kernel_eqZ_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAA_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAB_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAC_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAD_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAE_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAF_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAG_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAH_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAI_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAJ_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAK_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAL_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAM_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAN_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAO_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAP_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAQ_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAR_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAS_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAT_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAU_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAV_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAW_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAX_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAY_f2c_kernel.cu"
#include "d2fdxz_kernel_eqAZ_f2c_kernel.cu"
#include "d2fdxz_kernel_eqBA_f2c_kernel.cu"
#include "d2fdxz_kernel_eqBB_f2c_kernel.cu"
#include "d2fdxz_kernel_eqBC_f2c_kernel.cu"
#include "d2fdxz_kernel_eqBD_f2c_kernel.cu"
#include "d2fdxz_kernel_scaling_f2c_kernel.cu"
#include "d2fdyz_kernel_null_f2c_kernel.cu"
#include "d2fdyz_kernel_interior_f2c_kernel.cu"
#include "d2fdyz_kernel_eqA_f2c_kernel.cu"
#include "d2fdyz_kernel_eqB_f2c_kernel.cu"
#include "d2fdyz_kernel_eqC_f2c_kernel.cu"
#include "d2fdyz_kernel_eqD_f2c_kernel.cu"
#include "d2fdyz_kernel_eqE_f2c_kernel.cu"
#include "d2fdyz_kernel_eqF_f2c_kernel.cu"
#include "d2fdyz_kernel_eqG_f2c_kernel.cu"
#include "d2fdyz_kernel_eqH_f2c_kernel.cu"
#include "d2fdyz_kernel_eqI_f2c_kernel.cu"
#include "d2fdyz_kernel_eqJ_f2c_kernel.cu"
#include "d2fdyz_kernel_eqK_f2c_kernel.cu"
#include "d2fdyz_kernel_eqL_f2c_kernel.cu"
#include "d2fdyz_kernel_eqM_f2c_kernel.cu"
#include "d2fdyz_kernel_eqN_f2c_kernel.cu"
#include "d2fdyz_kernel_eqO_f2c_kernel.cu"
#include "d2fdyz_kernel_eqP_f2c_kernel.cu"
#include "d2fdyz_kernel_eqQ_f2c_kernel.cu"
#include "d2fdyz_kernel_eqR_f2c_kernel.cu"
#include "d2fdyz_kernel_eqS_f2c_kernel.cu"
#include "d2fdyz_kernel_eqT_f2c_kernel.cu"
#include "d2fdyz_kernel_eqU_f2c_kernel.cu"
#include "d2fdyz_kernel_eqV_f2c_kernel.cu"
#include "d2fdyz_kernel_eqW_f2c_kernel.cu"
#include "d2fdyz_kernel_eqX_f2c_kernel.cu"
#include "d2fdyz_kernel_eqY_f2c_kernel.cu"
#include "d2fdyz_kernel_eqZ_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAA_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAB_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAC_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAD_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAE_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAF_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAG_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAH_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAI_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAJ_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAK_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAL_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAM_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAN_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAO_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAP_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAQ_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAR_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAS_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAT_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAU_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAV_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAW_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAX_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAY_f2c_kernel.cu"
#include "d2fdyz_kernel_eqAZ_f2c_kernel.cu"
#include "d2fdyz_kernel_eqBA_f2c_kernel.cu"
#include "d2fdyz_kernel_eqBB_f2c_kernel.cu"
#include "d2fdyz_kernel_eqBC_f2c_kernel.cu"
#include "d2fdyz_kernel_eqBD_f2c_kernel.cu"
#include "d2fdyz_kernel_scaling_f2c_kernel.cu"
#include "boundary_kernel_CPandGAS_xdir_f2c_kernel.cu"
#include "boundary_kernel_CPandGAS_ydir_f2c_kernel.cu"
#include "boundary_kernel_CPandGAS_zdir_f2c_kernel.cu"
#include "maths_kernel_eqX_f2c_kernel.cu"
#include "maths_kernel_eqBR_xdir_f2c_kernel.cu"
#include "maths_kernel_eqBR_ydir_f2c_kernel.cu"
#include "maths_kernel_eqBR_zdir_f2c_kernel.cu"
#include "maths_kernel_eqT_f2c_kernel.cu"
#include "boundary_kernel_internalenergy_xdir_f2c_kernel.cu"
#include "boundary_kernel_internalenergy_ydir_f2c_kernel.cu"
#include "boundary_kernel_internalenergy_zdir_f2c_kernel.cu"
#include "maths_kernel_eqW_f2c_kernel.cu"
#include "maths_kernel_eqZ_f2c_kernel.cu"
#include "maths_kernel_eqAS_f2c_kernel.cu"
#include "boundary_kernel_temperature_xdir_f2c_kernel.cu"
#include "boundary_kernel_temperature_ydir_f2c_kernel.cu"
#include "boundary_kernel_temperature_zdir_f2c_kernel.cu"
#include "maths_kernel_eqAV_f2c_kernel.cu"
#include "maths_kernel_eqBC_f2c_kernel.cu"
#include "maths_kernel_eqBD_f2c_kernel.cu"
#include "maths_kernel_eqBE_f2c_kernel.cu"
#include "maths_kernel_eqAH_f2c_kernel.cu"
#include "hf_kernel_eqA_f2c_kernel.cu"
#include "hf_kernel_eqB_f2c_kernel.cu"
#include "hf_kernel_eqC_f2c_kernel.cu"
#include "hf_kernel_eqD_f2c_kernel.cu"
#include "hf_kernel_eqE_f2c_kernel.cu"
#include "hf_kernel_eqF_f2c_kernel.cu"
#include "maths_kernel_eqAF_f2c_kernel.cu"
#include "copy_kernel_f2c_kernel.cu"
#include "copy_kernel_xdir_f2c_kernel.cu"
#include "copy_kernel_ydir_f2c_kernel.cu"
#include "copy_kernel_zdir_f2c_kernel.cu"
#include "maths_kernel_eqA_f2c_kernel.cu"
#include "maths_kernel_eqAP_f2c_kernel.cu"
#include "maths_kernel_eqAQ_f2c_kernel.cu"
#include "boundary_kernel_mass_xdir_f2c_kernel.cu"
#include "boundary_kernel_mass_ydir_f2c_kernel.cu"
#include "boundary_kernel_mass_zdir_f2c_kernel.cu"
#include "maths_kernel_eqAT_f2c_kernel.cu"
#include "maths_kernel_eqH_f2c_kernel.cu"
#include "copy_kernel_sdim_to_mdim_f2c_kernel.cu"
#include "maths_kernel_eqBFG_f2c_kernel.cu"
#include "maths_kernel_eqBH_f2c_kernel.cu"
#include "maths_kernel_eqBS_f2c_kernel.cu"
#include "maths_kernel_eqAA_fused_f2c_kernel.cu"
#include "maths_kernel_eqBL_f2c_kernel.cu"
#include "boundary_kernel_speciesH_xdir_f2c_kernel.cu"
#include "boundary_kernel_speciesH_ydir_f2c_kernel.cu"
#include "boundary_kernel_speciesH_zdir_f2c_kernel.cu"
#include "maths_kernel_eqL_f2c_kernel.cu"
#include "maths_kernel_eqBA_f2c_kernel.cu"
#include "maths_kernel_eqAI_f2c_kernel.cu"
#include "hf_kernel_eqS_f2c_kernel.cu"
#include "hf_kernel_eqT_f2c_kernel.cu"
#include "hf_kernel_eqU_f2c_kernel.cu"
#include "hf_kernel_eqV_f2c_kernel.cu"
#include "hf_kernel_eqW_f2c_kernel.cu"
#include "hf_kernel_eqX_f2c_kernel.cu"
#include "maths_kernel_eqBB_f2c_kernel.cu"
#include "maths_kernel_eqAG_f2c_kernel.cu"
#include "maths_kernel_eqO_f2c_kernel.cu"
#include "maths_kernel_eqY_f2c_kernel.cu"
#include "maths_kernel_eqM_f2c_kernel.cu"
#include "hf_kernel_eqM_f2c_kernel.cu"
#include "hf_kernel_eqN_f2c_kernel.cu"
#include "hf_kernel_eqO_f2c_kernel.cu"
#include "hf_kernel_eqP_f2c_kernel.cu"
#include "hf_kernel_eqQ_f2c_kernel.cu"
#include "hf_kernel_eqR_f2c_kernel.cu"
#include "hf_kernel_eqG_f2c_kernel.cu"
#include "hf_kernel_eqH_f2c_kernel.cu"
#include "hf_kernel_eqI_f2c_kernel.cu"
#include "hf_kernel_eqJ_f2c_kernel.cu"
#include "hf_kernel_eqK_f2c_kernel.cu"
#include "hf_kernel_eqL_f2c_kernel.cu"
#include "maths_kernel_eqAM_f2c_kernel.cu"
#include "maths_kernel_eqBIJK_f2c_kernel.cu"
#include "maths_kernel_eqAN_f2c_kernel.cu"
#include "boundary_kernel_density_xdir_f2c_kernel.cu"
#include "boundary_kernel_density_ydir_f2c_kernel.cu"
#include "boundary_kernel_density_zdir_f2c_kernel.cu"
#include "maths_kernel_eqU_fused_f2c_kernel.cu"
#include "boundary_kernel_velcomp_xdir_f2c_kernel.cu"
#include "boundary_kernel_velcomp_ydir_f2c_kernel.cu"
#include "boundary_kernel_velcomp_zdir_f2c_kernel.cu"
#include "boundary_kernel_eqA_xdir_f2c_kernel.cu"
#include "boundary_kernel_eqA_ydir_f2c_kernel.cu"
#include "boundary_kernel_eqA_zdir_f2c_kernel.cu"
#include "maths_kernel_eqQ_f2c_kernel.cu"
#include "boundary_kernel_velderiv_xdir_f2c_kernel.cu"
#include "boundary_kernel_velderiv_ydir_f2c_kernel.cu"
#include "boundary_kernel_velderiv_zdir_f2c_kernel.cu"
#include "maths_kernel_eqAB_f2c_kernel.cu"
#include "maths_kernel_eqAO_f2c_kernel.cu"
#include "boundary_kernel_pressure_xdir_f2c_kernel.cu"
#include "boundary_kernel_pressure_ydir_f2c_kernel.cu"
#include "boundary_kernel_pressure_zdir_f2c_kernel.cu"
#include "maths_kernel_eqS_f2c_kernel.cu"
#include "maths_kernel_eqAL_f2c_kernel.cu"
#include "maths_kernel_eqAK_f2c_kernel.cu"
#include "maths_kernel_eqG_f2c_kernel.cu"
#include "maths_kernel_eqR_f2c_kernel.cu"
#include "maths_kernel_eqN_f2c_kernel.cu"
#include "maths_kernel_eqAD_f2c_kernel.cu"
#include "maths_kernel_eqV_f2c_kernel.cu"
#include "maths_kernel_eqAA_f2c_kernel.cu"
#include "maths_kernel_eqAR_f2c_kernel.cu"
#include "maths_kernel_eqI_f2c_kernel.cu"
#include "maths_kernel_eqAJ_f2c_kernel.cu"
#include "maths_kernel_eqAE_f2c_kernel.cu"
#include "maths_kernel_eqAC_f2c_kernel.cu"
#include "maths_kernel_eqtau_f2c_kernel.cu"
#include "maths_kernel_eqC_f2c_kernel.cu"
#include "lincom_kernel_main_f2c_kernel.cu"
#include "lincom_kernel_eqA_f2c_kernel.cu"
#include "lincom_kernel_eqB_f2c_kernel.cu"
#include "lincom_kernel_eqC_f2c_kernel.cu"
#include "lincom_kernel_eqD_f2c_kernel.cu"
#include "lincom_kernel_eqE_f2c_kernel.cu"
#include "lincom_kernel_eqF_f2c_kernel.cu"
#include "fincom_kernel_main_f2c_kernel.cu"
#include "adaptt_kernel_err_eval_f2c_kernel.cu"
#include "adaptt_kernel_err_eval_MD_f2c_kernel.cu"
#include "bcdt_kernel_xdir_f2c_kernel.cu"
#include "bcdt_kernel_ydir_f2c_kernel.cu"
#include "bcdt_kernel_zdir_f2c_kernel.cu"
#include "bcdt_kernel_xdir_eqA_f2c_kernel.cu"
#include "bcdt_kernel_ydir_eqA_f2c_kernel.cu"
#include "bcdt_kernel_zdir_eqA_f2c_kernel.cu"
#include "bcut_kernel_xdir_eqF_f2c_kernel.cu"
#include "bcut_kernel_xdir_eqG_f2c_kernel.cu"
#include "bcut_kernel_xdir_eqH_f2c_kernel.cu"
#include "bcut_kernel_xdir_eqI_f2c_kernel.cu"
#include "bcut_kernel_ydir_f2c_kernel.cu"
#include "bcut_kernel_zdir_f2c_kernel.cu"
#include "bcyt_kernel_xdir_eqA_f2c_kernel.cu"
#include "bcyt_kernel_xdir_eqB_f2c_kernel.cu"
#include "bcyt_kernel_xdir_eqC_f2c_kernel.cu"
#include "bcyt_kernel_xdir_eqD_f2c_kernel.cu"
#include "bcyt_kernel_ydir_f2c_kernel.cu"
#include "bcyt_kernel_zdir_f2c_kernel.cu"
#include "bounds_kernel_eqA_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqB_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqC_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqAA_xdir_f2c_kernel.cu"
#include "copy_kernel_xxdir_f2c_kernel.cu"
#include "bounds_kernel_eqD_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqF_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqH_xl_f2c_kernel.cu"
#include "bounds_kernel_eqP_xl_f2c_kernel.cu"
#include "bounds_kernel_eqQ_xl_f2c_kernel.cu"
#include "bounds_kernel_eqI_xl_f2c_kernel.cu"
#include "bounds_kernel_eqJ_xl_f2c_kernel.cu"
#include "bounds_kernel_eqR_xl_f2c_kernel.cu"
#include "bounds_kernel_eqS_xl_f2c_kernel.cu"
#include "bounds_kernel_eqE_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqG_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqK_xl_f2c_kernel.cu"
#include "bounds_kernel_eqT_xl_f2c_kernel.cu"
#include "bounds_kernel_eqL_xl_f2c_kernel.cu"
#include "bounds_kernel_eqU_xl_f2c_kernel.cu"
#include "bounds_kernel_eqV_xl_f2c_kernel.cu"
#include "bounds_kernel_eqAB_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqAC_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqAD_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqAE_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqAF_xl_f2c_kernel.cu"
#include "bounds_kernel_eqAG_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqAH_xdir_f2c_kernel.cu"
#include "bounds_kernel_eqM_xl_f2c_kernel.cu"
#include "bounds_kernel_eqW_xl_f2c_kernel.cu"
#include "bounds_kernel_eqX_xl_f2c_kernel.cu"
#include "bounds_kernel_eqN_xl_f2c_kernel.cu"
#include "bounds_kernel_eqO_xl_f2c_kernel.cu"
#include "bounds_kernel_eqY_xl_f2c_kernel.cu"
#include "bounds_kernel_eqZ_xl_f2c_kernel.cu"
#include "bounds_kernel_eqH_xr_f2c_kernel.cu"
#include "bounds_kernel_eqP_xr_f2c_kernel.cu"
#include "bounds_kernel_eqQ_xr_f2c_kernel.cu"
#include "bounds_kernel_eqI_xr_f2c_kernel.cu"
#include "bounds_kernel_eqJ_xr_f2c_kernel.cu"
#include "bounds_kernel_eqR_xr_f2c_kernel.cu"
#include "bounds_kernel_eqS_xr_f2c_kernel.cu"
#include "bounds_kernel_eqK_xr_f2c_kernel.cu"
#include "bounds_kernel_eqT_xr_f2c_kernel.cu"
#include "bounds_kernel_eqL_xr_f2c_kernel.cu"
#include "bounds_kernel_eqU_xr_f2c_kernel.cu"
#include "bounds_kernel_eqV_xr_f2c_kernel.cu"
#include "bounds_kernel_eqAF_xr_f2c_kernel.cu"
#include "bounds_kernel_eqM_xr_f2c_kernel.cu"
#include "bounds_kernel_eqW_xr_f2c_kernel.cu"
#include "bounds_kernel_eqX_xr_f2c_kernel.cu"
#include "bounds_kernel_eqN_xr_f2c_kernel.cu"
#include "bounds_kernel_eqO_xr_f2c_kernel.cu"
#include "bounds_kernel_eqY_xr_f2c_kernel.cu"
#include "bounds_kernel_eqZ_xr_f2c_kernel.cu"
#include "bounds_kernel_eqA_ydir_f2c_kernel.cu"
#include "bounds_kernel_eqB_ydir_f2c_kernel.cu"
#include "bounds_kernel_eqC_ydir_f2c_kernel.cu"
#include "bounds_kernel_eqAA_ydir_f2c_kernel.cu"
#include "copy_kernel_yydir_f2c_kernel.cu"
#include "bounds_kernel_eqD_ydir_f2c_kernel.cu"
#include "bounds_kernel_eqF_ydir_f2c_kernel.cu"
#include "bounds_kernel_eqH_yl_f2c_kernel.cu"
#include "bounds_kernel_eqP_yl_f2c_kernel.cu"
#include "bounds_kernel_eqQ_yl_f2c_kernel.cu"
#include "bounds_kernel_eqI_yl_f2c_kernel.cu"
#include "bounds_kernel_eqJ_yl_f2c_kernel.cu"
#include "bounds_kernel_eqR_yl_f2c_kernel.cu"
#include "bounds_kernel_eqS_yl_f2c_kernel.cu"
#include "bounds_kernel_eqE_ydir_f2c_kernel.cu"
#include "bounds_kernel_eqG_ydir_f2c_kernel.cu"
#include "bounds_kernel_eqK_yl_f2c_kernel.cu"
#include "bounds_kernel_eqT_yl_f2c_kernel.cu"
#include "bounds_kernel_eqL_yl_f2c_kernel.cu"
#include "bounds_kernel_eqU_yl_f2c_kernel.cu"
#include "bounds_kernel_eqV_yl_f2c_kernel.cu"
#include "bounds_kernel_eqAF_yl_f2c_kernel.cu"
#include "bounds_kernel_eqAG_ydir_f2c_kernel.cu"
#include "bounds_kernel_eqAH_ydir_f2c_kernel.cu"
#include "bounds_kernel_eqM_yl_f2c_kernel.cu"
#include "bounds_kernel_eqW_yl_f2c_kernel.cu"
#include "bounds_kernel_eqX_yl_f2c_kernel.cu"
#include "bounds_kernel_eqN_yl_f2c_kernel.cu"
#include "bounds_kernel_eqO_yl_f2c_kernel.cu"
#include "bounds_kernel_eqY_yl_f2c_kernel.cu"
#include "bounds_kernel_eqZ_yl_f2c_kernel.cu"
#include "bounds_kernel_eqH_yr_f2c_kernel.cu"
#include "bounds_kernel_eqP_yr_f2c_kernel.cu"
#include "bounds_kernel_eqQ_yr_f2c_kernel.cu"
#include "bounds_kernel_eqI_yr_f2c_kernel.cu"
#include "bounds_kernel_eqJ_yr_f2c_kernel.cu"
#include "bounds_kernel_eqR_yr_f2c_kernel.cu"
#include "bounds_kernel_eqS_yr_f2c_kernel.cu"
#include "bounds_kernel_eqK_yr_f2c_kernel.cu"
#include "bounds_kernel_eqT_yr_f2c_kernel.cu"
#include "bounds_kernel_eqL_yr_f2c_kernel.cu"
#include "bounds_kernel_eqU_yr_f2c_kernel.cu"
#include "bounds_kernel_eqV_yr_f2c_kernel.cu"
#include "bounds_kernel_eqAF_yr_f2c_kernel.cu"
#include "bounds_kernel_eqM_yr_f2c_kernel.cu"
#include "bounds_kernel_eqW_yr_f2c_kernel.cu"
#include "bounds_kernel_eqX_yr_f2c_kernel.cu"
#include "bounds_kernel_eqN_yr_f2c_kernel.cu"
#include "bounds_kernel_eqO_yr_f2c_kernel.cu"
#include "bounds_kernel_eqY_yr_f2c_kernel.cu"
#include "bounds_kernel_eqZ_yr_f2c_kernel.cu"
#include "bounds_kernel_eqA_zdir_f2c_kernel.cu"
#include "bounds_kernel_eqB_zdir_f2c_kernel.cu"
#include "bounds_kernel_eqC_zdir_f2c_kernel.cu"
#include "bounds_kernel_eqAA_zdir_f2c_kernel.cu"
#include "copy_kernel_zzdir_f2c_kernel.cu"
#include "bounds_kernel_eqD_zdir_f2c_kernel.cu"
#include "bounds_kernel_eqF_zdir_f2c_kernel.cu"
#include "bounds_kernel_eqH_zl_f2c_kernel.cu"
#include "bounds_kernel_eqP_zl_f2c_kernel.cu"
#include "bounds_kernel_eqQ_zl_f2c_kernel.cu"
#include "bounds_kernel_eqI_zl_f2c_kernel.cu"
#include "bounds_kernel_eqJ_zl_f2c_kernel.cu"
#include "bounds_kernel_eqR_zl_f2c_kernel.cu"
#include "bounds_kernel_eqS_zl_f2c_kernel.cu"
#include "bounds_kernel_eqE_zdir_f2c_kernel.cu"
#include "bounds_kernel_eqG_zdir_f2c_kernel.cu"
#include "bounds_kernel_eqK_zl_f2c_kernel.cu"
#include "bounds_kernel_eqT_zl_f2c_kernel.cu"
#include "bounds_kernel_eqL_zl_f2c_kernel.cu"
#include "bounds_kernel_eqU_zl_f2c_kernel.cu"
#include "bounds_kernel_eqV_zl_f2c_kernel.cu"
#include "bounds_kernel_eqAF_zl_f2c_kernel.cu"
#include "bounds_kernel_eqAG_zdir_f2c_kernel.cu"
#include "bounds_kernel_eqAH_zdir_f2c_kernel.cu"
#include "bounds_kernel_eqM_zl_f2c_kernel.cu"
#include "bounds_kernel_eqW_zl_f2c_kernel.cu"
#include "bounds_kernel_eqX_zl_f2c_kernel.cu"
#include "bounds_kernel_eqN_zl_f2c_kernel.cu"
#include "bounds_kernel_eqO_zl_f2c_kernel.cu"
#include "bounds_kernel_eqY_zl_f2c_kernel.cu"
#include "bounds_kernel_eqZ_zl_f2c_kernel.cu"
#include "bounds_kernel_eqH_zr_f2c_kernel.cu"
#include "bounds_kernel_eqP_zr_f2c_kernel.cu"
#include "bounds_kernel_eqQ_zr_f2c_kernel.cu"
#include "bounds_kernel_eqI_zr_f2c_kernel.cu"
#include "bounds_kernel_eqJ_zr_f2c_kernel.cu"
#include "bounds_kernel_eqR_zr_f2c_kernel.cu"
#include "bounds_kernel_eqS_zr_f2c_kernel.cu"
#include "bounds_kernel_eqK_zr_f2c_kernel.cu"
#include "bounds_kernel_eqT_zr_f2c_kernel.cu"
#include "bounds_kernel_eqL_zr_f2c_kernel.cu"
#include "bounds_kernel_eqU_zr_f2c_kernel.cu"
#include "bounds_kernel_eqV_zr_f2c_kernel.cu"
#include "bounds_kernel_eqAF_zr_f2c_kernel.cu"
#include "bounds_kernel_eqM_zr_f2c_kernel.cu"
#include "bounds_kernel_eqW_zr_f2c_kernel.cu"
#include "bounds_kernel_eqX_zr_f2c_kernel.cu"
#include "bounds_kernel_eqN_zr_f2c_kernel.cu"
#include "bounds_kernel_eqO_zr_f2c_kernel.cu"
#include "bounds_kernel_eqY_zr_f2c_kernel.cu"
#include "bounds_kernel_eqZ_zr_f2c_kernel.cu"
#include "boundt_kernel_eqE_xdir_f2c_kernel.cu"
#include "bountt_kernel_eqA_xdir_f2c_kernel.cu"
#include "bountt_kernel_eqB_xdir_f2c_kernel.cu"
#include "bountt_kernel_eqF_xdir_f2c_kernel.cu"
#include "bountt_kernel_eqD_f2c_kernel.cu"
#include "bountt_kernel_eqC_xdir_f2c_kernel.cu"
#include "bountt_kernel_eqE_xdir_f2c_kernel.cu"
#include "bountt_kernel_eqG_xyz_f2c_kernel.cu"
#include "boundt_kernel_eqG_xdir_f2c_kernel.cu"
#include "bountt_kernel_eqH_xdir_f2c_kernel.cu"
#include "boundt_kernel_eqE_ydir_f2c_kernel.cu"
#include "bountt_kernel_eqA_ydir_f2c_kernel.cu"
#include "bountt_kernel_eqB_ydir_f2c_kernel.cu"
#include "bountt_kernel_eqF_ydir_f2c_kernel.cu"
#include "bountt_kernel_eqC_ydir_f2c_kernel.cu"
#include "bountt_kernel_eqE_ydir_f2c_kernel.cu"
#include "boundt_kernel_eqG_ydir_f2c_kernel.cu"
#include "bountt_kernel_eqH_ydir_f2c_kernel.cu"
#include "boundt_kernel_eqE_zdir_f2c_kernel.cu"
#include "bountt_kernel_eqA_zdir_f2c_kernel.cu"
#include "bountt_kernel_eqB_zdir_f2c_kernel.cu"
#include "bountt_kernel_eqF_zdir_f2c_kernel.cu"
#include "bountt_kernel_eqC_zdir_f2c_kernel.cu"
#include "bountt_kernel_eqE_zdir_f2c_kernel.cu"
#include "boundt_kernel_eqG_zdir_f2c_kernel.cu"
#include "bountt_kernel_eqH_zdir_f2c_kernel.cu"
#include "boundt_kernel_eqA_xdir_f2c_kernel.cu"
#include "boundt_kernel_eqB_xdir_f2c_kernel.cu"
#include "boundt_kernel_eqF_xdir_f2c_kernel.cu"
#include "boundt_kernel_eqC_xdir_f2c_kernel.cu"
#include "boundt_kernel_eqD_xdir_f2c_kernel.cu"
#include "boundt_kernel_eqH_xyz_f2c_kernel.cu"
#include "boundt_kernel_eqA_ydir_f2c_kernel.cu"
#include "boundt_kernel_eqB_ydir_f2c_kernel.cu"
#include "boundt_kernel_eqF_ydir_f2c_kernel.cu"
#include "boundt_kernel_eqC_ydir_f2c_kernel.cu"
#include "boundt_kernel_eqD_ydir_f2c_kernel.cu"
#include "boundt_kernel_eqA_zdir_f2c_kernel.cu"
#include "boundt_kernel_eqB_zdir_f2c_kernel.cu"
#include "boundt_kernel_eqF_zdir_f2c_kernel.cu"
#include "boundt_kernel_eqC_zdir_f2c_kernel.cu"
#include "boundt_kernel_eqD_zdir_f2c_kernel.cu"
#include "radcal_kernel_meancoef_f2c_kernel.cu"
#include "radcal_kernel_addspecies_f2c_kernel.cu"
#include "radcal_kernel_addradiation_f2c_kernel.cu"
#include "maths_kernel_eqJ_f2c_kernel.cu"
#include "maths_kernel_eqAW_f2c_kernel.cu"
#include "maths_kernel_eqAX_f2c_kernel.cu"
#include "maths_kernel_eqAY_f2c_kernel.cu"
#include "maths_kernel_eqAZ_f2c_kernel.cu"
#include "maths_kernel_eqP_f2c_kernel.cu"
#include "maths_kernel_eqBM_f2c_kernel.cu"
#include "maths_kernel_eqBN_f2c_kernel.cu"
#include "maths_kernel_eqB_f2c_kernel.cu"
#include "maths_kernel_eqF_f2c_kernel.cu"
#include "maths_kernel_eqE_f2c_kernel.cu"
#include "maths_kernel_eqD_f2c_kernel.cu"
#include "maths_kernel_eqK_f2c_kernel.cu"
#include "maths_kernel_eqU_f2c_kernel.cu"
#include "maths_kernel_eqBO_f2c_kernel.cu"
#include "maths_kernel_eqBP_f2c_kernel.cu"
#include "maths_kernel_eqAU_f2c_kernel.cu"
#include "maths_kernel_eqV_fused_f2c_kernel.cu"
#include "turbin_kernel_eqA_f2c_kernel.cu"
#include "turbin_kernel_eqB_f2c_kernel.cu"
#include "turbin_kernel_eqC_f2c_kernel.cu"
#include "temper_kernel_eqA_f2c_kernel.cu"
#include "set_zero_kernel_MD5_f2c_kernel.cu"
#include "temper_kernel_eqB_f2c_kernel.cu"
#include "temper_kernel_eqC_f2c_kernel.cu"
#include "temper_kernel_eqD_f2c_kernel.cu"
#include "temper_kernel_eqE_f2c_kernel.cu"
#include "temper_kernel_eqF_f2c_kernel.cu"
#include "tempin_kernel_main_f2c_kernel.cu"
#include "copy_kernel_mdim_to_sdim_f2c_kernel.cu"


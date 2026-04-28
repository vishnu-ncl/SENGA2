// Auto-generated at 2026-04-28 18:43:38.388186 by ops-translator

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

//  global constants
int ncofmx;
int ntinmx;
int nspcmx;
int nssmax;
int nstpmx;
int ndcfmx;
int nvcfmx;
int nccfmx;
int nrkmax;
int ncbcsz;
int nbcprr;
int nspimx;
int ntbase;
int nintmx;
int nctmax;
int nctmm1;
int nrsmax;
int nbcpri;
int ncfrmx;
double acoffx;
double bcoffx;
double ccoffx;
double dcoffx;
double ecoffx;
double acof1x;
double bcof1x;
double ccof1x;
double dcof1x;
double acof2x;
double bcof2x;
double ccof2x;
double dcof2x;
double acof3x;
double bcof3x;
double acof4x;
double bcof4x;
double ccof4x;
double acof5x;
double bcof5x;
double ccof5x;
double dcof5x;
double ovdelx;
double acofsx;
double bcofsx;
double ccofsx;
double dcofsx;
double ecofsx;
double acfs1x;
double bcfs1x;
double ccfs1x;
double dcfs1x;
double ecfs1x;
double acfs2x;
double bcfs2x;
double ccfs2x;
double dcfs2x;
double ecfs2x;
double acfs3x;
double bcfs3x;
double acfs4x;
double bcfs4x;
double ccfs4x;
double acfs5x;
double bcfs5x;
double ccfs5x;
double dcfs5x;
double ovdlx2;
double acoffy;
double bcoffy;
double ccoffy;
double dcoffy;
double ecoffy;
double acof1y;
double bcof1y;
double ccof1y;
double dcof1y;
double acof2y;
double bcof2y;
double ccof2y;
double dcof2y;
double acof3y;
double bcof3y;
double acof4y;
double bcof4y;
double ccof4y;
double acof5y;
double bcof5y;
double ccof5y;
double dcof5y;
double ovdely;
double acofsy;
double bcofsy;
double ccofsy;
double dcofsy;
double ecofsy;
double acfs1y;
double bcfs1y;
double ccfs1y;
double dcfs1y;
double ecfs1y;
double acfs2y;
double bcfs2y;
double ccfs2y;
double dcfs2y;
double ecfs2y;
double acfs3y;
double bcfs3y;
double acfs4y;
double bcfs4y;
double ccfs4y;
double acfs5y;
double bcfs5y;
double ccfs5y;
double dcfs5y;
double ovdly2;
double acoffz;
double bcoffz;
double ccoffz;
double dcoffz;
double ecoffz;
double acof1z;
double bcof1z;
double ccof1z;
double dcof1z;
double acof2z;
double bcof2z;
double ccof2z;
double dcof2z;
double acof3z;
double bcof3z;
double acof4z;
double bcof4z;
double ccof4z;
double acof5z;
double bcof5z;
double ccof5z;
double dcof5z;
double ovdelz;
double acofsz;
double bcofsz;
double ccofsz;
double dcofsz;
double ecofsz;
double acfs1z;
double bcfs1z;
double ccfs1z;
double dcfs1z;
double ecfs1z;
double acfs2z;
double bcfs2z;
double ccfs2z;
double dcfs2z;
double ecfs2z;
double acfs3z;
double bcfs3z;
double acfs4z;
double bcfs4z;
double ccfs4z;
double acfs5z;
double bcfs5z;
double ccfs5z;
double dcfs5z;
double ovdlz2;
double acofx1;
double bcofx1;
double acofy1;
double bcofy1;
double acofz1;
double bcofz1;
double acofxy;
double bcofxy;
double ccofxy;
double dcofxy;
double ecofxy;
double acf1xy;
double bcf1xy;
double ccf1xy;
double dcf1xy;
double acf2xy;
double bcf2xy;
double ccf2xy;
double dcf2xy;
double acf3xy;
double bcf3xy;
double acf4xy;
double bcf4xy;
double ccf4xy;
double acf5xy;
double bcf5xy;
double ccf5xy;
double dcf5xy;
double acc1xy;
double bcc1xy;
double ccc1xy;
double dcc1xy;
double acc2xy;
double bcc2xy;
double ccc2xy;
double dcc2xy;
double acofxz;
double bcofxz;
double ccofxz;
double dcofxz;
double ecofxz;
double acf1xz;
double bcf1xz;
double ccf1xz;
double dcf1xz;
double acf2xz;
double bcf2xz;
double ccf2xz;
double dcf2xz;
double acf3xz;
double bcf3xz;
double acf4xz;
double bcf4xz;
double ccf4xz;
double acf5xz;
double bcf5xz;
double ccf5xz;
double dcf5xz;
double acc1xz;
double bcc1xz;
double ccc1xz;
double dcc1xz;
double acc2xz;
double bcc2xz;
double ccc2xz;
double dcc2xz;
double acofyz;
double bcofyz;
double ccofyz;
double dcofyz;
double ecofyz;
double acf1yz;
double bcf1yz;
double ccf1yz;
double dcf1yz;
double acf2yz;
double bcf2yz;
double ccf2yz;
double dcf2yz;
double acf3yz;
double bcf3yz;
double acf4yz;
double bcf4yz;
double ccf4yz;
double acf5yz;
double bcf5yz;
double ccf5yz;
double dcf5yz;
double acc1yz;
double bcc1yz;
double ccc1yz;
double dcc1yz;
double acc2yz;
double bcc2yz;
double ccc2yz;
double dcc2yz;
double foursb;
double trfrth;
double rlamda;
double alamda;

void ops_init_backend(){}

extern "C"
void ops_decl_const_f2c(char const *name, int dim, int size, char *dat) {

    ops_execute(OPS_instance::getOPSInstance());

    if(!strcmp(name, "ncofmx")) {
        ncofmx = *(int *)dat;
    }
    else
    if(!strcmp(name, "ntinmx")) {
        ntinmx = *(int *)dat;
    }
    else
    if(!strcmp(name, "nspcmx")) {
        nspcmx = *(int *)dat;
    }
    else
    if(!strcmp(name, "nssmax")) {
        nssmax = *(int *)dat;
    }
    else
    if(!strcmp(name, "nstpmx")) {
        nstpmx = *(int *)dat;
    }
    else
    if(!strcmp(name, "ndcfmx")) {
        ndcfmx = *(int *)dat;
    }
    else
    if(!strcmp(name, "nvcfmx")) {
        nvcfmx = *(int *)dat;
    }
    else
    if(!strcmp(name, "nccfmx")) {
        nccfmx = *(int *)dat;
    }
    else
    if(!strcmp(name, "nrkmax")) {
        nrkmax = *(int *)dat;
    }
    else
    if(!strcmp(name, "ncbcsz")) {
        ncbcsz = *(int *)dat;
    }
    else
    if(!strcmp(name, "nbcprr")) {
        nbcprr = *(int *)dat;
    }
    else
    if(!strcmp(name, "nspimx")) {
        nspimx = *(int *)dat;
    }
    else
    if(!strcmp(name, "ntbase")) {
        ntbase = *(int *)dat;
    }
    else
    if(!strcmp(name, "nintmx")) {
        nintmx = *(int *)dat;
    }
    else
    if(!strcmp(name, "nctmax")) {
        nctmax = *(int *)dat;
    }
    else
    if(!strcmp(name, "nctmm1")) {
        nctmm1 = *(int *)dat;
    }
    else
    if(!strcmp(name, "nrsmax")) {
        nrsmax = *(int *)dat;
    }
    else
    if(!strcmp(name, "nbcpri")) {
        nbcpri = *(int *)dat;
    }
    else
    if(!strcmp(name, "ncfrmx")) {
        ncfrmx = *(int *)dat;
    }
    else
    if(!strcmp(name, "acoffx")) {
        acoffx = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcoffx")) {
        bcoffx = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccoffx")) {
        ccoffx = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcoffx")) {
        dcoffx = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecoffx")) {
        ecoffx = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof1x")) {
        acof1x = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof1x")) {
        bcof1x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof1x")) {
        ccof1x = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcof1x")) {
        dcof1x = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof2x")) {
        acof2x = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof2x")) {
        bcof2x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof2x")) {
        ccof2x = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcof2x")) {
        dcof2x = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof3x")) {
        acof3x = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof3x")) {
        bcof3x = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof4x")) {
        acof4x = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof4x")) {
        bcof4x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof4x")) {
        ccof4x = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof5x")) {
        acof5x = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof5x")) {
        bcof5x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof5x")) {
        ccof5x = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcof5x")) {
        dcof5x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ovdelx")) {
        ovdelx = *(double *)dat;
    }
    else
    if(!strcmp(name, "acofsx")) {
        acofsx = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcofsx")) {
        bcofsx = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccofsx")) {
        ccofsx = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcofsx")) {
        dcofsx = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecofsx")) {
        ecofsx = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs1x")) {
        acfs1x = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs1x")) {
        bcfs1x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs1x")) {
        ccfs1x = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcfs1x")) {
        dcfs1x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecfs1x")) {
        ecfs1x = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs2x")) {
        acfs2x = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs2x")) {
        bcfs2x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs2x")) {
        ccfs2x = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcfs2x")) {
        dcfs2x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecfs2x")) {
        ecfs2x = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs3x")) {
        acfs3x = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs3x")) {
        bcfs3x = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs4x")) {
        acfs4x = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs4x")) {
        bcfs4x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs4x")) {
        ccfs4x = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs5x")) {
        acfs5x = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs5x")) {
        bcfs5x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs5x")) {
        ccfs5x = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcfs5x")) {
        dcfs5x = *(double *)dat;
    }
    else
    if(!strcmp(name, "ovdlx2")) {
        ovdlx2 = *(double *)dat;
    }
    else
    if(!strcmp(name, "acoffy")) {
        acoffy = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcoffy")) {
        bcoffy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccoffy")) {
        ccoffy = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcoffy")) {
        dcoffy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecoffy")) {
        ecoffy = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof1y")) {
        acof1y = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof1y")) {
        bcof1y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof1y")) {
        ccof1y = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcof1y")) {
        dcof1y = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof2y")) {
        acof2y = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof2y")) {
        bcof2y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof2y")) {
        ccof2y = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcof2y")) {
        dcof2y = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof3y")) {
        acof3y = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof3y")) {
        bcof3y = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof4y")) {
        acof4y = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof4y")) {
        bcof4y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof4y")) {
        ccof4y = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof5y")) {
        acof5y = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof5y")) {
        bcof5y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof5y")) {
        ccof5y = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcof5y")) {
        dcof5y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ovdely")) {
        ovdely = *(double *)dat;
    }
    else
    if(!strcmp(name, "acofsy")) {
        acofsy = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcofsy")) {
        bcofsy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccofsy")) {
        ccofsy = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcofsy")) {
        dcofsy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecofsy")) {
        ecofsy = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs1y")) {
        acfs1y = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs1y")) {
        bcfs1y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs1y")) {
        ccfs1y = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcfs1y")) {
        dcfs1y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecfs1y")) {
        ecfs1y = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs2y")) {
        acfs2y = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs2y")) {
        bcfs2y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs2y")) {
        ccfs2y = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcfs2y")) {
        dcfs2y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecfs2y")) {
        ecfs2y = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs3y")) {
        acfs3y = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs3y")) {
        bcfs3y = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs4y")) {
        acfs4y = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs4y")) {
        bcfs4y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs4y")) {
        ccfs4y = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs5y")) {
        acfs5y = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs5y")) {
        bcfs5y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs5y")) {
        ccfs5y = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcfs5y")) {
        dcfs5y = *(double *)dat;
    }
    else
    if(!strcmp(name, "ovdly2")) {
        ovdly2 = *(double *)dat;
    }
    else
    if(!strcmp(name, "acoffz")) {
        acoffz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcoffz")) {
        bcoffz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccoffz")) {
        ccoffz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcoffz")) {
        dcoffz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecoffz")) {
        ecoffz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof1z")) {
        acof1z = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof1z")) {
        bcof1z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof1z")) {
        ccof1z = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcof1z")) {
        dcof1z = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof2z")) {
        acof2z = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof2z")) {
        bcof2z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof2z")) {
        ccof2z = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcof2z")) {
        dcof2z = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof3z")) {
        acof3z = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof3z")) {
        bcof3z = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof4z")) {
        acof4z = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof4z")) {
        bcof4z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof4z")) {
        ccof4z = *(double *)dat;
    }
    else
    if(!strcmp(name, "acof5z")) {
        acof5z = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcof5z")) {
        bcof5z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccof5z")) {
        ccof5z = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcof5z")) {
        dcof5z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ovdelz")) {
        ovdelz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acofsz")) {
        acofsz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcofsz")) {
        bcofsz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccofsz")) {
        ccofsz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcofsz")) {
        dcofsz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecofsz")) {
        ecofsz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs1z")) {
        acfs1z = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs1z")) {
        bcfs1z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs1z")) {
        ccfs1z = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcfs1z")) {
        dcfs1z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecfs1z")) {
        ecfs1z = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs2z")) {
        acfs2z = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs2z")) {
        bcfs2z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs2z")) {
        ccfs2z = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcfs2z")) {
        dcfs2z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecfs2z")) {
        ecfs2z = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs3z")) {
        acfs3z = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs3z")) {
        bcfs3z = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs4z")) {
        acfs4z = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs4z")) {
        bcfs4z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs4z")) {
        ccfs4z = *(double *)dat;
    }
    else
    if(!strcmp(name, "acfs5z")) {
        acfs5z = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcfs5z")) {
        bcfs5z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccfs5z")) {
        ccfs5z = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcfs5z")) {
        dcfs5z = *(double *)dat;
    }
    else
    if(!strcmp(name, "ovdlz2")) {
        ovdlz2 = *(double *)dat;
    }
    else
    if(!strcmp(name, "acofx1")) {
        acofx1 = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcofx1")) {
        bcofx1 = *(double *)dat;
    }
    else
    if(!strcmp(name, "acofy1")) {
        acofy1 = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcofy1")) {
        bcofy1 = *(double *)dat;
    }
    else
    if(!strcmp(name, "acofz1")) {
        acofz1 = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcofz1")) {
        bcofz1 = *(double *)dat;
    }
    else
    if(!strcmp(name, "acofxy")) {
        acofxy = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcofxy")) {
        bcofxy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccofxy")) {
        ccofxy = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcofxy")) {
        dcofxy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecofxy")) {
        ecofxy = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf1xy")) {
        acf1xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf1xy")) {
        bcf1xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf1xy")) {
        ccf1xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcf1xy")) {
        dcf1xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf2xy")) {
        acf2xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf2xy")) {
        bcf2xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf2xy")) {
        ccf2xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcf2xy")) {
        dcf2xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf3xy")) {
        acf3xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf3xy")) {
        bcf3xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf4xy")) {
        acf4xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf4xy")) {
        bcf4xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf4xy")) {
        ccf4xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf5xy")) {
        acf5xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf5xy")) {
        bcf5xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf5xy")) {
        ccf5xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcf5xy")) {
        dcf5xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "acc1xy")) {
        acc1xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcc1xy")) {
        bcc1xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccc1xy")) {
        ccc1xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcc1xy")) {
        dcc1xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "acc2xy")) {
        acc2xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcc2xy")) {
        bcc2xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccc2xy")) {
        ccc2xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcc2xy")) {
        dcc2xy = *(double *)dat;
    }
    else
    if(!strcmp(name, "acofxz")) {
        acofxz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcofxz")) {
        bcofxz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccofxz")) {
        ccofxz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcofxz")) {
        dcofxz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecofxz")) {
        ecofxz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf1xz")) {
        acf1xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf1xz")) {
        bcf1xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf1xz")) {
        ccf1xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcf1xz")) {
        dcf1xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf2xz")) {
        acf2xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf2xz")) {
        bcf2xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf2xz")) {
        ccf2xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcf2xz")) {
        dcf2xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf3xz")) {
        acf3xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf3xz")) {
        bcf3xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf4xz")) {
        acf4xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf4xz")) {
        bcf4xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf4xz")) {
        ccf4xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf5xz")) {
        acf5xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf5xz")) {
        bcf5xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf5xz")) {
        ccf5xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcf5xz")) {
        dcf5xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acc1xz")) {
        acc1xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcc1xz")) {
        bcc1xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccc1xz")) {
        ccc1xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcc1xz")) {
        dcc1xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acc2xz")) {
        acc2xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcc2xz")) {
        bcc2xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccc2xz")) {
        ccc2xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcc2xz")) {
        dcc2xz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acofyz")) {
        acofyz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcofyz")) {
        bcofyz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccofyz")) {
        ccofyz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcofyz")) {
        dcofyz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ecofyz")) {
        ecofyz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf1yz")) {
        acf1yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf1yz")) {
        bcf1yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf1yz")) {
        ccf1yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcf1yz")) {
        dcf1yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf2yz")) {
        acf2yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf2yz")) {
        bcf2yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf2yz")) {
        ccf2yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcf2yz")) {
        dcf2yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf3yz")) {
        acf3yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf3yz")) {
        bcf3yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf4yz")) {
        acf4yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf4yz")) {
        bcf4yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf4yz")) {
        ccf4yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acf5yz")) {
        acf5yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcf5yz")) {
        bcf5yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccf5yz")) {
        ccf5yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcf5yz")) {
        dcf5yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acc1yz")) {
        acc1yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcc1yz")) {
        bcc1yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccc1yz")) {
        ccc1yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcc1yz")) {
        dcc1yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "acc2yz")) {
        acc2yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "bcc2yz")) {
        bcc2yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "ccc2yz")) {
        ccc2yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "dcc2yz")) {
        dcc2yz = *(double *)dat;
    }
    else
    if(!strcmp(name, "foursb")) {
        foursb = *(double *)dat;
    }
    else
    if(!strcmp(name, "trfrth")) {
        trfrth = *(double *)dat;
    }
    else
    if(!strcmp(name, "rlamda")) {
        rlamda = *(double *)dat;
    }
    else
    if(!strcmp(name, "alamda")) {
        alamda = *(double *)dat;
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


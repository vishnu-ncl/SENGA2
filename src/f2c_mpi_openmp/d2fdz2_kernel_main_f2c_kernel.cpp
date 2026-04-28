// Auto-generated at 2026-04-28 18:43:15.565430 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void d2fdz2_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void d2fdz2_kernel_main_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[7];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
    args[4] = desc->args[4];
    args[5] = desc->args[5];
    args[6] = desc->args[6];
#endif

//  ======
//  Timing
//  ======
    double __t1, __t2, __c1, __c2;

    ops_arg arg0 = args[0];
    ops_arg arg1 = args[1];
    ops_arg arg2 = args[2];
    ops_arg arg3 = args[3];
    ops_arg arg4 = args[4];
    ops_arg arg5 = args[5];
    ops_arg arg6 = args[6];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 7, range, 15)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 15, "d2fdz2_kernel_main");
        block->instance->OPS_kernels[15].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "d2fdz2_kernel_main");
#endif

//  =================================================
//  compute locally allocated range for the sub-block
//  =================================================
    int start_indx[3];
    int end_indx[3];
    int arg_idx[3];

//  Range here is in C-style while start_indx and end_indx is Fortran style
#if defined(OPS_MPI) && !defined(OPS_LAZY)
    if ( getRange(block, start_indx, end_indx, range) < 0 ) return;
#else
    for (int n = 0; n < 3; n++) {
        start_indx[n] = range[2*n] + 1;
        end_indx[n]   = range[2*n+1];
    }
#endif

#ifdef OPS_MPI
    getIdx(block, start_indx, arg_idx);
#else
    arg_idx[0] = start_indx[0];
    arg_idx[1] = start_indx[1];
    arg_idx[2] = start_indx[2];
#endif

//  ======================================================
//  Initialize global variable with the dimensions of dats
//  ======================================================
    int xdim0_d2fdz2_kernel_main = args[0].dat->size[0];
    int ydim0_d2fdz2_kernel_main = args[0].dat->size[1];
    int xdim1_d2fdz2_kernel_main = args[1].dat->size[0];
    int ydim1_d2fdz2_kernel_main = args[1].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ functn_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ fderiv_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int * __restrict__  nzglbl = (int *)args[2].data;
    int * __restrict__  nendzl = (int *)args[3].data;
    int * __restrict__  nendzr = (int *)args[4].data;
    int * __restrict__  nbound = (int *)args[5].data;

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 7);
    ops_halo_exchanges(args, 7, range);
    ops_H_D_exchanges_host(args, 7);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[15].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {
                int idx[] = {arg_idx[0] + n_x, arg_idx[1] + n_y, arg_idx[2] + n_z};

                const  ACC<double> functn(xdim0_d2fdz2_kernel_main, ydim0_d2fdz2_kernel_main, functn_p + (n_x * 1) + (n_y * xdim0_d2fdz2_kernel_main * 1) + (n_z * xdim0_d2fdz2_kernel_main * ydim0_d2fdz2_kernel_main * 1));

                 ACC<double> fderiv(xdim1_d2fdz2_kernel_main, ydim1_d2fdz2_kernel_main, fderiv_p + (n_x * 1) + (n_y * xdim1_d2fdz2_kernel_main * 1) + (n_z * xdim1_d2fdz2_kernel_main * ydim1_d2fdz2_kernel_main * 1));

    double fdifap;
    double fdifbp;
    double fdifcp;
    double fdifdp;
    double fdifep;
    double fdifam;
    double fdifbm;
    double fdifcm;
    double fdifdm;
    double fdifem;
    int kc;
    int kstart;
    int kfinis;

    kstart = 1;
    kfinis = nzglbl[0];
    if (nendzl[0] == nbound[0]) {
        kstart = 6;
    }
    if (nendzr[0] == nbound[0]) {
        kfinis = nzglbl[0] - 5;
    }
    kc = idx[2];
    if (kc >= kstart && kc <= kfinis) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(0, 0, -3);
        fdifdp = functn(0, 0, 4) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(0, 0, -4);
        fdifep = functn(0, 0, 5) - functn(0, 0, 0);
        fdifem = functn(0, 0, 0) - functn(0, 0, -5);
        fderiv(0, 0, 0) = acofsz * (fdifap - fdifam) + bcofsz * (fdifbp - fdifbm) + ccofsz * (fdifcp - fdifcm) + dcofsz * (fdifdp - fdifdm) + ecofsz * (fdifep - fdifem);
    } else if (kc == 1) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifdp = functn(0, 0, 4) - functn(0, 0, 0);
        fdifep = functn(0, 0, 5) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs1z * fdifap + bcfs1z * fdifbp + ccfs1z * fdifcp + dcfs1z * fdifdp + ecfs1z * fdifep;
    } else if (kc == 2) {
        fdifap = functn(0, 0, -1) - functn(0, 0, 0);
        fdifbp = functn(0, 0, 1) - functn(0, 0, 0);
        fdifcp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifdp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifep = functn(0, 0, 4) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs2z * fdifap + bcfs2z * fdifbp + ccfs2z * fdifcp + dcfs2z * fdifdp + ecfs2z * fdifep;
    } else if (kc == 3) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fderiv(0, 0, 0) = acfs3z * (fdifap - fdifam) + bcfs3z * (fdifbp - fdifbm);
    } else if (kc == 4) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(0, 0, -3);
        fderiv(0, 0, 0) = acfs4z * (fdifap - fdifam) + bcfs4z * (fdifbp - fdifbm) + ccfs4z * (fdifcp - fdifcm);
    } else if (kc == 5) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(0, 0, -3);
        fdifdp = functn(0, 0, 4) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(0, 0, -4);
        fderiv(0, 0, 0) = acfs5z * (fdifap - fdifam) + bcfs5z * (fdifbp - fdifbm) + ccfs5z * (fdifcp - fdifcm) + dcfs5z * (fdifdp - fdifdm);
    } else if (kc == nzglbl[0] - 4) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(0, 0, -3);
        fdifdp = functn(0, 0, 4) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(0, 0, -4);
        fderiv(0, 0, 0) = acfs5z * (fdifap - fdifam) + bcfs5z * (fdifbp - fdifbm) + ccfs5z * (fdifcp - fdifcm) + dcfs5z * (fdifdp - fdifdm);
    } else if (kc == nzglbl[0] - 3) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(0, 0, -3);
        fderiv(0, 0, 0) = acfs4z * (fdifap - fdifam) + bcfs4z * (fdifbp - fdifbm) + ccfs4z * (fdifcp - fdifcm);
    } else if (kc == nzglbl[0] - 2) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fderiv(0, 0, 0) = acfs3z * (fdifap - fdifam) + bcfs3z * (fdifbp - fdifbm);
    } else if (kc == nzglbl[0] - 1) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifbp = functn(0, 0, -1) - functn(0, 0, 0);
        fdifcp = functn(0, 0, -2) - functn(0, 0, 0);
        fdifdp = functn(0, 0, -3) - functn(0, 0, 0);
        fdifep = functn(0, 0, -4) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs2z * fdifap + bcfs2z * fdifbp + ccfs2z * fdifcp + dcfs2z * fdifdp + ecfs2z * fdifep;
    } else if (kc == nzglbl[0]) {
        fdifap = functn(0, 0, -1) - functn(0, 0, 0);
        fdifbp = functn(0, 0, -2) - functn(0, 0, 0);
        fdifcp = functn(0, 0, -3) - functn(0, 0, 0);
        fdifdp = functn(0, 0, -4) - functn(0, 0, 0);
        fdifep = functn(0, 0, -5) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs1z * fdifap + bcfs1z * fdifbp + ccfs1z * fdifcp + dcfs1z * fdifdp + ecfs1z * fdifep;
    }
    fderiv(0, 0, 0) = fderiv(0, 0, 0) * ovdlz2;

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[15].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 7);
    ops_set_halo_dirtybit3(&args[1], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[15].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[15].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[15].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void d2fdz2_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("d2fdz2_kernel_main", args, 7, 15, dim, 0, range, block, d2fdz2_kernel_main_execute);
}
#endif

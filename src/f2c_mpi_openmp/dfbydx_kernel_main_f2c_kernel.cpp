// Auto-generated at 2026-04-28 18:43:14.692135 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void dfbydx_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void dfbydx_kernel_main_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 7, range, 6)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 6, "dfbydx_kernel_main");
        block->instance->OPS_kernels[6].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "dfbydx_kernel_main");
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
    int xdim0_dfbydx_kernel_main = args[0].dat->size[0];
    int ydim0_dfbydx_kernel_main = args[0].dat->size[1];
    int xdim1_dfbydx_kernel_main = args[1].dat->size[0];
    int ydim1_dfbydx_kernel_main = args[1].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ functn_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ fderiv_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int * __restrict__  nxglbl = (int *)args[2].data;
    int * __restrict__  nendxl = (int *)args[3].data;
    int * __restrict__  nendxr = (int *)args[4].data;
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
        block->instance->OPS_kernels[6].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {
                int idx[] = {arg_idx[0] + n_x, arg_idx[1] + n_y, arg_idx[2] + n_z};

                const  ACC<double> functn(xdim0_dfbydx_kernel_main, ydim0_dfbydx_kernel_main, functn_p + (n_x * 1) + (n_y * xdim0_dfbydx_kernel_main * 1) + (n_z * xdim0_dfbydx_kernel_main * ydim0_dfbydx_kernel_main * 1));

                 ACC<double> fderiv(xdim1_dfbydx_kernel_main, ydim1_dfbydx_kernel_main, fderiv_p + (n_x * 1) + (n_y * xdim1_dfbydx_kernel_main * 1) + (n_z * xdim1_dfbydx_kernel_main * ydim1_dfbydx_kernel_main * 1));

    double fdiffa;
    double fdiffb;
    double fdiffc;
    double fdiffd;
    double fdiffe;
    int ic;
    int istart;
    int ifinis;

    istart = 1;
    ifinis = nxglbl[0];
    if (nendxl[0] == nbound[0]) {
        istart = 6;
    }
    if (nendxr[0] == nbound[0]) {
        ifinis = nxglbl[0] - 5;
    }
    ic = idx[0];
    if (ic >= istart && ic <= ifinis) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(-3, 0, 0);
        fdiffd = functn(4, 0, 0) - functn(-4, 0, 0);
        fdiffe = functn(5, 0, 0) - functn(-5, 0, 0);
        fderiv(0, 0, 0) = acoffx * fdiffa + bcoffx * fdiffb + ccoffx * fdiffc + dcoffx * fdiffd + ecoffx * fdiffe;
    } else if (ic == 1) {
        fdiffa = functn(1, 0, 0) - functn(0, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(0, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(0, 0, 0);
        fdiffd = functn(4, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acof1x * fdiffa + bcof1x * fdiffb + ccof1x * fdiffc + dcof1x * fdiffd;
    } else if (ic == 2) {
        fdiffa = functn(-1, 0, 0) - functn(0, 0, 0);
        fdiffb = functn(1, 0, 0) - functn(0, 0, 0);
        fdiffc = functn(2, 0, 0) - functn(0, 0, 0);
        fdiffd = functn(3, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acof2x * fdiffa + bcof2x * fdiffb + ccof2x * fdiffc + dcof2x * fdiffd;
    } else if (ic == 3) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fderiv(0, 0, 0) = acof3x * fdiffa + bcof3x * fdiffb;
    } else if (ic == 4) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(-3, 0, 0);
        fderiv(0, 0, 0) = acof4x * fdiffa + bcof4x * fdiffb + ccof4x * fdiffc;
    } else if (ic == 5) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(-3, 0, 0);
        fdiffd = functn(4, 0, 0) - functn(-4, 0, 0);
        fderiv(0, 0, 0) = acof5x * fdiffa + bcof5x * fdiffb + ccof5x * fdiffc + dcof5x * fdiffd;
    } else if (ic == nxglbl[0] - 4) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(-3, 0, 0);
        fdiffd = functn(4, 0, 0) - functn(-4, 0, 0);
        fderiv(0, 0, 0) = acof5x * fdiffa + bcof5x * fdiffb + ccof5x * fdiffc + dcof5x * fdiffd;
    } else if (ic == nxglbl[0] - 3) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(-3, 0, 0);
        fderiv(0, 0, 0) = acof4x * fdiffa + bcof4x * fdiffb + ccof4x * fdiffc;
    } else if (ic == nxglbl[0] - 2) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fderiv(0, 0, 0) = acof3x * fdiffa + bcof3x * fdiffb;
    } else if (ic == nxglbl[0] - 1) {
        fdiffa = functn(0, 0, 0) - functn(1, 0, 0);
        fdiffb = functn(0, 0, 0) - functn(-1, 0, 0);
        fdiffc = functn(0, 0, 0) - functn(-2, 0, 0);
        fdiffd = functn(0, 0, 0) - functn(-3, 0, 0);
        fderiv(0, 0, 0) = acof2x * fdiffa + bcof2x * fdiffb + ccof2x * fdiffc + dcof2x * fdiffd;
    } else if (ic == nxglbl[0]) {
        fdiffa = functn(0, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(0, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(0, 0, 0) - functn(-3, 0, 0);
        fdiffd = functn(0, 0, 0) - functn(-4, 0, 0);
        fderiv(0, 0, 0) = acof1x * fdiffa + bcof1x * fdiffb + ccof1x * fdiffc + dcof1x * fdiffd;
    }
    fderiv(0, 0, 0) = fderiv(0, 0, 0) * ovdelx;

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[6].time += __t1 - __t2;
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
        block->instance->OPS_kernels[6].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[6].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[6].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void dfbydx_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("dfbydx_kernel_main", args, 7, 6, dim, 0, range, block, dfbydx_kernel_main_execute);
}
#endif

// Auto-generated at 2026-04-28 18:43:17.206174 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void d2fdxy_kernel_eqAG_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void d2fdxy_kernel_eqAG_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[4];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
#endif

//  ======
//  Timing
//  ======
    double __t1, __t2, __c1, __c2;

    ops_arg arg0 = args[0];
    ops_arg arg1 = args[1];
    ops_arg arg2 = args[2];
    ops_arg arg3 = args[3];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 4, range, 50)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 50, "d2fdxy_kernel_eqAG");
        block->instance->OPS_kernels[50].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "d2fdxy_kernel_eqAG");
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
    int xdim0_d2fdxy_kernel_eqAG = args[0].dat->size[0];
    int ydim0_d2fdxy_kernel_eqAG = args[0].dat->size[1];
    int xdim1_d2fdxy_kernel_eqAG = args[1].dat->size[0];
    int ydim1_d2fdxy_kernel_eqAG = args[1].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ functn_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ fderiv_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int * __restrict__  nxglbl = (int *)args[2].data;

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 4);
    ops_halo_exchanges(args, 4, range);
    ops_H_D_exchanges_host(args, 4);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[50].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {
                int idx[] = {arg_idx[0] + n_x, arg_idx[1] + n_y, arg_idx[2] + n_z};

                const  ACC<double> functn(xdim0_d2fdxy_kernel_eqAG, ydim0_d2fdxy_kernel_eqAG, functn_p + (n_x * 1) + (n_y * xdim0_d2fdxy_kernel_eqAG * 1) + (n_z * xdim0_d2fdxy_kernel_eqAG * ydim0_d2fdxy_kernel_eqAG * 1));

                 ACC<double> fderiv(xdim1_d2fdxy_kernel_eqAG, ydim1_d2fdxy_kernel_eqAG, fderiv_p + (n_x * 1) + (n_y * xdim1_d2fdxy_kernel_eqAG * 1) + (n_z * xdim1_d2fdxy_kernel_eqAG * ydim1_d2fdxy_kernel_eqAG * 1));

    double fdiffa;
    double fdiffb;
    double fdiffc;
    double fdiffd;
    double fstora;
    double fstorb;
    double fstorc;
    int ic;
    int jc;

    ic = idx[0];
    jc = idx[1];
    if (ic >= nxglbl[0] - 4 && ic <= nxglbl[0] - 2 && jc >= 3 && jc <= 5) {
        fdiffa = functn(1, 1, 0) - functn(1, -1, 0) - functn(-1, 1, 0) + functn(-1, -1, 0);
        fdiffb = functn(2, 2, 0) - functn(2, -2, 0) - functn(-2, 2, 0) + functn(-2, -2, 0);
        fderiv(0, 0, 0) = acf3xy * fdiffa + bcf3xy * fdiffb;
        fstora = fdiffa;
        fstorb = fdiffb;
    }
    if (ic >= nxglbl[0] - 4 && ic <= nxglbl[0] - 3 && jc >= 4 && jc <= 5) {
        fdiffc = functn(3, 3, 0) - functn(3, -3, 0) - functn(-3, 3, 0) + functn(-3, -3, 0);
        fderiv(0, 0, 0) = acf4xy * fstora + bcf4xy * fstorb + ccf4xy * fdiffc;
        fstorc = fdiffc;
    }
    if (ic == nxglbl[0] - 4 && jc == 5) {
        fdiffd = functn(4, 4, 0) - functn(4, -4, 0) - functn(-4, 4, 0) + functn(-4, -4, 0);
        fderiv(0, 0, 0) = acf5xy * fstora + bcf5xy * fstorb + ccf5xy * fstorc + dcf5xy * fdiffd;
    }

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[50].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 4);
    ops_set_halo_dirtybit3(&args[1], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[50].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[50].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[50].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void d2fdxy_kernel_eqAG_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("d2fdxy_kernel_eqAG", args, 4, 50, dim, 0, range, block, d2fdxy_kernel_eqAG_execute);
}
#endif

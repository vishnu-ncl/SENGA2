// Auto-generated at 2026-04-28 18:43:37.676302 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void set_zero_kernel_MD5_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void set_zero_kernel_MD5_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[1];
    args[0] = desc->args[0];
#endif

//  ======
//  Timing
//  ======
    double __t1, __t2, __c1, __c2;

    ops_arg arg0 = args[0];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 1, range, 559)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 559, "set_zero_kernel_MD5");
        block->instance->OPS_kernels[559].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "set_zero_kernel_MD5");
#endif

//  =================================================
//  compute locally allocated range for the sub-block
//  =================================================
    int start_indx[3];
    int end_indx[3];

//  Range here is in C-style while start_indx and end_indx is Fortran style
#if defined(OPS_MPI) && !defined(OPS_LAZY)
    if ( getRange(block, start_indx, end_indx, range) < 0 ) return;
#else
    for (int n = 0; n < 3; n++) {
        start_indx[n] = range[2*n] + 1;
        end_indx[n]   = range[2*n+1];
    }
#endif

//  ======================================================
//  Initialize global variable with the dimensions of dats
//  ======================================================
    int xdim0_set_zero_kernel_MD5 = args[0].dat->size[0];
    int ydim0_set_zero_kernel_MD5 = args[0].dat->size[1];
    int zdim0_set_zero_kernel_MD5 = args[0].dat->size[2];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int multi_d0 = getDatDimFromOpsArg(&args[0]);
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, multi_d0);
    double * __restrict__ farray_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 1);
    ops_halo_exchanges(args, 1, range);
    ops_H_D_exchanges_host(args, 1);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[559].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

#ifdef OPS_SOA
                 ACC<double> farray(5,xdim0_set_zero_kernel_MD5, ydim0_set_zero_kernel_MD5, zdim0_set_zero_kernel_MD5, farray_p + (n_x * 1) + (n_y * xdim0_set_zero_kernel_MD5 * 1) + (n_z * xdim0_set_zero_kernel_MD5 * ydim0_set_zero_kernel_MD5 * 1));
#else
                 ACC<double> farray(5,xdim0_set_zero_kernel_MD5, ydim0_set_zero_kernel_MD5, zdim0_set_zero_kernel_MD5, farray_p + 5 * ((n_x * 1) + (n_y * xdim0_set_zero_kernel_MD5 * 1) + (n_z * xdim0_set_zero_kernel_MD5 * ydim0_set_zero_kernel_MD5 * 1)));
#endif

    int icp;

    for (icp = 1; icp <= nctmax; ++icp) {
        farray(icp-1, 0, 0, 0) = 0.0;
    }

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[559].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 1);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[559].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[559].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
    }
}

#ifdef OPS_LAZY
extern "C"
void set_zero_kernel_MD5_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("set_zero_kernel_MD5", args, 1, 559, dim, 0, range, block, set_zero_kernel_MD5_execute);
}
#endif

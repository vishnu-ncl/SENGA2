// Auto-generated at 2026-04-28 18:43:30.577751 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bounds_kernel_eqT_xr_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bounds_kernel_eqT_xr_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 4, range, 378)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 378, "bounds_kernel_eqT_xr");
        block->instance->OPS_kernels[378].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "bounds_kernel_eqT_xr");
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
    int xdim0_bounds_kernel_eqT_xr = args[0].dat->size[0];
    int ydim0_bounds_kernel_eqT_xr = args[0].dat->size[1];
    int xdim1_bounds_kernel_eqT_xr = args[1].dat->size[0];
    int ydim1_bounds_kernel_eqT_xr = args[1].dat->size[1];
    int xdim2_bounds_kernel_eqT_xr = args[2].dat->size[0];
    int ydim2_bounds_kernel_eqT_xr = args[2].dat->size[1];
    int xdim3_bounds_kernel_eqT_xr = args[3].dat->size[0];
    int ydim3_bounds_kernel_eqT_xr = args[3].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ bcl2xr_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ bcl1xr_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ ova2xr_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

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
        block->instance->OPS_kernels[378].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                 ACC<double> drhs(xdim0_bounds_kernel_eqT_xr, ydim0_bounds_kernel_eqT_xr, drhs_p + (n_x * 1) + (n_y * xdim0_bounds_kernel_eqT_xr * 1) + (n_z * xdim0_bounds_kernel_eqT_xr * ydim0_bounds_kernel_eqT_xr * 1));

                const  ACC<double> bcl2xr(xdim1_bounds_kernel_eqT_xr, ydim1_bounds_kernel_eqT_xr, bcl2xr_p + (n_x * 0) + (n_y * xdim1_bounds_kernel_eqT_xr * 1) + (n_z * xdim1_bounds_kernel_eqT_xr * ydim1_bounds_kernel_eqT_xr * 1));

                const  ACC<double> bcl1xr(xdim2_bounds_kernel_eqT_xr, ydim2_bounds_kernel_eqT_xr, bcl1xr_p + (n_x * 0) + (n_y * xdim2_bounds_kernel_eqT_xr * 1) + (n_z * xdim2_bounds_kernel_eqT_xr * ydim2_bounds_kernel_eqT_xr * 1));

                const  ACC<double> ova2xr(xdim3_bounds_kernel_eqT_xr, ydim3_bounds_kernel_eqT_xr, ova2xr_p + (n_x * 0) + (n_y * xdim3_bounds_kernel_eqT_xr * 1) + (n_z * xdim3_bounds_kernel_eqT_xr * ydim3_bounds_kernel_eqT_xr * 1));

    drhs(0, 0, 0) = drhs(0, 0, 0) - bcl2xr(0, 0, 0) - bcl1xr(0, 0, 0) * ova2xr(0, 0, 0);

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[378].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 4);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[378].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[378].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[378].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[378].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[378].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
    }
}

#ifdef OPS_LAZY
extern "C"
void bounds_kernel_eqT_xr_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bounds_kernel_eqT_xr", args, 4, 378, dim, 0, range, block, bounds_kernel_eqT_xr_execute);
}
#endif

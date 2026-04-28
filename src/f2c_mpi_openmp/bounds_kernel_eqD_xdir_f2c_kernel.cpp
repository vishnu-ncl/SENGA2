// Auto-generated at 2026-04-28 18:43:28.738978 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bounds_kernel_eqD_xdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bounds_kernel_eqD_xdir_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[3];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
#endif

//  ======
//  Timing
//  ======
    double __t1, __t2, __c1, __c2;

    ops_arg arg0 = args[0];
    ops_arg arg1 = args[1];
    ops_arg arg2 = args[2];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 3, range, 340)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 340, "bounds_kernel_eqD_xdir");
        block->instance->OPS_kernels[340].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "bounds_kernel_eqD_xdir");
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
    int xdim0_bounds_kernel_eqD_xdir = args[0].dat->size[0];
    int ydim0_bounds_kernel_eqD_xdir = args[0].dat->size[1];
    int xdim1_bounds_kernel_eqD_xdir = args[1].dat->size[0];
    int ydim1_bounds_kernel_eqD_xdir = args[1].dat->size[1];
    int xdim2_bounds_kernel_eqD_xdir = args[2].dat->size[0];
    int ydim2_bounds_kernel_eqD_xdir = args[2].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ sorpx_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ strhx_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ ratex_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 3);
    ops_halo_exchanges(args, 3, range);
    ops_H_D_exchanges_host(args, 3);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[340].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                 ACC<double> sorpx(xdim0_bounds_kernel_eqD_xdir, ydim0_bounds_kernel_eqD_xdir, sorpx_p + (n_x * 0) + (n_y * xdim0_bounds_kernel_eqD_xdir * 1) + (n_z * xdim0_bounds_kernel_eqD_xdir * ydim0_bounds_kernel_eqD_xdir * 1));

                const  ACC<double> strhx(xdim1_bounds_kernel_eqD_xdir, ydim1_bounds_kernel_eqD_xdir, strhx_p + (n_x * 0) + (n_y * xdim1_bounds_kernel_eqD_xdir * 1) + (n_z * xdim1_bounds_kernel_eqD_xdir * ydim1_bounds_kernel_eqD_xdir * 1));

                const  ACC<double> ratex(xdim2_bounds_kernel_eqD_xdir, ydim2_bounds_kernel_eqD_xdir, ratex_p + (n_x * 0) + (n_y * xdim2_bounds_kernel_eqD_xdir * 1) + (n_z * xdim2_bounds_kernel_eqD_xdir * ydim2_bounds_kernel_eqD_xdir * 1));

    sorpx(0, 0, 0) = sorpx(0, 0, 0) + strhx(0, 0, 0) * ratex(0, 0, 0);

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[340].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 3);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[340].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[340].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[340].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[340].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
    }
}

#ifdef OPS_LAZY
extern "C"
void bounds_kernel_eqD_xdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bounds_kernel_eqD_xdir", args, 3, 340, dim, 0, range, block, bounds_kernel_eqD_xdir_execute);
}
#endif

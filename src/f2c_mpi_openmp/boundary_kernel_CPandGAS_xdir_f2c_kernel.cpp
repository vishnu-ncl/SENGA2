// Auto-generated at 2026-04-28 18:43:24.768070 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void boundary_kernel_CPandGAS_xdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void boundary_kernel_CPandGAS_xdir_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[5];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
    args[4] = desc->args[4];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 5, range, 193)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 193, "boundary_kernel_CPandGAS_xdir");
        block->instance->OPS_kernels[193].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "boundary_kernel_CPandGAS_xdir");
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
    int xdim0_boundary_kernel_CPandGAS_xdir = args[0].dat->size[0];
    int ydim0_boundary_kernel_CPandGAS_xdir = args[0].dat->size[1];
    int xdim1_boundary_kernel_CPandGAS_xdir = args[1].dat->size[0];
    int ydim1_boundary_kernel_CPandGAS_xdir = args[1].dat->size[1];
    int xdim2_boundary_kernel_CPandGAS_xdir = args[2].dat->size[0];
    int ydim2_boundary_kernel_CPandGAS_xdir = args[2].dat->size[1];
    int xdim3_boundary_kernel_CPandGAS_xdir = args[3].dat->size[0];
    int ydim3_boundary_kernel_CPandGAS_xdir = args[3].dat->size[1];
    int xdim4_boundary_kernel_CPandGAS_xdir = args[4].dat->size[0];
    int ydim4_boundary_kernel_CPandGAS_xdir = args[4].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ transp_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ store7_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ strgx_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ strrx_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 5);
    ops_halo_exchanges(args, 5, range);
    ops_H_D_exchanges_host(args, 5);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[193].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                const  ACC<double> transp(xdim0_boundary_kernel_CPandGAS_xdir, ydim0_boundary_kernel_CPandGAS_xdir, transp_p + (n_x * 1) + (n_y * xdim0_boundary_kernel_CPandGAS_xdir * 1) + (n_z * xdim0_boundary_kernel_CPandGAS_xdir * ydim0_boundary_kernel_CPandGAS_xdir * 1));

                const  ACC<double> store7(xdim1_boundary_kernel_CPandGAS_xdir, ydim1_boundary_kernel_CPandGAS_xdir, store7_p + (n_x * 1) + (n_y * xdim1_boundary_kernel_CPandGAS_xdir * 1) + (n_z * xdim1_boundary_kernel_CPandGAS_xdir * ydim1_boundary_kernel_CPandGAS_xdir * 1));

                const  ACC<double> drhs(xdim2_boundary_kernel_CPandGAS_xdir, ydim2_boundary_kernel_CPandGAS_xdir, drhs_p + (n_x * 1) + (n_y * xdim2_boundary_kernel_CPandGAS_xdir * 1) + (n_z * xdim2_boundary_kernel_CPandGAS_xdir * ydim2_boundary_kernel_CPandGAS_xdir * 1));

                 ACC<double> strgx(xdim3_boundary_kernel_CPandGAS_xdir, ydim3_boundary_kernel_CPandGAS_xdir, strgx_p + (n_x * 0) + (n_y * xdim3_boundary_kernel_CPandGAS_xdir * 1) + (n_z * xdim3_boundary_kernel_CPandGAS_xdir * ydim3_boundary_kernel_CPandGAS_xdir * 1));

                 ACC<double> strrx(xdim4_boundary_kernel_CPandGAS_xdir, ydim4_boundary_kernel_CPandGAS_xdir, strrx_p + (n_x * 0) + (n_y * xdim4_boundary_kernel_CPandGAS_xdir * 1) + (n_z * xdim4_boundary_kernel_CPandGAS_xdir * ydim4_boundary_kernel_CPandGAS_xdir * 1));

    strgx(0, 0, 0) = transp(0, 0, 0);
    strrx(0, 0, 0) = store7(0, 0, 0) / drhs(0, 0, 0);

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[193].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 5);
    ops_set_halo_dirtybit3(&args[3], range);
    ops_set_halo_dirtybit3(&args[4], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[193].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[193].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[193].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[193].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[193].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[193].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
    }
}

#ifdef OPS_LAZY
extern "C"
void boundary_kernel_CPandGAS_xdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("boundary_kernel_CPandGAS_xdir", args, 5, 193, dim, 0, range, block, boundary_kernel_CPandGAS_xdir_execute);
}
#endif

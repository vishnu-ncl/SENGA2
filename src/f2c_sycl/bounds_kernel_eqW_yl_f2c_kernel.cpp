// Auto-generated at 2026-04-28 18:44:53.380012 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bounds_kernel_eqW_yl_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bounds_kernel_eqW_yl_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[8];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
    args[4] = desc->args[4];
    args[5] = desc->args[5];
    args[6] = desc->args[6];
    args[7] = desc->args[7];
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
    ops_arg arg7 = args[7];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 8, range, 415)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 415, "bounds_kernel_eqW_yl");
        block->instance->OPS_kernels[415].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "bounds_kernel_eqW_yl");
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
    int xdim0_bounds_kernel_eqW_yl = args[0].dat->size[0];
    int ydim0_bounds_kernel_eqW_yl = args[0].dat->size[1];
    int xdim1_bounds_kernel_eqW_yl = args[1].dat->size[0];
    int ydim1_bounds_kernel_eqW_yl = args[1].dat->size[1];
    int xdim2_bounds_kernel_eqW_yl = args[2].dat->size[0];
    int ydim2_bounds_kernel_eqW_yl = args[2].dat->size[1];
    int xdim3_bounds_kernel_eqW_yl = args[3].dat->size[0];
    int ydim3_bounds_kernel_eqW_yl = args[3].dat->size[1];
    int xdim4_bounds_kernel_eqW_yl = args[4].dat->size[0];
    int ydim4_bounds_kernel_eqW_yl = args[4].dat->size[1];
    int xdim5_bounds_kernel_eqW_yl = args[5].dat->size[0];
    int ydim5_bounds_kernel_eqW_yl = args[5].dat->size[1];
    int xdim6_bounds_kernel_eqW_yl = args[6].dat->size[0];
    int ydim6_bounds_kernel_eqW_yl = args[6].dat->size[1];
    int xdim7_bounds_kernel_eqW_yl = args[7].dat->size[0];
    int ydim7_bounds_kernel_eqW_yl = args[7].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ urhs_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ vrhs_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ wrhs_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ bcl2yl_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ bcl5yl_p = (double *)(args[6].data_d) + base6 - 1; // Subtracting 1 to convert to C-style

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ ova2yl_p = (double *)(args[7].data_d) + base7 - 1; // Subtracting 1 to convert to C-style

//  Subtracting 1 here as start_indx and end_indx is in Fortran style - converting it to c-style range
    int start_0 = start_indx[0]-1;
    int end_0 = end_indx[0];
    int start_1 = start_indx[1]-1;
    int end_1 = end_indx[1];
    int start_2 = start_indx[2]-1;
    int end_2 = end_indx[2];

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 8);
    ops_halo_exchanges(args, 8, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[415].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            cgh.parallel_for<class bounds_kernel_eqW_yl_sycl>(
                cl::sycl::nd_range<3>(
                    cl::sycl::range<3>(
                        ((end_2-start_2-1)/block->instance->OPS_block_size_z+1)*block->instance->OPS_block_size_z,
                        ((end_1-start_1-1)/block->instance->OPS_block_size_y+1)*block->instance->OPS_block_size_y,
                        ((end_0-start_0-1)/block->instance->OPS_block_size_x+1)*block->instance->OPS_block_size_x),
                    cl::sycl::range<3>(
                        block->instance->OPS_block_size_z,
                        block->instance->OPS_block_size_y,
                        block->instance->OPS_block_size_x)
            )
            , [=](cl::sycl::nd_item<3> item
            ) [[intel::kernel_args_restrict]] {

                int n_z = item.get_global_id()[0];
                int n_y = item.get_global_id()[1];
                int n_x = item.get_global_id()[2];

                 ACC<double> drhs(xdim0_bounds_kernel_eqW_yl, ydim0_bounds_kernel_eqW_yl, drhs_p + (n_x * 1) + (n_y * xdim0_bounds_kernel_eqW_yl * 1) + (n_z * xdim0_bounds_kernel_eqW_yl * ydim0_bounds_kernel_eqW_yl * 1));
                 ACC<double> urhs(xdim1_bounds_kernel_eqW_yl, ydim1_bounds_kernel_eqW_yl, urhs_p + (n_x * 1) + (n_y * xdim1_bounds_kernel_eqW_yl * 1) + (n_z * xdim1_bounds_kernel_eqW_yl * ydim1_bounds_kernel_eqW_yl * 1));
                 ACC<double> vrhs(xdim2_bounds_kernel_eqW_yl, ydim2_bounds_kernel_eqW_yl, vrhs_p + (n_x * 1) + (n_y * xdim2_bounds_kernel_eqW_yl * 1) + (n_z * xdim2_bounds_kernel_eqW_yl * ydim2_bounds_kernel_eqW_yl * 1));
                 ACC<double> wrhs(xdim3_bounds_kernel_eqW_yl, ydim3_bounds_kernel_eqW_yl, wrhs_p + (n_x * 1) + (n_y * xdim3_bounds_kernel_eqW_yl * 1) + (n_z * xdim3_bounds_kernel_eqW_yl * ydim3_bounds_kernel_eqW_yl * 1));
                 ACC<double> erhs(xdim4_bounds_kernel_eqW_yl, ydim4_bounds_kernel_eqW_yl, erhs_p + (n_x * 1) + (n_y * xdim4_bounds_kernel_eqW_yl * 1) + (n_z * xdim4_bounds_kernel_eqW_yl * ydim4_bounds_kernel_eqW_yl * 1));
                const  ACC<double> bcl2yl(xdim5_bounds_kernel_eqW_yl, ydim5_bounds_kernel_eqW_yl, bcl2yl_p + (n_x * 1) + (n_y * xdim5_bounds_kernel_eqW_yl * 0) + (n_z * xdim5_bounds_kernel_eqW_yl * ydim5_bounds_kernel_eqW_yl * 1));
                const  ACC<double> bcl5yl(xdim6_bounds_kernel_eqW_yl, ydim6_bounds_kernel_eqW_yl, bcl5yl_p + (n_x * 1) + (n_y * xdim6_bounds_kernel_eqW_yl * 0) + (n_z * xdim6_bounds_kernel_eqW_yl * ydim6_bounds_kernel_eqW_yl * 1));
                const  ACC<double> ova2yl(xdim7_bounds_kernel_eqW_yl, ydim7_bounds_kernel_eqW_yl, ova2yl_p + (n_x * 1) + (n_y * xdim7_bounds_kernel_eqW_yl * 0) + (n_z * xdim7_bounds_kernel_eqW_yl * ydim7_bounds_kernel_eqW_yl * 1));

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    drhs(0, 0, 0) = drhs(0, 0, 0) - bcl2yl(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0);
    urhs(0, 0, 0) = 0.0;
    vrhs(0, 0, 0) = 0.0;
    wrhs(0, 0, 0) = 0.0;
    erhs(0, 0, 0) = 0.0;

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[415].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 8);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
    ops_set_halo_dirtybit3(&args[3], range);
    ops_set_halo_dirtybit3(&args[4], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[415].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[415].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[415].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[415].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[415].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[415].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[415].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[415].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[415].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
    }
}

#ifdef OPS_LAZY
extern "C"
void bounds_kernel_eqW_yl_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bounds_kernel_eqW_yl", args, 8, 415, dim, 1, range, block, bounds_kernel_eqW_yl_execute);
}
#endif

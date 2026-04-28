// Auto-generated at 2026-04-28 18:44:54.800876 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bounds_kernel_eqAA_zdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bounds_kernel_eqAA_zdir_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[10];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
    args[4] = desc->args[4];
    args[5] = desc->args[5];
    args[6] = desc->args[6];
    args[7] = desc->args[7];
    args[8] = desc->args[8];
    args[9] = desc->args[9];
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
    ops_arg arg8 = args[8];
    ops_arg arg9 = args[9];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 10, range, 444)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 444, "bounds_kernel_eqAA_zdir");
        block->instance->OPS_kernels[444].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "bounds_kernel_eqAA_zdir");
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
    int xdim0_bounds_kernel_eqAA_zdir = args[0].dat->size[0];
    int ydim0_bounds_kernel_eqAA_zdir = args[0].dat->size[1];
    int xdim1_bounds_kernel_eqAA_zdir = args[1].dat->size[0];
    int ydim1_bounds_kernel_eqAA_zdir = args[1].dat->size[1];
    int xdim2_bounds_kernel_eqAA_zdir = args[2].dat->size[0];
    int ydim2_bounds_kernel_eqAA_zdir = args[2].dat->size[1];
    int xdim3_bounds_kernel_eqAA_zdir = args[3].dat->size[0];
    int ydim3_bounds_kernel_eqAA_zdir = args[3].dat->size[1];
    int xdim4_bounds_kernel_eqAA_zdir = args[4].dat->size[0];
    int ydim4_bounds_kernel_eqAA_zdir = args[4].dat->size[1];
    int xdim5_bounds_kernel_eqAA_zdir = args[5].dat->size[0];
    int ydim5_bounds_kernel_eqAA_zdir = args[5].dat->size[1];
    int xdim6_bounds_kernel_eqAA_zdir = args[6].dat->size[0];
    int ydim6_bounds_kernel_eqAA_zdir = args[6].dat->size[1];
    int xdim7_bounds_kernel_eqAA_zdir = args[7].dat->size[0];
    int ydim7_bounds_kernel_eqAA_zdir = args[7].dat->size[1];
    int xdim8_bounds_kernel_eqAA_zdir = args[8].dat->size[0];
    int ydim8_bounds_kernel_eqAA_zdir = args[8].dat->size[1];
    int xdim9_bounds_kernel_eqAA_zdir = args[9].dat->size[0];
    int ydim9_bounds_kernel_eqAA_zdir = args[9].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ tt1z_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ tt2z_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ tt5z_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ t1bz_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ t2bz_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ t51bz_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ t52bz_p = (double *)(args[6].data_d) + base6 - 1; // Subtracting 1 to convert to C-style

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ gam1z_p = (double *)(args[7].data_d) + base7 - 1; // Subtracting 1 to convert to C-style

    int base8 = getDatBaseFromOpsArg3D(&args[8], start_indx, 1);
    double * __restrict__ strdz_p = (double *)(args[8].data_d) + base8 - 1; // Subtracting 1 to convert to C-style

    int base9 = getDatBaseFromOpsArg3D(&args[9], start_indx, 1);
    double * __restrict__ acouz_p = (double *)(args[9].data_d) + base9 - 1; // Subtracting 1 to convert to C-style

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
    ops_H_D_exchanges_device(args, 10);
    ops_halo_exchanges(args, 10, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[444].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            cgh.parallel_for<class bounds_kernel_eqAA_zdir_sycl>(
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

                 ACC<double> tt1z(xdim0_bounds_kernel_eqAA_zdir, ydim0_bounds_kernel_eqAA_zdir, tt1z_p + (n_x * 1) + (n_y * xdim0_bounds_kernel_eqAA_zdir * 1) + (n_z * xdim0_bounds_kernel_eqAA_zdir * ydim0_bounds_kernel_eqAA_zdir * 0));
                 ACC<double> tt2z(xdim1_bounds_kernel_eqAA_zdir, ydim1_bounds_kernel_eqAA_zdir, tt2z_p + (n_x * 1) + (n_y * xdim1_bounds_kernel_eqAA_zdir * 1) + (n_z * xdim1_bounds_kernel_eqAA_zdir * ydim1_bounds_kernel_eqAA_zdir * 0));
                 ACC<double> tt5z(xdim2_bounds_kernel_eqAA_zdir, ydim2_bounds_kernel_eqAA_zdir, tt5z_p + (n_x * 1) + (n_y * xdim2_bounds_kernel_eqAA_zdir * 1) + (n_z * xdim2_bounds_kernel_eqAA_zdir * ydim2_bounds_kernel_eqAA_zdir * 0));
                const  ACC<double> t1bz(xdim3_bounds_kernel_eqAA_zdir, ydim3_bounds_kernel_eqAA_zdir, t1bz_p + (n_x * 1) + (n_y * xdim3_bounds_kernel_eqAA_zdir * 1) + (n_z * xdim3_bounds_kernel_eqAA_zdir * ydim3_bounds_kernel_eqAA_zdir * 0));
                const  ACC<double> t2bz(xdim4_bounds_kernel_eqAA_zdir, ydim4_bounds_kernel_eqAA_zdir, t2bz_p + (n_x * 1) + (n_y * xdim4_bounds_kernel_eqAA_zdir * 1) + (n_z * xdim4_bounds_kernel_eqAA_zdir * ydim4_bounds_kernel_eqAA_zdir * 0));
                const  ACC<double> t51bz(xdim5_bounds_kernel_eqAA_zdir, ydim5_bounds_kernel_eqAA_zdir, t51bz_p + (n_x * 1) + (n_y * xdim5_bounds_kernel_eqAA_zdir * 1) + (n_z * xdim5_bounds_kernel_eqAA_zdir * ydim5_bounds_kernel_eqAA_zdir * 0));
                const  ACC<double> t52bz(xdim6_bounds_kernel_eqAA_zdir, ydim6_bounds_kernel_eqAA_zdir, t52bz_p + (n_x * 1) + (n_y * xdim6_bounds_kernel_eqAA_zdir * 1) + (n_z * xdim6_bounds_kernel_eqAA_zdir * ydim6_bounds_kernel_eqAA_zdir * 0));
                const  ACC<double> gam1z(xdim7_bounds_kernel_eqAA_zdir, ydim7_bounds_kernel_eqAA_zdir, gam1z_p + (n_x * 1) + (n_y * xdim7_bounds_kernel_eqAA_zdir * 1) + (n_z * xdim7_bounds_kernel_eqAA_zdir * ydim7_bounds_kernel_eqAA_zdir * 0));
                const  ACC<double> strdz(xdim8_bounds_kernel_eqAA_zdir, ydim8_bounds_kernel_eqAA_zdir, strdz_p + (n_x * 1) + (n_y * xdim8_bounds_kernel_eqAA_zdir * 1) + (n_z * xdim8_bounds_kernel_eqAA_zdir * ydim8_bounds_kernel_eqAA_zdir * 0));
                const  ACC<double> acouz(xdim9_bounds_kernel_eqAA_zdir, ydim9_bounds_kernel_eqAA_zdir, acouz_p + (n_x * 1) + (n_y * xdim9_bounds_kernel_eqAA_zdir * 1) + (n_z * xdim9_bounds_kernel_eqAA_zdir * ydim9_bounds_kernel_eqAA_zdir * 0));

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    tt1z(0, 0, 0) = t51bz(0, 0, 0) + t52bz(0, 0, 0) * (gam1z(0, 0, 0) + 1.0) - strdz(0, 0, 0) * acouz(0, 0, 0) * t2bz(0, 0, 0);
    tt2z(0, 0, 0) = acouz(0, 0, 0) * acouz(0, 0, 0) * t1bz(0, 0, 0) - t51bz(0, 0, 0) - (gam1z(0, 0, 0) + 1.0) * t52bz(0, 0, 0);
    tt5z(0, 0, 0) = t51bz(0, 0, 0) + t52bz(0, 0, 0) * (gam1z(0, 0, 0) + 1.0) + strdz(0, 0, 0) * acouz(0, 0, 0) * t2bz(0, 0, 0);

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[444].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 10);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[444].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[444].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[444].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[444].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[444].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[444].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[444].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[444].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[444].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
        block->instance->OPS_kernels[444].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg8);
        block->instance->OPS_kernels[444].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg9);
    }
}

#ifdef OPS_LAZY
extern "C"
void bounds_kernel_eqAA_zdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bounds_kernel_eqAA_zdir", args, 10, 444, dim, 1, range, block, bounds_kernel_eqAA_zdir_execute);
}
#endif

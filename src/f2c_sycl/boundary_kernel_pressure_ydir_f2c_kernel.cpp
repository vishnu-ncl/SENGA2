// Auto-generated at 2026-04-28 18:44:48.273743 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void boundary_kernel_pressure_ydir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void boundary_kernel_pressure_ydir_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[15];
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
    args[10] = desc->args[10];
    args[11] = desc->args[11];
    args[12] = desc->args[12];
    args[13] = desc->args[13];
    args[14] = desc->args[14];
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
    ops_arg arg10 = args[10];
    ops_arg arg11 = args[11];
    ops_arg arg12 = args[12];
    ops_arg arg13 = args[13];
    ops_arg arg14 = args[14];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 15, range, 289)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 289, "boundary_kernel_pressure_ydir");
        block->instance->OPS_kernels[289].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "boundary_kernel_pressure_ydir");
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
    int xdim0_boundary_kernel_pressure_ydir = args[0].dat->size[0];
    int ydim0_boundary_kernel_pressure_ydir = args[0].dat->size[1];
    int xdim1_boundary_kernel_pressure_ydir = args[1].dat->size[0];
    int ydim1_boundary_kernel_pressure_ydir = args[1].dat->size[1];
    int xdim2_boundary_kernel_pressure_ydir = args[2].dat->size[0];
    int ydim2_boundary_kernel_pressure_ydir = args[2].dat->size[1];
    int xdim3_boundary_kernel_pressure_ydir = args[3].dat->size[0];
    int ydim3_boundary_kernel_pressure_ydir = args[3].dat->size[1];
    int xdim4_boundary_kernel_pressure_ydir = args[4].dat->size[0];
    int ydim4_boundary_kernel_pressure_ydir = args[4].dat->size[1];
    int xdim5_boundary_kernel_pressure_ydir = args[5].dat->size[0];
    int ydim5_boundary_kernel_pressure_ydir = args[5].dat->size[1];
    int xdim6_boundary_kernel_pressure_ydir = args[6].dat->size[0];
    int ydim6_boundary_kernel_pressure_ydir = args[6].dat->size[1];
    int xdim7_boundary_kernel_pressure_ydir = args[7].dat->size[0];
    int ydim7_boundary_kernel_pressure_ydir = args[7].dat->size[1];
    int xdim8_boundary_kernel_pressure_ydir = args[8].dat->size[0];
    int ydim8_boundary_kernel_pressure_ydir = args[8].dat->size[1];
    int xdim9_boundary_kernel_pressure_ydir = args[9].dat->size[0];
    int ydim9_boundary_kernel_pressure_ydir = args[9].dat->size[1];
    int xdim10_boundary_kernel_pressure_ydir = args[10].dat->size[0];
    int ydim10_boundary_kernel_pressure_ydir = args[10].dat->size[1];
    int xdim11_boundary_kernel_pressure_ydir = args[11].dat->size[0];
    int ydim11_boundary_kernel_pressure_ydir = args[11].dat->size[1];
    int xdim12_boundary_kernel_pressure_ydir = args[12].dat->size[0];
    int ydim12_boundary_kernel_pressure_ydir = args[12].dat->size[1];
    int xdim13_boundary_kernel_pressure_ydir = args[13].dat->size[0];
    int ydim13_boundary_kernel_pressure_ydir = args[13].dat->size[1];
    int xdim14_boundary_kernel_pressure_ydir = args[14].dat->size[0];
    int ydim14_boundary_kernel_pressure_ydir = args[14].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ strpy_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ bcl5y_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ t3by_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ t4by_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ t51by_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ t52by_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[6].data_d) + base6 - 1; // Subtracting 1 to convert to C-style

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ utmp_p = (double *)(args[7].data_d) + base7 - 1; // Subtracting 1 to convert to C-style

    int base8 = getDatBaseFromOpsArg3D(&args[8], start_indx, 1);
    double * __restrict__ wtmp_p = (double *)(args[8].data_d) + base8 - 1; // Subtracting 1 to convert to C-style

    int base9 = getDatBaseFromOpsArg3D(&args[9], start_indx, 1);
    double * __restrict__ prun_p = (double *)(args[9].data_d) + base9 - 1; // Subtracting 1 to convert to C-style

    int base10 = getDatBaseFromOpsArg3D(&args[10], start_indx, 1);
    double * __restrict__ store1_p = (double *)(args[10].data_d) + base10 - 1; // Subtracting 1 to convert to C-style

    int base11 = getDatBaseFromOpsArg3D(&args[11], start_indx, 1);
    double * __restrict__ store3_p = (double *)(args[11].data_d) + base11 - 1; // Subtracting 1 to convert to C-style

    int base12 = getDatBaseFromOpsArg3D(&args[12], start_indx, 1);
    double * __restrict__ store4_p = (double *)(args[12].data_d) + base12 - 1; // Subtracting 1 to convert to C-style

    int base13 = getDatBaseFromOpsArg3D(&args[13], start_indx, 1);
    double * __restrict__ store5_p = (double *)(args[13].data_d) + base13 - 1; // Subtracting 1 to convert to C-style

    int base14 = getDatBaseFromOpsArg3D(&args[14], start_indx, 1);
    double * __restrict__ store6_p = (double *)(args[14].data_d) + base14 - 1; // Subtracting 1 to convert to C-style

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
    ops_H_D_exchanges_device(args, 15);
    ops_halo_exchanges(args, 15, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[289].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            cgh.parallel_for<class boundary_kernel_pressure_ydir_sycl>(
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

                 ACC<double> strpy(xdim0_boundary_kernel_pressure_ydir, ydim0_boundary_kernel_pressure_ydir, strpy_p + (n_x * 1) + (n_y * xdim0_boundary_kernel_pressure_ydir * 0) + (n_z * xdim0_boundary_kernel_pressure_ydir * ydim0_boundary_kernel_pressure_ydir * 1));
                 ACC<double> bcl5y(xdim1_boundary_kernel_pressure_ydir, ydim1_boundary_kernel_pressure_ydir, bcl5y_p + (n_x * 1) + (n_y * xdim1_boundary_kernel_pressure_ydir * 0) + (n_z * xdim1_boundary_kernel_pressure_ydir * ydim1_boundary_kernel_pressure_ydir * 1));
                 ACC<double> t3by(xdim2_boundary_kernel_pressure_ydir, ydim2_boundary_kernel_pressure_ydir, t3by_p + (n_x * 1) + (n_y * xdim2_boundary_kernel_pressure_ydir * 0) + (n_z * xdim2_boundary_kernel_pressure_ydir * ydim2_boundary_kernel_pressure_ydir * 1));
                 ACC<double> t4by(xdim3_boundary_kernel_pressure_ydir, ydim3_boundary_kernel_pressure_ydir, t4by_p + (n_x * 1) + (n_y * xdim3_boundary_kernel_pressure_ydir * 0) + (n_z * xdim3_boundary_kernel_pressure_ydir * ydim3_boundary_kernel_pressure_ydir * 1));
                 ACC<double> t51by(xdim4_boundary_kernel_pressure_ydir, ydim4_boundary_kernel_pressure_ydir, t51by_p + (n_x * 1) + (n_y * xdim4_boundary_kernel_pressure_ydir * 0) + (n_z * xdim4_boundary_kernel_pressure_ydir * ydim4_boundary_kernel_pressure_ydir * 1));
                 ACC<double> t52by(xdim5_boundary_kernel_pressure_ydir, ydim5_boundary_kernel_pressure_ydir, t52by_p + (n_x * 1) + (n_y * xdim5_boundary_kernel_pressure_ydir * 0) + (n_z * xdim5_boundary_kernel_pressure_ydir * ydim5_boundary_kernel_pressure_ydir * 1));
                const  ACC<double> drhs(xdim6_boundary_kernel_pressure_ydir, ydim6_boundary_kernel_pressure_ydir, drhs_p + (n_x * 1) + (n_y * xdim6_boundary_kernel_pressure_ydir * 1) + (n_z * xdim6_boundary_kernel_pressure_ydir * ydim6_boundary_kernel_pressure_ydir * 1));
                const  ACC<double> utmp(xdim7_boundary_kernel_pressure_ydir, ydim7_boundary_kernel_pressure_ydir, utmp_p + (n_x * 1) + (n_y * xdim7_boundary_kernel_pressure_ydir * 1) + (n_z * xdim7_boundary_kernel_pressure_ydir * ydim7_boundary_kernel_pressure_ydir * 1));
                const  ACC<double> wtmp(xdim8_boundary_kernel_pressure_ydir, ydim8_boundary_kernel_pressure_ydir, wtmp_p + (n_x * 1) + (n_y * xdim8_boundary_kernel_pressure_ydir * 1) + (n_z * xdim8_boundary_kernel_pressure_ydir * ydim8_boundary_kernel_pressure_ydir * 1));
                const  ACC<double> prun(xdim9_boundary_kernel_pressure_ydir, ydim9_boundary_kernel_pressure_ydir, prun_p + (n_x * 1) + (n_y * xdim9_boundary_kernel_pressure_ydir * 1) + (n_z * xdim9_boundary_kernel_pressure_ydir * ydim9_boundary_kernel_pressure_ydir * 1));
                const  ACC<double> store1(xdim10_boundary_kernel_pressure_ydir, ydim10_boundary_kernel_pressure_ydir, store1_p + (n_x * 1) + (n_y * xdim10_boundary_kernel_pressure_ydir * 1) + (n_z * xdim10_boundary_kernel_pressure_ydir * ydim10_boundary_kernel_pressure_ydir * 1));
                const  ACC<double> store3(xdim11_boundary_kernel_pressure_ydir, ydim11_boundary_kernel_pressure_ydir, store3_p + (n_x * 1) + (n_y * xdim11_boundary_kernel_pressure_ydir * 1) + (n_z * xdim11_boundary_kernel_pressure_ydir * ydim11_boundary_kernel_pressure_ydir * 1));
                const  ACC<double> store4(xdim12_boundary_kernel_pressure_ydir, ydim12_boundary_kernel_pressure_ydir, store4_p + (n_x * 1) + (n_y * xdim12_boundary_kernel_pressure_ydir * 1) + (n_z * xdim12_boundary_kernel_pressure_ydir * ydim12_boundary_kernel_pressure_ydir * 1));
                const  ACC<double> store5(xdim13_boundary_kernel_pressure_ydir, ydim13_boundary_kernel_pressure_ydir, store5_p + (n_x * 1) + (n_y * xdim13_boundary_kernel_pressure_ydir * 1) + (n_z * xdim13_boundary_kernel_pressure_ydir * ydim13_boundary_kernel_pressure_ydir * 1));
                const  ACC<double> store6(xdim14_boundary_kernel_pressure_ydir, ydim14_boundary_kernel_pressure_ydir, store6_p + (n_x * 1) + (n_y * xdim14_boundary_kernel_pressure_ydir * 1) + (n_z * xdim14_boundary_kernel_pressure_ydir * ydim14_boundary_kernel_pressure_ydir * 1));

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    strpy(0, 0, 0) = prun(0, 0, 0);
    bcl5y(0, 0, 0) = store5(0, 0, 0);
    t4by(0, 0, 0) = t4by(0, 0, 0) - store6(0, 0, 0) / drhs(0, 0, 0);
    t3by(0, 0, 0) = t3by(0, 0, 0) - store4(0, 0, 0) / drhs(0, 0, 0);
    t51by(0, 0, 0) = -utmp(0, 0, 0) * store4(0, 0, 0) - wtmp(0, 0, 0) * store6(0, 0, 0);
    t52by(0, 0, 0) = -prun(0, 0, 0) * (store1(0, 0, 0) + store3(0, 0, 0));

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[289].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 15);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
    ops_set_halo_dirtybit3(&args[3], range);
    ops_set_halo_dirtybit3(&args[4], range);
    ops_set_halo_dirtybit3(&args[5], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[289].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg8);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg9);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg10);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg11);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg12);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg13);
        block->instance->OPS_kernels[289].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg14);
    }
}

#ifdef OPS_LAZY
extern "C"
void boundary_kernel_pressure_ydir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("boundary_kernel_pressure_ydir", args, 15, 289, dim, 1, range, block, boundary_kernel_pressure_ydir_execute);
}
#endif

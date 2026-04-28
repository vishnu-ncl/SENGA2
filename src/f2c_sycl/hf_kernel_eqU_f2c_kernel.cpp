// Auto-generated at 2026-04-28 18:44:46.861609 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void hf_kernel_eqU_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void hf_kernel_eqU_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 7, range, 248)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 248, "hf_kernel_eqU");
        block->instance->OPS_kernels[248].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "hf_kernel_eqU");
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
    int xdim0_hf_kernel_eqU = args[0].dat->size[0];
    int ydim0_hf_kernel_eqU = args[0].dat->size[1];
    int xdim1_hf_kernel_eqU = args[1].dat->size[0];
    int ydim1_hf_kernel_eqU = args[1].dat->size[1];
    int xdim2_hf_kernel_eqU = args[2].dat->size[0];
    int ydim2_hf_kernel_eqU = args[2].dat->size[1];
    int xdim3_hf_kernel_eqU = args[3].dat->size[0];
    int ydim3_hf_kernel_eqU = args[3].dat->size[1];
    int xdim4_hf_kernel_eqU = args[4].dat->size[0];
    int ydim4_hf_kernel_eqU = args[4].dat->size[1];
    int xdim5_hf_kernel_eqU = args[5].dat->size[0];
    int ydim5_hf_kernel_eqU = args[5].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ rate_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ vtmp_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ utmp_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ store7_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ arr1_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style

    double deltay_val = *(double *)args[6].data;

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
    ops_H_D_exchanges_device(args, 7);
    ops_halo_exchanges(args, 7, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[248].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            cgh.parallel_for<class hf_kernel_eqU_sycl>(
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

                 ACC<double> rate(xdim0_hf_kernel_eqU, ydim0_hf_kernel_eqU, rate_p + (n_x * 1) + (n_y * xdim0_hf_kernel_eqU * 1) + (n_z * xdim0_hf_kernel_eqU * ydim0_hf_kernel_eqU * 1));
                 ACC<double> erhs(xdim1_hf_kernel_eqU, ydim1_hf_kernel_eqU, erhs_p + (n_x * 1) + (n_y * xdim1_hf_kernel_eqU * 1) + (n_z * xdim1_hf_kernel_eqU * ydim1_hf_kernel_eqU * 1));
                 ACC<double> vtmp(xdim2_hf_kernel_eqU, ydim2_hf_kernel_eqU, vtmp_p + (n_x * 1) + (n_y * xdim2_hf_kernel_eqU * 1) + (n_z * xdim2_hf_kernel_eqU * ydim2_hf_kernel_eqU * 1));
                const  ACC<double> utmp(xdim3_hf_kernel_eqU, ydim3_hf_kernel_eqU, utmp_p + (n_x * 1) + (n_y * xdim3_hf_kernel_eqU * 1) + (n_z * xdim3_hf_kernel_eqU * ydim3_hf_kernel_eqU * 1));
                const  ACC<double> store7(xdim4_hf_kernel_eqU, ydim4_hf_kernel_eqU, store7_p + (n_x * 1) + (n_y * xdim4_hf_kernel_eqU * 1) + (n_z * xdim4_hf_kernel_eqU * ydim4_hf_kernel_eqU * 1));
                const  ACC<double> arr1(xdim5_hf_kernel_eqU, ydim5_hf_kernel_eqU, arr1_p + (n_x * 1) + (n_y * xdim5_hf_kernel_eqU * 1) + (n_z * xdim5_hf_kernel_eqU * ydim5_hf_kernel_eqU * 1));

                const double *deltay = &deltay_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double fornow;
    double fornow0;
    double fornow1;
    double fornow2;
    double fornow3;
    double fornow4;

    fornow0 = store7(0, 0, 0) * arr1(0, 0, 0);
    fornow1 = store7(0, 1, 0) * arr1(0, 1, 0);
    fornow2 = store7(0, 2, 0) * arr1(0, 2, 0);
    fornow3 = store7(0, 3, 0) * arr1(0, 3, 0);
    fornow4 = store7(0, 4, 0) * arr1(0, 4, 0);
    fornow = 0.0;
    fornow = ((-25.0 / 12.0) * fornow0 + (4.0) * fornow1 - (3.0) * fornow2 + (4.0 / 3.0) * fornow3 - (1.0 / 4.0) * fornow4) / deltay[0];
    rate(0, 0, 0) = rate(0, 0, 0) + fornow;
    vtmp(0, 0, 0) = vtmp(0, 0, 0) + fornow;
    erhs(0, 0, 0) = erhs(0, 0, 0) + fornow * utmp(0, 0, 0);

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[248].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 7);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[248].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[248].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[248].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[248].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[248].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[248].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[248].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
    }
}

#ifdef OPS_LAZY
extern "C"
void hf_kernel_eqU_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("hf_kernel_eqU", args, 7, 248, dim, 1, range, block, hf_kernel_eqU_execute);
}
#endif

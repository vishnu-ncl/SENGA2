// Auto-generated at 2026-04-28 18:44:59.185060 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void turbin_kernel_eqA_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void turbin_kernel_eqA_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 7, range, 555)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 555, "turbin_kernel_eqA");
        block->instance->OPS_kernels[555].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "turbin_kernel_eqA");
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
    int xdim0_turbin_kernel_eqA = args[0].dat->size[0];
    int ydim0_turbin_kernel_eqA = args[0].dat->size[1];
    int xdim1_turbin_kernel_eqA = args[1].dat->size[0];
    int ydim1_turbin_kernel_eqA = args[1].dat->size[1];
    int xdim2_turbin_kernel_eqA = args[2].dat->size[0];
    int ydim2_turbin_kernel_eqA = args[2].dat->size[1];
    int xdim3_turbin_kernel_eqA = args[3].dat->size[0];
    int ydim3_turbin_kernel_eqA = args[3].dat->size[1];
    int xdim4_turbin_kernel_eqA = args[4].dat->size[0];
    int ydim4_turbin_kernel_eqA = args[4].dat->size[1];
    int xdim5_turbin_kernel_eqA = args[5].dat->size[0];
    int ydim5_turbin_kernel_eqA = args[5].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ urun_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ utmp_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ vrun_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ vtmp_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ wrun_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ wtmp_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style

#ifdef OPS_MPI
    double * __restrict__ p_a6 = (double *)(((ops_reduction)args[6].data)->data + ((ops_reduction)args[6].data)->size * block->index);
#else //OPS_MPI
    double * __restrict__ p_a6 = (double *)((ops_reduction)args[6].data)->data;
#endif //OPS_MPI

//  Subtracting 1 here as start_indx and end_indx is in Fortran style - converting it to c-style range
    int start_0 = start_indx[0]-1;
    int end_0 = end_indx[0];
    int start_1 = start_indx[1]-1;
    int end_1 = end_indx[1];
    int start_2 = start_indx[2]-1;
    int end_2 = end_indx[2];

    int maxblocks = (end_0-start_0-1)/block->instance->OPS_block_size_x+1;
    maxblocks *= (end_1-start_1-1)/block->instance->OPS_block_size_y+1;
    maxblocks *= (end_2-start_2-1)/block->instance->OPS_block_size_z+1;
    int reduct_bytes = 0;
    size_t reduct_size = 0;

    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    reduct_size = MAX(reduct_size,1*sizeof(double));

    reallocReductArrays(block->instance, reduct_bytes);
    reduct_bytes = 0;

    arg6.data = block->instance->OPS_reduct_h + reduct_bytes;
    double *arg6_data_d = (double*)(block->instance->OPS_reduct_d + reduct_bytes);
    for (int b = 0; b < maxblocks; b++) {
        for (int d = 0; d < 1; d++)   ((double *)arg6.data)[d+b*1] = ZERO_double;
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));

    mvReductArraysToDevice(block->instance, reduct_bytes);

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 7);
    ops_halo_exchanges(args, 7, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[555].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            cl::sycl::accessor<char, 1, cl::sycl::access::mode::read_write, 
            cl::sycl::access::target::local> local_mem(reduct_size * cl::sycl::range<1>(block->instance->OPS_block_size_x * block->instance->OPS_block_size_y * block->instance->OPS_block_size_z), cgh);

            cgh.parallel_for<class turbin_kernel_eqA_sycl>(
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

                const  ACC<double> urun(xdim0_turbin_kernel_eqA, ydim0_turbin_kernel_eqA, urun_p + (n_x * 1) + (n_y * xdim0_turbin_kernel_eqA * 1) + (n_z * xdim0_turbin_kernel_eqA * ydim0_turbin_kernel_eqA * 1));
                const  ACC<double> utmp(xdim1_turbin_kernel_eqA, ydim1_turbin_kernel_eqA, utmp_p + (n_x * 1) + (n_y * xdim1_turbin_kernel_eqA * 1) + (n_z * xdim1_turbin_kernel_eqA * ydim1_turbin_kernel_eqA * 1));
                const  ACC<double> vrun(xdim2_turbin_kernel_eqA, ydim2_turbin_kernel_eqA, vrun_p + (n_x * 1) + (n_y * xdim2_turbin_kernel_eqA * 1) + (n_z * xdim2_turbin_kernel_eqA * ydim2_turbin_kernel_eqA * 1));
                const  ACC<double> vtmp(xdim3_turbin_kernel_eqA, ydim3_turbin_kernel_eqA, vtmp_p + (n_x * 1) + (n_y * xdim3_turbin_kernel_eqA * 1) + (n_z * xdim3_turbin_kernel_eqA * ydim3_turbin_kernel_eqA * 1));
                const  ACC<double> wrun(xdim4_turbin_kernel_eqA, ydim4_turbin_kernel_eqA, wrun_p + (n_x * 1) + (n_y * xdim4_turbin_kernel_eqA * 1) + (n_z * xdim4_turbin_kernel_eqA * ydim4_turbin_kernel_eqA * 1));
                const  ACC<double> wtmp(xdim5_turbin_kernel_eqA, ydim5_turbin_kernel_eqA, wtmp_p + (n_x * 1) + (n_y * xdim5_turbin_kernel_eqA * 1) + (n_z * xdim5_turbin_kernel_eqA * ydim5_turbin_kernel_eqA * 1));

                double tket[1];
                tket[0] = ZERO_double;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    *tket = *tket + urun(0, 0, 0) * urun(0, 0, 0) + utmp(0, 0, 0) * utmp(0, 0, 0) + vrun(0, 0, 0) * vrun(0, 0, 0) + vtmp(0, 0, 0) * vtmp(0, 0, 0) + wrun(0, 0, 0) * wrun(0, 0, 0) + wtmp(0, 0, 0) * wtmp(0, 0, 0);

                }

                int group_size = item.get_local_range(0);
                    group_size *= item.get_local_range(1);
                    group_size *= item.get_local_range(2);
                for (int d = 0; d < 1; d++) {
                    ops_reduction_sycl<OPS_INC>(arg6_data_d + d + item.get_group_linear_id()*1, tket[d], (double*)&local_mem[0], item, group_size);
                }
            });
        });
    }

//  ==============================
//  Reduction across blocks
//  ==============================
    mvReductArraysToHost(block->instance, reduct_bytes);

    for (int b = 0; b < maxblocks; b++)
        for (int d = 0; d < 1; d++)
            p_a6[d] = p_a6[d] + ((double *)arg6.data)[d+b*1];

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[555].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 7);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[555].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
    }
}

#ifdef OPS_LAZY
extern "C"
void turbin_kernel_eqA_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("turbin_kernel_eqA", args, 7, 555, dim, 1, range, block, turbin_kernel_eqA_execute);
}
#endif

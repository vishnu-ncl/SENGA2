// Auto-generated at 2026-04-28 18:44:48.995091 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void adaptt_kernel_err_eval_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void adaptt_kernel_err_eval_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 4, range, 315)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 315, "adaptt_kernel_err_eval");
        block->instance->OPS_kernels[315].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "adaptt_kernel_err_eval");
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
    int xdim0_adaptt_kernel_err_eval = args[0].dat->size[0];
    int ydim0_adaptt_kernel_err_eval = args[0].dat->size[1];
    int xdim1_adaptt_kernel_err_eval = args[1].dat->size[0];
    int ydim1_adaptt_kernel_err_eval = args[1].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ err_arr_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ run_arr_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    double ernrm_val = *(double *)args[2].data;

#ifdef OPS_MPI
    double * __restrict__ p_a3 = (double *)(((ops_reduction)args[3].data)->data + ((ops_reduction)args[3].data)->size * block->index);
#else //OPS_MPI
    double * __restrict__ p_a3 = (double *)((ops_reduction)args[3].data)->data;
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

    arg3.data = block->instance->OPS_reduct_h + reduct_bytes;
    double *arg3_data_d = (double*)(block->instance->OPS_reduct_d + reduct_bytes);
    for (int b = 0; b < maxblocks; b++) {
        for (int d = 0; d < 1; d++)   ((double *)arg3.data)[d+b*1] = -INFINITY_double;
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));

    mvReductArraysToDevice(block->instance, reduct_bytes);

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 4);
    ops_halo_exchanges(args, 4, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[315].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            cl::sycl::accessor<char, 1, cl::sycl::access::mode::read_write, 
            cl::sycl::access::target::local> local_mem(reduct_size * cl::sycl::range<1>(block->instance->OPS_block_size_x * block->instance->OPS_block_size_y * block->instance->OPS_block_size_z), cgh);

            cgh.parallel_for<class adaptt_kernel_err_eval_sycl>(
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

                const  ACC<double> err_arr(xdim0_adaptt_kernel_err_eval, ydim0_adaptt_kernel_err_eval, err_arr_p + (n_x * 1) + (n_y * xdim0_adaptt_kernel_err_eval * 1) + (n_z * xdim0_adaptt_kernel_err_eval * ydim0_adaptt_kernel_err_eval * 1));
                const  ACC<double> run_arr(xdim1_adaptt_kernel_err_eval, ydim1_adaptt_kernel_err_eval, run_arr_p + (n_x * 1) + (n_y * xdim1_adaptt_kernel_err_eval * 1) + (n_z * xdim1_adaptt_kernel_err_eval * ydim1_adaptt_kernel_err_eval * 1));

                const double *ernrm = &ernrm_val;
                double ertot[1];
                ertot[0] = -INFINITY_double;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double fornow;

    fornow = f2c::abs(err_arr(0, 0, 0)) / (f2c::abs(run_arr(0, 0, 0)) + ernrm[0]);
    *ertot = f2c::max(*ertot, fornow);

                }

                int group_size = item.get_local_range(0);
                    group_size *= item.get_local_range(1);
                    group_size *= item.get_local_range(2);
                for (int d = 0; d < 1; d++) {
                    ops_reduction_sycl<OPS_MAX>(arg3_data_d + d + item.get_group_linear_id()*1, ertot[d], (double*)&local_mem[0], item, group_size);
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
            p_a3[d] = MAX(p_a3[d], ((double *)arg3.data)[d+b*1]);

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[315].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 4);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[315].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[315].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[315].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void adaptt_kernel_err_eval_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("adaptt_kernel_err_eval", args, 4, 315, dim, 1, range, block, adaptt_kernel_err_eval_execute);
}
#endif

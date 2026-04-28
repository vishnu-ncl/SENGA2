// Auto-generated at 2026-04-28 18:44:59.033448 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void maths_kernel_eqBO_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void maths_kernel_eqBO_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[6];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
    args[4] = desc->args[4];
    args[5] = desc->args[5];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 6, range, 551)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 551, "maths_kernel_eqBO");
        block->instance->OPS_kernels[551].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "maths_kernel_eqBO");
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
    int xdim0_maths_kernel_eqBO = args[0].dat->size[0];
    int ydim0_maths_kernel_eqBO = args[0].dat->size[1];
    int xdim1_maths_kernel_eqBO = args[1].dat->size[0];
    int ydim1_maths_kernel_eqBO = args[1].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    int * __restrict__ itndex_p = (int *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ trun_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    double *arg2h = (double *)args[2].data;

    int *arg3h = (int *)args[3].data;

    int ipower_val = *(int *)args[4].data;

    int ispec_val = *(int *)args[5].data;

//  Subtracting 1 here as start_indx and end_indx is in Fortran style - converting it to c-style range
    int start_0 = start_indx[0]-1;
    int end_0 = end_indx[0];
    int start_1 = start_indx[1]-1;
    int end_1 = end_indx[1];
    int start_2 = start_indx[2]-1;
    int end_2 = end_indx[2];

    int consts_bytes = 0;

    consts_bytes += ROUND_UP(arg2.dim*sizeof(double));

    consts_bytes += ROUND_UP(arg3.dim*sizeof(int));

    reallocConstArrays(block->instance, consts_bytes);
    consts_bytes = 0;

    arg2.data = block->instance->OPS_consts_h + consts_bytes;
    double* arg2_data_d = (double*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg2.dim; d++)     ((double *)arg2.data)[d] = arg2h[d];
    consts_bytes += ROUND_UP(arg2.dim*sizeof(double));
    arg3.data = block->instance->OPS_consts_h + consts_bytes;
    int* arg3_data_d = (int*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg3.dim; d++)     ((int *)arg3.data)[d] = arg3h[d];
    consts_bytes += ROUND_UP(arg3.dim*sizeof(int));

    mvConstArraysToDevice(block->instance, consts_bytes);

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 6);
    ops_halo_exchanges(args, 6, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[551].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto ntinmx_sycl = (*ntinmx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ntbase_sycl = (*ntbase_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class maths_kernel_eqBO_sycl>(
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

                 ACC<int> itndex(xdim0_maths_kernel_eqBO, ydim0_maths_kernel_eqBO, itndex_p + (n_x * 1) + (n_y * xdim0_maths_kernel_eqBO * 1) + (n_z * xdim0_maths_kernel_eqBO * ydim0_maths_kernel_eqBO * 1));
                const  ACC<double> trun(xdim1_maths_kernel_eqBO, ydim1_maths_kernel_eqBO, trun_p + (n_x * 1) + (n_y * xdim1_maths_kernel_eqBO * 1) + (n_z * xdim1_maths_kernel_eqBO * ydim1_maths_kernel_eqBO * 1));

                const double *tinthi = arg2_data_d;
                const int *ntint = arg3_data_d;
                const int *ipower = &ipower_val;
                const int *ispec = &ispec_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    int itint;

    itint = 1;
    while (trun(0, 0, 0) > tinthi[(itint-1)+(ispec[0]-1)*ntinmx_sycl[0]] && itint < ntint[ispec[0]-1]) {
        itint = itint + 1;
    }
    itndex(0, 0, 0) = itndex(0, 0, 0) + (itint - 1) * cl::sycl::pow(ntbase_sycl[0], ipower[0]);

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[551].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 6);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[551].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[551].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[551].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void maths_kernel_eqBO_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("maths_kernel_eqBO", args, 6, 551, dim, 1, range, block, maths_kernel_eqBO_execute);
}
#endif

// Auto-generated at 2026-04-28 18:44:45.750125 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void maths_kernel_eqBD_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void maths_kernel_eqBD_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 10, range, 212)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 212, "maths_kernel_eqBD");
        block->instance->OPS_kernels[212].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "maths_kernel_eqBD");
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
    int xdim0_maths_kernel_eqBD = args[0].dat->size[0];
    int ydim0_maths_kernel_eqBD = args[0].dat->size[1];
    int xdim1_maths_kernel_eqBD = args[1].dat->size[0];
    int ydim1_maths_kernel_eqBD = args[1].dat->size[1];
    int xdim2_maths_kernel_eqBD = args[2].dat->size[0];
    int ydim2_maths_kernel_eqBD = args[2].dat->size[1];
    int xdim3_maths_kernel_eqBD = args[3].dat->size[0];
    int ydim3_maths_kernel_eqBD = args[3].dat->size[1];
    int xdim4_maths_kernel_eqBD = args[4].dat->size[0];
    int ydim4_maths_kernel_eqBD = args[4].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ combo1_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ combo2_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ combo3_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ transp_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ yrhs_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    double *arg5h = (double *)args[5].data;

    double *arg6h = (double *)args[6].data;

    int ncocon_val = *(int *)args[7].data;

    int ncocm1_val = *(int *)args[8].data;

    int ispec_val = *(int *)args[9].data;

//  Subtracting 1 here as start_indx and end_indx is in Fortran style - converting it to c-style range
    int start_0 = start_indx[0]-1;
    int end_0 = end_indx[0];
    int start_1 = start_indx[1]-1;
    int end_1 = end_indx[1];
    int start_2 = start_indx[2]-1;
    int end_2 = end_indx[2];

    int consts_bytes = 0;

    consts_bytes += ROUND_UP(arg5.dim*sizeof(double));

    consts_bytes += ROUND_UP(arg6.dim*sizeof(double));

    reallocConstArrays(block->instance, consts_bytes);
    consts_bytes = 0;

    arg5.data = block->instance->OPS_consts_h + consts_bytes;
    double* arg5_data_d = (double*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg5.dim; d++)     ((double *)arg5.data)[d] = arg5h[d];
    consts_bytes += ROUND_UP(arg5.dim*sizeof(double));
    arg6.data = block->instance->OPS_consts_h + consts_bytes;
    double* arg6_data_d = (double*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg6.dim; d++)     ((double *)arg6.data)[d] = arg6h[d];
    consts_bytes += ROUND_UP(arg6.dim*sizeof(double));

    mvConstArraysToDevice(block->instance, consts_bytes);

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 10);
    ops_halo_exchanges(args, 10, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[212].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto nccfmx_sycl = (*nccfmx_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class maths_kernel_eqBD_sycl>(
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

                 ACC<double> combo1(xdim0_maths_kernel_eqBD, ydim0_maths_kernel_eqBD, combo1_p + (n_x * 1) + (n_y * xdim0_maths_kernel_eqBD * 1) + (n_z * xdim0_maths_kernel_eqBD * ydim0_maths_kernel_eqBD * 1));
                 ACC<double> combo2(xdim1_maths_kernel_eqBD, ydim1_maths_kernel_eqBD, combo2_p + (n_x * 1) + (n_y * xdim1_maths_kernel_eqBD * 1) + (n_z * xdim1_maths_kernel_eqBD * ydim1_maths_kernel_eqBD * 1));
                 ACC<double> combo3(xdim2_maths_kernel_eqBD, ydim2_maths_kernel_eqBD, combo3_p + (n_x * 1) + (n_y * xdim2_maths_kernel_eqBD * 1) + (n_z * xdim2_maths_kernel_eqBD * ydim2_maths_kernel_eqBD * 1));
                const  ACC<double> transp(xdim3_maths_kernel_eqBD, ydim3_maths_kernel_eqBD, transp_p + (n_x * 1) + (n_y * xdim3_maths_kernel_eqBD * 1) + (n_z * xdim3_maths_kernel_eqBD * ydim3_maths_kernel_eqBD * 1));
                const  ACC<double> yrhs(xdim4_maths_kernel_eqBD, ydim4_maths_kernel_eqBD, yrhs_p + (n_x * 1) + (n_y * xdim4_maths_kernel_eqBD * 1) + (n_z * xdim4_maths_kernel_eqBD * ydim4_maths_kernel_eqBD * 1));

                const double *condco = arg5_data_d;
                const double *ovwmol = arg6_data_d;
                const int *ncocon = &ncocon_val;
                const int *ncocm1 = &ncocm1_val;
                const int *ispec = &ispec_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double fornow;
    double ctrans;
    int icp;

    fornow = condco[(ncocon[0]-1)+(ispec[0]-1)*nccfmx_sycl[0]];
    for (icp = ncocm1[0]; icp >= 1; icp -= 1) {
        fornow = fornow * transp(0, 0, 0) + condco[(icp-1)+(ispec[0]-1)*nccfmx_sycl[0]];
    }
    ctrans = cl::sycl::exp(fornow);
    fornow = yrhs(0, 0, 0) * ovwmol[ispec[0]-1];
    combo1(0, 0, 0) = combo1(0, 0, 0) + fornow * ctrans;
    combo2(0, 0, 0) = combo2(0, 0, 0) + fornow / ctrans;
    combo3(0, 0, 0) = combo3(0, 0, 0) + fornow;

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[212].time += __t1 - __t2;
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
        block->instance->OPS_kernels[212].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[212].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[212].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[212].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[212].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[212].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
    }
}

#ifdef OPS_LAZY
extern "C"
void maths_kernel_eqBD_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("maths_kernel_eqBD", args, 10, 212, dim, 1, range, block, maths_kernel_eqBD_execute);
}
#endif

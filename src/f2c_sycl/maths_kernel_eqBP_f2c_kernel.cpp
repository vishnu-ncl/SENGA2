// Auto-generated at 2026-04-28 18:44:59.080725 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void maths_kernel_eqBP_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void maths_kernel_eqBP_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[12];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 12, range, 552)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 552, "maths_kernel_eqBP");
        block->instance->OPS_kernels[552].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "maths_kernel_eqBP");
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
    int xdim0_maths_kernel_eqBP = args[0].dat->size[0];
    int ydim0_maths_kernel_eqBP = args[0].dat->size[1];
    int xdim1_maths_kernel_eqBP = args[1].dat->size[0];
    int ydim1_maths_kernel_eqBP = args[1].dat->size[1];
    int xdim2_maths_kernel_eqBP = args[2].dat->size[0];
    int ydim2_maths_kernel_eqBP = args[2].dat->size[1];
    int xdim3_maths_kernel_eqBP = args[3].dat->size[0];
    int ydim3_maths_kernel_eqBP = args[3].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ erun_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ trun_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    int * __restrict__ itndex_p = (int *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ yrun_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    double *arg4h = (double *)args[4].data;

    int *arg5h = (int *)args[5].data;

    int *arg6h = (int *)args[6].data;

    int *arg7h = (int *)args[7].data;

    int ipower_val = *(int *)args[8].data;

    int icoef1_val = *(int *)args[9].data;

    int icoef2_val = *(int *)args[10].data;

    int ispec_val = *(int *)args[11].data;

//  Subtracting 1 here as start_indx and end_indx is in Fortran style - converting it to c-style range
    int start_0 = start_indx[0]-1;
    int end_0 = end_indx[0];
    int start_1 = start_indx[1]-1;
    int end_1 = end_indx[1];
    int start_2 = start_indx[2]-1;
    int end_2 = end_indx[2];

    int consts_bytes = 0;

    consts_bytes += ROUND_UP(arg4.dim*sizeof(double));

    consts_bytes += ROUND_UP(arg5.dim*sizeof(int));

    consts_bytes += ROUND_UP(arg6.dim*sizeof(int));

    consts_bytes += ROUND_UP(arg7.dim*sizeof(int));

    reallocConstArrays(block->instance, consts_bytes);
    consts_bytes = 0;

    arg4.data = block->instance->OPS_consts_h + consts_bytes;
    double* arg4_data_d = (double*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg4.dim; d++)     ((double *)arg4.data)[d] = arg4h[d];
    consts_bytes += ROUND_UP(arg4.dim*sizeof(double));
    arg5.data = block->instance->OPS_consts_h + consts_bytes;
    int* arg5_data_d = (int*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg5.dim; d++)     ((int *)arg5.data)[d] = arg5h[d];
    consts_bytes += ROUND_UP(arg5.dim*sizeof(int));
    arg6.data = block->instance->OPS_consts_h + consts_bytes;
    int* arg6_data_d = (int*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg6.dim; d++)     ((int *)arg6.data)[d] = arg6h[d];
    consts_bytes += ROUND_UP(arg6.dim*sizeof(int));
    arg7.data = block->instance->OPS_consts_h + consts_bytes;
    int* arg7_data_d = (int*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg7.dim; d++)     ((int *)arg7.data)[d] = arg7h[d];
    consts_bytes += ROUND_UP(arg7.dim*sizeof(int));

    mvConstArraysToDevice(block->instance, consts_bytes);

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 12);
    ops_halo_exchanges(args, 12, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[552].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto ncofmx_sycl = (*ncofmx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ntinmx_sycl = (*ntinmx_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class maths_kernel_eqBP_sycl>(
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

                 ACC<double> erun(xdim0_maths_kernel_eqBP, ydim0_maths_kernel_eqBP, erun_p + (n_x * 1) + (n_y * xdim0_maths_kernel_eqBP * 1) + (n_z * xdim0_maths_kernel_eqBP * ydim0_maths_kernel_eqBP * 1));
                const  ACC<double> trun(xdim1_maths_kernel_eqBP, ydim1_maths_kernel_eqBP, trun_p + (n_x * 1) + (n_y * xdim1_maths_kernel_eqBP * 1) + (n_z * xdim1_maths_kernel_eqBP * ydim1_maths_kernel_eqBP * 1));
                const  ACC<int> itndex(xdim2_maths_kernel_eqBP, ydim2_maths_kernel_eqBP, itndex_p + (n_x * 1) + (n_y * xdim2_maths_kernel_eqBP * 1) + (n_z * xdim2_maths_kernel_eqBP * ydim2_maths_kernel_eqBP * 1));
                const  ACC<double> yrun(xdim3_maths_kernel_eqBP, ydim3_maths_kernel_eqBP, yrun_p + (n_x * 1) + (n_y * xdim3_maths_kernel_eqBP * 1) + (n_z * xdim3_maths_kernel_eqBP * ydim3_maths_kernel_eqBP * 1));

                const double *amasch = arg4_data_d;
                const int *ncpoly = arg5_data_d;
                const int *ncpom1 = arg6_data_d;
                const int *ncenth = arg7_data_d;
                const int *ipower = &ipower_val;
                const int *icoef1 = &icoef1_val;
                const int *icoef2 = &icoef2_val;
                const int *ispec = &ispec_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double fornow;
    int icp;
    int itint;

    itint = 1 + f2c::mod(itndex(0, 0, 0), icoef1[0]) / icoef2[0];
    fornow = amasch[(ncpoly[(itint-1)+(ispec[0]-1)*ntinmx_sycl[0]]-1)+(itint-1)*ncofmx_sycl[0]+(ispec[0]-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
    for (icp = ncpom1[(itint-1)+(ispec[0]-1)*ntinmx_sycl[0]]; icp >= 1; icp -= 1) {
        fornow = fornow * trun(0, 0, 0) + amasch[(icp-1)+(itint-1)*ncofmx_sycl[0]+(ispec[0]-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
    }
    fornow = amasch[(ncenth[(itint-1)+(ispec[0]-1)*ntinmx_sycl[0]]-1)+(itint-1)*ncofmx_sycl[0]+(ispec[0]-1)*ncofmx_sycl[0]*ntinmx_sycl[0]] + fornow * trun(0, 0, 0);
    erun(0, 0, 0) = erun(0, 0, 0) + fornow * yrun(0, 0, 0);

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[552].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 12);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[552].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[552].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[552].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[552].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[552].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
    }
}

#ifdef OPS_LAZY
extern "C"
void maths_kernel_eqBP_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("maths_kernel_eqBP", args, 12, 552, dim, 1, range, block, maths_kernel_eqBP_execute);
}
#endif

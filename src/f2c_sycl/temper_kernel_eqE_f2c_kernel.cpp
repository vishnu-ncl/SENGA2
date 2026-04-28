// Auto-generated at 2026-04-28 18:44:59.462963 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void temper_kernel_eqE_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void temper_kernel_eqE_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[11];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 11, range, 563)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 563, "temper_kernel_eqE");
        block->instance->OPS_kernels[563].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "temper_kernel_eqE");
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
    int xdim0_temper_kernel_eqE = args[0].dat->size[0];
    int ydim0_temper_kernel_eqE = args[0].dat->size[1];
    int xdim1_temper_kernel_eqE = args[1].dat->size[0];
    int ydim1_temper_kernel_eqE = args[1].dat->size[1];
    int xdim2_temper_kernel_eqE = args[2].dat->size[0];
    int ydim2_temper_kernel_eqE = args[2].dat->size[1];
    int xdim3_temper_kernel_eqE = args[3].dat->size[0];
    int ydim3_temper_kernel_eqE = args[3].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ transp_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    int * __restrict__ itndex_p = (int *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ yrhs_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ trun_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    double *arg4h = (double *)args[4].data;

    int *arg5h = (int *)args[5].data;

    int *arg6h = (int *)args[6].data;

    double *arg7h = (double *)args[7].data;

    int *arg8h = (int *)args[8].data;

    int ispec_val = *(int *)args[9].data;

    int ipower_val = *(int *)args[10].data;

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

    consts_bytes += ROUND_UP(arg7.dim*sizeof(double));

    consts_bytes += ROUND_UP(arg8.dim*sizeof(int));

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
    double* arg7_data_d = (double*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg7.dim; d++)     ((double *)arg7.data)[d] = arg7h[d];
    consts_bytes += ROUND_UP(arg7.dim*sizeof(double));
    arg8.data = block->instance->OPS_consts_h + consts_bytes;
    int* arg8_data_d = (int*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg8.dim; d++)     ((int *)arg8.data)[d] = arg8h[d];
    consts_bytes += ROUND_UP(arg8.dim*sizeof(int));

    mvConstArraysToDevice(block->instance, consts_bytes);

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 11);
    ops_halo_exchanges(args, 11, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[563].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto ncofmx_sycl = (*ncofmx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ntinmx_sycl = (*ntinmx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ntbase_sycl = (*ntbase_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class temper_kernel_eqE_sycl>(
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

                 ACC<double> transp(xdim0_temper_kernel_eqE, ydim0_temper_kernel_eqE, transp_p + (n_x * 1) + (n_y * xdim0_temper_kernel_eqE * 1) + (n_z * xdim0_temper_kernel_eqE * ydim0_temper_kernel_eqE * 1));
                 ACC<int> itndex(xdim1_temper_kernel_eqE, ydim1_temper_kernel_eqE, itndex_p + (n_x * 1) + (n_y * xdim1_temper_kernel_eqE * 1) + (n_z * xdim1_temper_kernel_eqE * ydim1_temper_kernel_eqE * 1));
                const  ACC<double> yrhs(xdim2_temper_kernel_eqE, ydim2_temper_kernel_eqE, yrhs_p + (n_x * 1) + (n_y * xdim2_temper_kernel_eqE * 1) + (n_z * xdim2_temper_kernel_eqE * ydim2_temper_kernel_eqE * 1));
                const  ACC<double> trun(xdim3_temper_kernel_eqE, ydim3_temper_kernel_eqE, trun_p + (n_x * 1) + (n_y * xdim3_temper_kernel_eqE * 1) + (n_z * xdim3_temper_kernel_eqE * ydim3_temper_kernel_eqE * 1));

                const double *amascp = arg4_data_d;
                const int *ncpoly = arg5_data_d;
                const int *ncpom1 = arg6_data_d;
                const double *tinthi = arg7_data_d;
                const int *ntint = arg8_data_d;
                const int *ispec = &ispec_val;
                const int *ipower = &ipower_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double cpfory;
    int itint;
    int icp;

    itint = 1;
    while (trun(0, 0, 0) > tinthi[(itint-1)+(ispec[0]-1)*ntinmx_sycl[0]] && itint < ntint[ispec[0]-1]) {
        itint = itint + 1;
    }
    itndex(0, 0, 0) = itndex(0, 0, 0) + (itint - 1) * cl::sycl::pow(ntbase_sycl[0], ipower[0]);
    cpfory = amascp[(ncpoly[(itint-1)+(ispec[0]-1)*ntinmx_sycl[0]]-1)+(itint-1)*ncofmx_sycl[0]+(ispec[0]-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
    for (icp = ncpom1[(itint-1)+(ispec[0]-1)*ntinmx_sycl[0]]; icp >= 1; icp -= 1) {
        cpfory = cpfory * trun(0, 0, 0) + amascp[(icp-1)+(itint-1)*ncofmx_sycl[0]+(ispec[0]-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
    }
    transp(0, 0, 0) = transp(0, 0, 0) + yrhs(0, 0, 0) * cpfory;

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[563].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 11);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[563].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[563].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[563].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[563].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[563].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
    }
}

#ifdef OPS_LAZY
extern "C"
void temper_kernel_eqE_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("temper_kernel_eqE", args, 11, 563, dim, 1, range, block, temper_kernel_eqE_execute);
}
#endif

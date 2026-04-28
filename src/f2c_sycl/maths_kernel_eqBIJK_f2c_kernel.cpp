// Auto-generated at 2026-04-28 18:44:47.797949 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void maths_kernel_eqBIJK_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void maths_kernel_eqBIJK_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[9];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
    args[4] = desc->args[4];
    args[5] = desc->args[5];
    args[6] = desc->args[6];
    args[7] = desc->args[7];
    args[8] = desc->args[8];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 9, range, 270)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 270, "maths_kernel_eqBIJK");
        block->instance->OPS_kernels[270].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "maths_kernel_eqBIJK");
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
    int xdim0_maths_kernel_eqBIJK = args[0].dat->size[0];
    int ydim0_maths_kernel_eqBIJK = args[0].dat->size[1];
    int xdim1_maths_kernel_eqBIJK = args[1].dat->size[0];
    int ydim1_maths_kernel_eqBIJK = args[1].dat->size[1];
    int xdim2_maths_kernel_eqBIJK = args[2].dat->size[0];
    int ydim2_maths_kernel_eqBIJK = args[2].dat->size[1];
    int zdim2_maths_kernel_eqBIJK = args[2].dat->size[2];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ difmix_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ transp_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int multi_d2 = getDatDimFromOpsArg(&args[2]);
    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, multi_d2);
    double * __restrict__ yrhs_md_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    double *arg3h = (double *)args[3].data;

    double *arg4h = (double *)args[4].data;

    double *arg5h = (double *)args[5].data;

    double *arg6h = (double *)args[6].data;

    int ncovis_val = *(int *)args[7].data;

    int ncovm1_val = *(int *)args[8].data;

//  Subtracting 1 here as start_indx and end_indx is in Fortran style - converting it to c-style range
    int start_0 = start_indx[0]-1;
    int end_0 = end_indx[0];
    int start_1 = start_indx[1]-1;
    int end_1 = end_indx[1];
    int start_2 = start_indx[2]-1;
    int end_2 = end_indx[2];

    int consts_bytes = 0;

    consts_bytes += ROUND_UP(arg3.dim*sizeof(double));

    consts_bytes += ROUND_UP(arg4.dim*sizeof(double));

    consts_bytes += ROUND_UP(arg5.dim*sizeof(double));

    consts_bytes += ROUND_UP(arg6.dim*sizeof(double));

    reallocConstArrays(block->instance, consts_bytes);
    consts_bytes = 0;

    arg3.data = block->instance->OPS_consts_h + consts_bytes;
    double* arg3_data_d = (double*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg3.dim; d++)     ((double *)arg3.data)[d] = arg3h[d];
    consts_bytes += ROUND_UP(arg3.dim*sizeof(double));
    arg4.data = block->instance->OPS_consts_h + consts_bytes;
    double* arg4_data_d = (double*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg4.dim; d++)     ((double *)arg4.data)[d] = arg4h[d];
    consts_bytes += ROUND_UP(arg4.dim*sizeof(double));
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
    ops_H_D_exchanges_device(args, 9);
    ops_halo_exchanges(args, 9, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[270].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto nspcmx_sycl = (*nspcmx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto nvcfmx_sycl = (*nvcfmx_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class maths_kernel_eqBIJK_sycl>(
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

                 ACC<double> difmix(xdim0_maths_kernel_eqBIJK, ydim0_maths_kernel_eqBIJK, difmix_p + (n_x * 1) + (n_y * xdim0_maths_kernel_eqBIJK * 1) + (n_z * xdim0_maths_kernel_eqBIJK * ydim0_maths_kernel_eqBIJK * 1));
                const  ACC<double> transp(xdim1_maths_kernel_eqBIJK, ydim1_maths_kernel_eqBIJK, transp_p + (n_x * 1) + (n_y * xdim1_maths_kernel_eqBIJK * 1) + (n_z * xdim1_maths_kernel_eqBIJK * ydim1_maths_kernel_eqBIJK * 1));
#ifdef OPS_SOA
                const  ACC<double> yrhs_md(22,xdim2_maths_kernel_eqBIJK, ydim2_maths_kernel_eqBIJK, zdim2_maths_kernel_eqBIJK, yrhs_md_p + (n_x * 1) + (n_y * xdim2_maths_kernel_eqBIJK * 1) + (n_z * xdim2_maths_kernel_eqBIJK * ydim2_maths_kernel_eqBIJK * 1));
#else
                const  ACC<double> yrhs_md(22,xdim2_maths_kernel_eqBIJK, ydim2_maths_kernel_eqBIJK, zdim2_maths_kernel_eqBIJK, yrhs_md_p + 22 * ((n_x * 1) + (n_y * xdim2_maths_kernel_eqBIJK * 1) + (n_z * xdim2_maths_kernel_eqBIJK * ydim2_maths_kernel_eqBIJK * 1)));
#endif

                const double *viscco = arg3_data_d;
                const double *wilko1 = arg4_data_d;
                const double *wilko2 = arg5_data_d;
                const double *ovwmol = arg6_data_d;
                const int *ncovis = &ncovis_val;
                const int *ncovm1 = &ncovm1_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double fornow;
    double combo1;
    double combo2;
    double ctrans[(22)];
    int ispec;
    int jspec;
    int icp;

    for (ispec = 1; ispec <= nspcmx_sycl[0]; ++ispec) {
        fornow = viscco[(ncovis[0]-1)+(ispec-1)*nvcfmx_sycl[0]];
        for (icp = ncovm1[0]; icp >= 1; icp -= 1) {
            fornow = fornow * transp(0, 0, 0) + viscco[(icp-1)+(ispec-1)*nvcfmx_sycl[0]];
        }
        ctrans[ispec-1] = cl::sycl::exp(fornow);
    }
    combo1 = 0.0;
    for (ispec = 1; ispec <= nspcmx_sycl[0]; ++ispec) {
        combo2 = 0.0;
        for (jspec = 1; jspec <= nspcmx_sycl[0]; ++jspec) {
            fornow = cl::sycl::sqrt(ctrans[ispec-1] / ctrans[jspec-1]);
            fornow = 1.0 + fornow * wilko2[(jspec-1)+(ispec-1)*nspcmx_sycl[0]];
            fornow = wilko1[(jspec-1)+(ispec-1)*nspcmx_sycl[0]] * fornow * fornow;
            combo2 = combo2 + yrhs_md(jspec-1, 0, 0, 0) * ovwmol[jspec-1] * fornow;
        }
        fornow = ctrans[ispec-1] / combo2;
        combo1 = combo1 + yrhs_md(ispec-1, 0, 0, 0) * ovwmol[ispec-1] * fornow;
    }
    difmix(0, 0, 0) = combo1;

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[270].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 9);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[270].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[270].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[270].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[270].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
    }
}

#ifdef OPS_LAZY
extern "C"
void maths_kernel_eqBIJK_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("maths_kernel_eqBIJK", args, 9, 270, dim, 1, range, block, maths_kernel_eqBIJK_execute);
}
#endif

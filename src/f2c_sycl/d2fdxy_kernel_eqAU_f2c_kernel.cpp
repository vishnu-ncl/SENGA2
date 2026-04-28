// Auto-generated at 2026-04-28 18:44:37.780702 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void d2fdxy_kernel_eqAU_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void d2fdxy_kernel_eqAU_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 4, range, 64)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 64, "d2fdxy_kernel_eqAU");
        block->instance->OPS_kernels[64].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "d2fdxy_kernel_eqAU");
#endif

//  =================================================
//  compute locally allocated range for the sub-block
//  =================================================
    int start_indx[3];
    int end_indx[3];
    int arg_idx[3];

//  Range here is in C-style while start_indx and end_indx is Fortran style
#if defined(OPS_MPI) && !defined(OPS_LAZY)
    if ( getRange(block, start_indx, end_indx, range) < 0 ) return;
#else
    for (int n = 0; n < 3; n++) {
        start_indx[n] = range[2*n] + 1;
        end_indx[n]   = range[2*n+1];
    }
#endif

#ifdef OPS_MPI
    getIdx(block, start_indx, arg_idx);
#else
    arg_idx[0] = start_indx[0];
    arg_idx[1] = start_indx[1];
    arg_idx[2] = start_indx[2];
#endif

//  ======================================================
//  Initialize global variable with the dimensions of dats
//  ======================================================
    int xdim0_d2fdxy_kernel_eqAU = args[0].dat->size[0];
    int ydim0_d2fdxy_kernel_eqAU = args[0].dat->size[1];
    int xdim1_d2fdxy_kernel_eqAU = args[1].dat->size[0];
    int ydim1_d2fdxy_kernel_eqAU = args[1].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ functn_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ fderiv_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int nyglbl_val = *(int *)args[2].data;

//  Subtracting 1 here as start_indx and end_indx is in Fortran style - converting it to c-style range
    int start_0 = start_indx[0]-1;
    int end_0 = end_indx[0];
    int arg_idx_0 = arg_idx[0];
    int start_1 = start_indx[1]-1;
    int end_1 = end_indx[1];
    int arg_idx_1 = arg_idx[1];
    int start_2 = start_indx[2]-1;
    int end_2 = end_indx[2];
    int arg_idx_2 = arg_idx[2];

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 4);
    ops_halo_exchanges(args, 4, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[64].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto acf3xy_sycl = (*acf3xy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcf3xy_sycl = (*bcf3xy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acf4xy_sycl = (*acf4xy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcf4xy_sycl = (*bcf4xy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccf4xy_sycl = (*ccf4xy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acf5xy_sycl = (*acf5xy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcf5xy_sycl = (*bcf5xy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccf5xy_sycl = (*ccf5xy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcf5xy_sycl = (*dcf5xy_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class d2fdxy_kernel_eqAU_sycl>(
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

                int idx[] = {arg_idx_0+n_x, arg_idx_1+n_y, arg_idx_2+n_z};

                const  ACC<double> functn(xdim0_d2fdxy_kernel_eqAU, ydim0_d2fdxy_kernel_eqAU, functn_p + (n_x * 1) + (n_y * xdim0_d2fdxy_kernel_eqAU * 1) + (n_z * xdim0_d2fdxy_kernel_eqAU * ydim0_d2fdxy_kernel_eqAU * 1));
                 ACC<double> fderiv(xdim1_d2fdxy_kernel_eqAU, ydim1_d2fdxy_kernel_eqAU, fderiv_p + (n_x * 1) + (n_y * xdim1_d2fdxy_kernel_eqAU * 1) + (n_z * xdim1_d2fdxy_kernel_eqAU * ydim1_d2fdxy_kernel_eqAU * 1));

                const int *nyglbl = &nyglbl_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double fdiffa;
    double fdiffb;
    double fdiffc;
    double fdiffd;
    double fstora;
    double fstorb;
    double fstorc;
    int ic;
    int jc;

    ic = idx[0];
    jc = idx[1];
    if (ic >= 3 && ic <= 5 && jc >= nyglbl[0] - 4 && jc <= nyglbl[0] - 2) {
        fdiffa = functn(1, 1, 0) - functn(1, -1, 0) - functn(-1, 1, 0) + functn(-1, -1, 0);
        fdiffb = functn(2, 2, 0) - functn(2, -2, 0) - functn(-2, 2, 0) + functn(-2, -2, 0);
        fderiv(0, 0, 0) = acf3xy_sycl[0] * fdiffa + bcf3xy_sycl[0] * fdiffb;
        fstora = fdiffa;
        fstorb = fdiffb;
    }
    if (ic >= 4 && ic <= 5 && jc >= nyglbl[0] - 4 && jc <= nyglbl[0] - 3) {
        fdiffc = functn(3, 3, 0) - functn(3, -3, 0) - functn(-3, 3, 0) + functn(-3, -3, 0);
        fderiv(0, 0, 0) = acf4xy_sycl[0] * fstora + bcf4xy_sycl[0] * fstorb + ccf4xy_sycl[0] * fdiffc;
        fstorc = fdiffc;
    }
    if (ic == 5 && jc == nyglbl[0] - 4) {
        fdiffd = functn(4, 4, 0) - functn(4, -4, 0) - functn(-4, 4, 0) + functn(-4, -4, 0);
        fderiv(0, 0, 0) = acf5xy_sycl[0] * fstora + bcf5xy_sycl[0] * fstorb + ccf5xy_sycl[0] * fstorc + dcf5xy_sycl[0] * fdiffd;
    }

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[64].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 4);
    ops_set_halo_dirtybit3(&args[1], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[64].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[64].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[64].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void d2fdxy_kernel_eqAU_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("d2fdxy_kernel_eqAU", args, 4, 64, dim, 1, range, block, d2fdxy_kernel_eqAU_execute);
}
#endif

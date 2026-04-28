// Auto-generated at 2026-04-28 18:44:34.404400 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void dfbydz_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void dfbydz_kernel_main_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 7, range, 11)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 11, "dfbydz_kernel_main");
        block->instance->OPS_kernels[11].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "dfbydz_kernel_main");
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
    int xdim0_dfbydz_kernel_main = args[0].dat->size[0];
    int ydim0_dfbydz_kernel_main = args[0].dat->size[1];
    int xdim1_dfbydz_kernel_main = args[1].dat->size[0];
    int ydim1_dfbydz_kernel_main = args[1].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ functn_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ fderiv_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int nzglbl_val = *(int *)args[2].data;

    int nendzl_val = *(int *)args[3].data;

    int nendzr_val = *(int *)args[4].data;

    int nbound_val = *(int *)args[5].data;

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
    ops_H_D_exchanges_device(args, 7);
    ops_halo_exchanges(args, 7, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[11].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto acoffz_sycl = (*acoffz_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcoffz_sycl = (*bcoffz_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccoffz_sycl = (*ccoffz_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcoffz_sycl = (*dcoffz_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ecoffz_sycl = (*ecoffz_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof1z_sycl = (*acof1z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof1z_sycl = (*bcof1z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof1z_sycl = (*ccof1z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcof1z_sycl = (*dcof1z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof2z_sycl = (*acof2z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof2z_sycl = (*bcof2z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof2z_sycl = (*ccof2z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcof2z_sycl = (*dcof2z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof3z_sycl = (*acof3z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof3z_sycl = (*bcof3z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof4z_sycl = (*acof4z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof4z_sycl = (*bcof4z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof4z_sycl = (*ccof4z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof5z_sycl = (*acof5z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof5z_sycl = (*bcof5z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof5z_sycl = (*ccof5z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcof5z_sycl = (*dcof5z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ovdelz_sycl = (*ovdelz_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class dfbydz_kernel_main_sycl>(
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

                const  ACC<double> functn(xdim0_dfbydz_kernel_main, ydim0_dfbydz_kernel_main, functn_p + (n_x * 1) + (n_y * xdim0_dfbydz_kernel_main * 1) + (n_z * xdim0_dfbydz_kernel_main * ydim0_dfbydz_kernel_main * 1));
                 ACC<double> fderiv(xdim1_dfbydz_kernel_main, ydim1_dfbydz_kernel_main, fderiv_p + (n_x * 1) + (n_y * xdim1_dfbydz_kernel_main * 1) + (n_z * xdim1_dfbydz_kernel_main * ydim1_dfbydz_kernel_main * 1));

                const int *nzglbl = &nzglbl_val;
                const int *nendzl = &nendzl_val;
                const int *nendzr = &nendzr_val;
                const int *nbound = &nbound_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double fdiffa;
    double fdiffb;
    double fdiffc;
    double fdiffd;
    double fdiffe;
    int kc;
    int kstart;
    int kfinis;

    kstart = 1;
    kfinis = nzglbl[0];
    if (nendzl[0] == nbound[0]) {
        kstart = 6;
    }
    if (nendzr[0] == nbound[0]) {
        kfinis = nzglbl[0] - 5;
    }
    kc = idx[2];
    if (kc >= kstart && kc <= kfinis) {
        fdiffa = functn(0, 0, 1) - functn(0, 0, -1);
        fdiffb = functn(0, 0, 2) - functn(0, 0, -2);
        fdiffc = functn(0, 0, 3) - functn(0, 0, -3);
        fdiffd = functn(0, 0, 4) - functn(0, 0, -4);
        fdiffe = functn(0, 0, 5) - functn(0, 0, -5);
        fderiv(0, 0, 0) = acoffz_sycl[0] * fdiffa + bcoffz_sycl[0] * fdiffb + ccoffz_sycl[0] * fdiffc + dcoffz_sycl[0] * fdiffd + ecoffz_sycl[0] * fdiffe;
    } else if (kc == 1) {
        fdiffa = functn(0, 0, 1) - functn(0, 0, 0);
        fdiffb = functn(0, 0, 2) - functn(0, 0, 0);
        fdiffc = functn(0, 0, 3) - functn(0, 0, 0);
        fdiffd = functn(0, 0, 4) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acof1z_sycl[0] * fdiffa + bcof1z_sycl[0] * fdiffb + ccof1z_sycl[0] * fdiffc + dcof1z_sycl[0] * fdiffd;
    } else if (kc == 2) {
        fdiffa = functn(0, 0, -1) - functn(0, 0, 0);
        fdiffb = functn(0, 0, 1) - functn(0, 0, 0);
        fdiffc = functn(0, 0, 2) - functn(0, 0, 0);
        fdiffd = functn(0, 0, 3) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acof2z_sycl[0] * fdiffa + bcof2z_sycl[0] * fdiffb + ccof2z_sycl[0] * fdiffc + dcof2z_sycl[0] * fdiffd;
    } else if (kc == 3) {
        fdiffa = functn(0, 0, 1) - functn(0, 0, -1);
        fdiffb = functn(0, 0, 2) - functn(0, 0, -2);
        fderiv(0, 0, 0) = acof3z_sycl[0] * fdiffa + bcof3z_sycl[0] * fdiffb;
    } else if (kc == 4) {
        fdiffa = functn(0, 0, 1) - functn(0, 0, -1);
        fdiffb = functn(0, 0, 2) - functn(0, 0, -2);
        fdiffc = functn(0, 0, 3) - functn(0, 0, -3);
        fderiv(0, 0, 0) = acof4z_sycl[0] * fdiffa + bcof4z_sycl[0] * fdiffb + ccof4z_sycl[0] * fdiffc;
    } else if (kc == 5) {
        fdiffa = functn(0, 0, 1) - functn(0, 0, -1);
        fdiffb = functn(0, 0, 2) - functn(0, 0, -2);
        fdiffc = functn(0, 0, 3) - functn(0, 0, -3);
        fdiffd = functn(0, 0, 4) - functn(0, 0, -4);
        fderiv(0, 0, 0) = acof5z_sycl[0] * fdiffa + bcof5z_sycl[0] * fdiffb + ccof5z_sycl[0] * fdiffc + dcof5z_sycl[0] * fdiffd;
    } else if (kc == nzglbl[0] - 4) {
        fdiffa = functn(0, 0, 1) - functn(0, 0, -1);
        fdiffb = functn(0, 0, 2) - functn(0, 0, -2);
        fdiffc = functn(0, 0, 3) - functn(0, 0, -3);
        fdiffd = functn(0, 0, 4) - functn(0, 0, -4);
        fderiv(0, 0, 0) = acof5z_sycl[0] * fdiffa + bcof5z_sycl[0] * fdiffb + ccof5z_sycl[0] * fdiffc + dcof5z_sycl[0] * fdiffd;
    } else if (kc == nzglbl[0] - 3) {
        fdiffa = functn(0, 0, 1) - functn(0, 0, -1);
        fdiffb = functn(0, 0, 2) - functn(0, 0, -2);
        fdiffc = functn(0, 0, 3) - functn(0, 0, -3);
        fderiv(0, 0, 0) = acof4z_sycl[0] * fdiffa + bcof4z_sycl[0] * fdiffb + ccof4z_sycl[0] * fdiffc;
    } else if (kc == nzglbl[0] - 2) {
        fdiffa = functn(0, 0, 1) - functn(0, 0, -1);
        fdiffb = functn(0, 0, 2) - functn(0, 0, -2);
        fderiv(0, 0, 0) = acof3z_sycl[0] * fdiffa + bcof3z_sycl[0] * fdiffb;
    } else if (kc == nzglbl[0] - 1) {
        fdiffa = functn(0, 0, 0) - functn(0, 0, 1);
        fdiffb = functn(0, 0, 0) - functn(0, 0, -1);
        fdiffc = functn(0, 0, 0) - functn(0, 0, -2);
        fdiffd = functn(0, 0, 0) - functn(0, 0, -3);
        fderiv(0, 0, 0) = acof2z_sycl[0] * fdiffa + bcof2z_sycl[0] * fdiffb + ccof2z_sycl[0] * fdiffc + dcof2z_sycl[0] * fdiffd;
    } else if (kc == nzglbl[0]) {
        fdiffa = functn(0, 0, 0) - functn(0, 0, -1);
        fdiffb = functn(0, 0, 0) - functn(0, 0, -2);
        fdiffc = functn(0, 0, 0) - functn(0, 0, -3);
        fdiffd = functn(0, 0, 0) - functn(0, 0, -4);
        fderiv(0, 0, 0) = acof1z_sycl[0] * fdiffa + bcof1z_sycl[0] * fdiffb + ccof1z_sycl[0] * fdiffc + dcof1z_sycl[0] * fdiffd;
    }
    fderiv(0, 0, 0) = fderiv(0, 0, 0) * ovdelz_sycl[0];

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[11].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 7);
    ops_set_halo_dirtybit3(&args[1], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[11].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[11].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[11].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void dfbydz_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("dfbydz_kernel_main", args, 7, 11, dim, 1, range, block, dfbydz_kernel_main_execute);
}
#endif

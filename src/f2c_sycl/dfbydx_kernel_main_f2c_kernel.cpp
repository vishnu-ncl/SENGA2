// Auto-generated at 2026-04-28 18:44:33.508183 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void dfbydx_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void dfbydx_kernel_main_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 7, range, 6)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 6, "dfbydx_kernel_main");
        block->instance->OPS_kernels[6].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "dfbydx_kernel_main");
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
    int xdim0_dfbydx_kernel_main = args[0].dat->size[0];
    int ydim0_dfbydx_kernel_main = args[0].dat->size[1];
    int xdim1_dfbydx_kernel_main = args[1].dat->size[0];
    int ydim1_dfbydx_kernel_main = args[1].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ functn_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ fderiv_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int nxglbl_val = *(int *)args[2].data;

    int nendxl_val = *(int *)args[3].data;

    int nendxr_val = *(int *)args[4].data;

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
        block->instance->OPS_kernels[6].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto acoffx_sycl = (*acoffx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcoffx_sycl = (*bcoffx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccoffx_sycl = (*ccoffx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcoffx_sycl = (*dcoffx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ecoffx_sycl = (*ecoffx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof1x_sycl = (*acof1x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof1x_sycl = (*bcof1x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof1x_sycl = (*ccof1x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcof1x_sycl = (*dcof1x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof2x_sycl = (*acof2x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof2x_sycl = (*bcof2x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof2x_sycl = (*ccof2x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcof2x_sycl = (*dcof2x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof3x_sycl = (*acof3x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof3x_sycl = (*bcof3x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof4x_sycl = (*acof4x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof4x_sycl = (*bcof4x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof4x_sycl = (*ccof4x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof5x_sycl = (*acof5x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof5x_sycl = (*bcof5x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof5x_sycl = (*ccof5x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcof5x_sycl = (*dcof5x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ovdelx_sycl = (*ovdelx_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class dfbydx_kernel_main_sycl>(
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

                const  ACC<double> functn(xdim0_dfbydx_kernel_main, ydim0_dfbydx_kernel_main, functn_p + (n_x * 1) + (n_y * xdim0_dfbydx_kernel_main * 1) + (n_z * xdim0_dfbydx_kernel_main * ydim0_dfbydx_kernel_main * 1));
                 ACC<double> fderiv(xdim1_dfbydx_kernel_main, ydim1_dfbydx_kernel_main, fderiv_p + (n_x * 1) + (n_y * xdim1_dfbydx_kernel_main * 1) + (n_z * xdim1_dfbydx_kernel_main * ydim1_dfbydx_kernel_main * 1));

                const int *nxglbl = &nxglbl_val;
                const int *nendxl = &nendxl_val;
                const int *nendxr = &nendxr_val;
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
    int ic;
    int istart;
    int ifinis;

    istart = 1;
    ifinis = nxglbl[0];
    if (nendxl[0] == nbound[0]) {
        istart = 6;
    }
    if (nendxr[0] == nbound[0]) {
        ifinis = nxglbl[0] - 5;
    }
    ic = idx[0];
    if (ic >= istart && ic <= ifinis) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(-3, 0, 0);
        fdiffd = functn(4, 0, 0) - functn(-4, 0, 0);
        fdiffe = functn(5, 0, 0) - functn(-5, 0, 0);
        fderiv(0, 0, 0) = acoffx_sycl[0] * fdiffa + bcoffx_sycl[0] * fdiffb + ccoffx_sycl[0] * fdiffc + dcoffx_sycl[0] * fdiffd + ecoffx_sycl[0] * fdiffe;
    } else if (ic == 1) {
        fdiffa = functn(1, 0, 0) - functn(0, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(0, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(0, 0, 0);
        fdiffd = functn(4, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acof1x_sycl[0] * fdiffa + bcof1x_sycl[0] * fdiffb + ccof1x_sycl[0] * fdiffc + dcof1x_sycl[0] * fdiffd;
    } else if (ic == 2) {
        fdiffa = functn(-1, 0, 0) - functn(0, 0, 0);
        fdiffb = functn(1, 0, 0) - functn(0, 0, 0);
        fdiffc = functn(2, 0, 0) - functn(0, 0, 0);
        fdiffd = functn(3, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acof2x_sycl[0] * fdiffa + bcof2x_sycl[0] * fdiffb + ccof2x_sycl[0] * fdiffc + dcof2x_sycl[0] * fdiffd;
    } else if (ic == 3) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fderiv(0, 0, 0) = acof3x_sycl[0] * fdiffa + bcof3x_sycl[0] * fdiffb;
    } else if (ic == 4) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(-3, 0, 0);
        fderiv(0, 0, 0) = acof4x_sycl[0] * fdiffa + bcof4x_sycl[0] * fdiffb + ccof4x_sycl[0] * fdiffc;
    } else if (ic == 5) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(-3, 0, 0);
        fdiffd = functn(4, 0, 0) - functn(-4, 0, 0);
        fderiv(0, 0, 0) = acof5x_sycl[0] * fdiffa + bcof5x_sycl[0] * fdiffb + ccof5x_sycl[0] * fdiffc + dcof5x_sycl[0] * fdiffd;
    } else if (ic == nxglbl[0] - 4) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(-3, 0, 0);
        fdiffd = functn(4, 0, 0) - functn(-4, 0, 0);
        fderiv(0, 0, 0) = acof5x_sycl[0] * fdiffa + bcof5x_sycl[0] * fdiffb + ccof5x_sycl[0] * fdiffc + dcof5x_sycl[0] * fdiffd;
    } else if (ic == nxglbl[0] - 3) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(3, 0, 0) - functn(-3, 0, 0);
        fderiv(0, 0, 0) = acof4x_sycl[0] * fdiffa + bcof4x_sycl[0] * fdiffb + ccof4x_sycl[0] * fdiffc;
    } else if (ic == nxglbl[0] - 2) {
        fdiffa = functn(1, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(2, 0, 0) - functn(-2, 0, 0);
        fderiv(0, 0, 0) = acof3x_sycl[0] * fdiffa + bcof3x_sycl[0] * fdiffb;
    } else if (ic == nxglbl[0] - 1) {
        fdiffa = functn(0, 0, 0) - functn(1, 0, 0);
        fdiffb = functn(0, 0, 0) - functn(-1, 0, 0);
        fdiffc = functn(0, 0, 0) - functn(-2, 0, 0);
        fdiffd = functn(0, 0, 0) - functn(-3, 0, 0);
        fderiv(0, 0, 0) = acof2x_sycl[0] * fdiffa + bcof2x_sycl[0] * fdiffb + ccof2x_sycl[0] * fdiffc + dcof2x_sycl[0] * fdiffd;
    } else if (ic == nxglbl[0]) {
        fdiffa = functn(0, 0, 0) - functn(-1, 0, 0);
        fdiffb = functn(0, 0, 0) - functn(-2, 0, 0);
        fdiffc = functn(0, 0, 0) - functn(-3, 0, 0);
        fdiffd = functn(0, 0, 0) - functn(-4, 0, 0);
        fderiv(0, 0, 0) = acof1x_sycl[0] * fdiffa + bcof1x_sycl[0] * fdiffb + ccof1x_sycl[0] * fdiffc + dcof1x_sycl[0] * fdiffd;
    }
    fderiv(0, 0, 0) = fderiv(0, 0, 0) * ovdelx_sycl[0];

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[6].time += __t1 - __t2;
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
        block->instance->OPS_kernels[6].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[6].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[6].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void dfbydx_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("dfbydx_kernel_main", args, 7, 6, dim, 1, range, block, dfbydx_kernel_main_execute);
}
#endif

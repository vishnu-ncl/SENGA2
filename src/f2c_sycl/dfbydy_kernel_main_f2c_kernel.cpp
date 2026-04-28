// Auto-generated at 2026-04-28 18:44:34.217285 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void dfbydy_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void dfbydy_kernel_main_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 7, range, 9)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 9, "dfbydy_kernel_main");
        block->instance->OPS_kernels[9].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "dfbydy_kernel_main");
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
    int xdim0_dfbydy_kernel_main = args[0].dat->size[0];
    int ydim0_dfbydy_kernel_main = args[0].dat->size[1];
    int xdim1_dfbydy_kernel_main = args[1].dat->size[0];
    int ydim1_dfbydy_kernel_main = args[1].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ functn_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ fderiv_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int nyglbl_val = *(int *)args[2].data;

    int nendyl_val = *(int *)args[3].data;

    int nendyr_val = *(int *)args[4].data;

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
        block->instance->OPS_kernels[9].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto acoffy_sycl = (*acoffy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcoffy_sycl = (*bcoffy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccoffy_sycl = (*ccoffy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcoffy_sycl = (*dcoffy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ecoffy_sycl = (*ecoffy_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof1y_sycl = (*acof1y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof1y_sycl = (*bcof1y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof1y_sycl = (*ccof1y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcof1y_sycl = (*dcof1y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof2y_sycl = (*acof2y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof2y_sycl = (*bcof2y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof2y_sycl = (*ccof2y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcof2y_sycl = (*dcof2y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof3y_sycl = (*acof3y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof3y_sycl = (*bcof3y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof4y_sycl = (*acof4y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof4y_sycl = (*bcof4y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof4y_sycl = (*ccof4y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acof5y_sycl = (*acof5y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcof5y_sycl = (*bcof5y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccof5y_sycl = (*ccof5y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcof5y_sycl = (*dcof5y_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ovdely_sycl = (*ovdely_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class dfbydy_kernel_main_sycl>(
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

                const  ACC<double> functn(xdim0_dfbydy_kernel_main, ydim0_dfbydy_kernel_main, functn_p + (n_x * 1) + (n_y * xdim0_dfbydy_kernel_main * 1) + (n_z * xdim0_dfbydy_kernel_main * ydim0_dfbydy_kernel_main * 1));
                 ACC<double> fderiv(xdim1_dfbydy_kernel_main, ydim1_dfbydy_kernel_main, fderiv_p + (n_x * 1) + (n_y * xdim1_dfbydy_kernel_main * 1) + (n_z * xdim1_dfbydy_kernel_main * ydim1_dfbydy_kernel_main * 1));

                const int *nyglbl = &nyglbl_val;
                const int *nendyl = &nendyl_val;
                const int *nendyr = &nendyr_val;
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
    int jc;
    int jstart;
    int jfinis;

    jstart = 1;
    jfinis = nyglbl[0];
    if (nendyl[0] == nbound[0]) {
        jstart = 6;
    }
    if (nendyr[0] == nbound[0]) {
        jfinis = nyglbl[0] - 5;
    }
    jc = idx[1];
    if (jc >= jstart && jc <= jfinis) {
        fdiffa = functn(0, 1, 0) - functn(0, -1, 0);
        fdiffb = functn(0, 2, 0) - functn(0, -2, 0);
        fdiffc = functn(0, 3, 0) - functn(0, -3, 0);
        fdiffd = functn(0, 4, 0) - functn(0, -4, 0);
        fdiffe = functn(0, 5, 0) - functn(0, -5, 0);
        fderiv(0, 0, 0) = acoffy_sycl[0] * fdiffa + bcoffy_sycl[0] * fdiffb + ccoffy_sycl[0] * fdiffc + dcoffy_sycl[0] * fdiffd + ecoffy_sycl[0] * fdiffe;
    } else if (jc == 1) {
        fdiffa = functn(0, 1, 0) - functn(0, 0, 0);
        fdiffb = functn(0, 2, 0) - functn(0, 0, 0);
        fdiffc = functn(0, 3, 0) - functn(0, 0, 0);
        fdiffd = functn(0, 4, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acof1y_sycl[0] * fdiffa + bcof1y_sycl[0] * fdiffb + ccof1y_sycl[0] * fdiffc + dcof1y_sycl[0] * fdiffd;
    } else if (jc == 2) {
        fdiffa = functn(0, -1, 0) - functn(0, 0, 0);
        fdiffb = functn(0, 1, 0) - functn(0, 0, 0);
        fdiffc = functn(0, 2, 0) - functn(0, 0, 0);
        fdiffd = functn(0, 3, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acof2y_sycl[0] * fdiffa + bcof2y_sycl[0] * fdiffb + ccof2y_sycl[0] * fdiffc + dcof2y_sycl[0] * fdiffd;
    } else if (jc == 3) {
        fdiffa = functn(0, 1, 0) - functn(0, -1, 0);
        fdiffb = functn(0, 2, 0) - functn(0, -2, 0);
        fderiv(0, 0, 0) = acof3y_sycl[0] * fdiffa + bcof3y_sycl[0] * fdiffb;
    } else if (jc == 4) {
        fdiffa = functn(0, 1, 0) - functn(0, -1, 0);
        fdiffb = functn(0, 2, 0) - functn(0, -2, 0);
        fdiffc = functn(0, 3, 0) - functn(0, -3, 0);
        fderiv(0, 0, 0) = acof4y_sycl[0] * fdiffa + bcof4y_sycl[0] * fdiffb + ccof4y_sycl[0] * fdiffc;
    } else if (jc == 5) {
        fdiffa = functn(0, 1, 0) - functn(0, -1, 0);
        fdiffb = functn(0, 2, 0) - functn(0, -2, 0);
        fdiffc = functn(0, 3, 0) - functn(0, -3, 0);
        fdiffd = functn(0, 4, 0) - functn(0, -4, 0);
        fderiv(0, 0, 0) = acof5y_sycl[0] * fdiffa + bcof5y_sycl[0] * fdiffb + ccof5y_sycl[0] * fdiffc + dcof5y_sycl[0] * fdiffd;
    } else if (jc == nyglbl[0] - 4) {
        fdiffa = functn(0, 1, 0) - functn(0, -1, 0);
        fdiffb = functn(0, 2, 0) - functn(0, -2, 0);
        fdiffc = functn(0, 3, 0) - functn(0, -3, 0);
        fdiffd = functn(0, 4, 0) - functn(0, -4, 0);
        fderiv(0, 0, 0) = acof5y_sycl[0] * fdiffa + bcof5y_sycl[0] * fdiffb + ccof5y_sycl[0] * fdiffc + dcof5y_sycl[0] * fdiffd;
    } else if (jc == nyglbl[0] - 3) {
        fdiffa = functn(0, 1, 0) - functn(0, -1, 0);
        fdiffb = functn(0, 2, 0) - functn(0, -2, 0);
        fdiffc = functn(0, 3, 0) - functn(0, -3, 0);
        fderiv(0, 0, 0) = acof4y_sycl[0] * fdiffa + bcof4y_sycl[0] * fdiffb + ccof4y_sycl[0] * fdiffc;
    } else if (jc == nyglbl[0] - 2) {
        fdiffa = functn(0, 1, 0) - functn(0, -1, 0);
        fdiffb = functn(0, 2, 0) - functn(0, -2, 0);
        fderiv(0, 0, 0) = acof3y_sycl[0] * fdiffa + bcof3y_sycl[0] * fdiffb;
    } else if (jc == nyglbl[0] - 1) {
        fdiffa = functn(0, 0, 0) - functn(0, 1, 0);
        fdiffb = functn(0, 0, 0) - functn(0, -1, 0);
        fdiffc = functn(0, 0, 0) - functn(0, -2, 0);
        fdiffd = functn(0, 0, 0) - functn(0, -3, 0);
        fderiv(0, 0, 0) = acof2y_sycl[0] * fdiffa + bcof2y_sycl[0] * fdiffb + ccof2y_sycl[0] * fdiffc + dcof2y_sycl[0] * fdiffd;
    } else if (jc == nyglbl[0]) {
        fdiffa = functn(0, 0, 0) - functn(0, -1, 0);
        fdiffb = functn(0, 0, 0) - functn(0, -2, 0);
        fdiffc = functn(0, 0, 0) - functn(0, -3, 0);
        fdiffd = functn(0, 0, 0) - functn(0, -4, 0);
        fderiv(0, 0, 0) = acof1y_sycl[0] * fdiffa + bcof1y_sycl[0] * fdiffb + ccof1y_sycl[0] * fdiffc + dcof1y_sycl[0] * fdiffd;
    }
    fderiv(0, 0, 0) = fderiv(0, 0, 0) * ovdely_sycl[0];

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[9].time += __t1 - __t2;
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
        block->instance->OPS_kernels[9].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[9].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[9].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void dfbydy_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("dfbydy_kernel_main", args, 7, 9, dim, 1, range, block, dfbydy_kernel_main_execute);
}
#endif

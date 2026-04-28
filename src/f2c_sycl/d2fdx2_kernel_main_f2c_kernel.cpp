// Auto-generated at 2026-04-28 18:44:34.020178 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void d2fdx2_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void d2fdx2_kernel_main_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 7, range, 7)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 7, "d2fdx2_kernel_main");
        block->instance->OPS_kernels[7].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "d2fdx2_kernel_main");
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
    int xdim0_d2fdx2_kernel_main = args[0].dat->size[0];
    int ydim0_d2fdx2_kernel_main = args[0].dat->size[1];
    int xdim1_d2fdx2_kernel_main = args[1].dat->size[0];
    int ydim1_d2fdx2_kernel_main = args[1].dat->size[1];

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
        block->instance->OPS_kernels[7].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto acofsx_sycl = (*acofsx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcofsx_sycl = (*bcofsx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccofsx_sycl = (*ccofsx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcofsx_sycl = (*dcofsx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ecofsx_sycl = (*ecofsx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acfs1x_sycl = (*acfs1x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcfs1x_sycl = (*bcfs1x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccfs1x_sycl = (*ccfs1x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcfs1x_sycl = (*dcfs1x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ecfs1x_sycl = (*ecfs1x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acfs2x_sycl = (*acfs2x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcfs2x_sycl = (*bcfs2x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccfs2x_sycl = (*ccfs2x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcfs2x_sycl = (*dcfs2x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ecfs2x_sycl = (*ecfs2x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acfs3x_sycl = (*acfs3x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcfs3x_sycl = (*bcfs3x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acfs4x_sycl = (*acfs4x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcfs4x_sycl = (*bcfs4x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccfs4x_sycl = (*ccfs4x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acfs5x_sycl = (*acfs5x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcfs5x_sycl = (*bcfs5x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccfs5x_sycl = (*ccfs5x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcfs5x_sycl = (*dcfs5x_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ovdlx2_sycl = (*ovdlx2_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class d2fdx2_kernel_main_sycl>(
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

                const  ACC<double> functn(xdim0_d2fdx2_kernel_main, ydim0_d2fdx2_kernel_main, functn_p + (n_x * 1) + (n_y * xdim0_d2fdx2_kernel_main * 1) + (n_z * xdim0_d2fdx2_kernel_main * ydim0_d2fdx2_kernel_main * 1));
                 ACC<double> fderiv(xdim1_d2fdx2_kernel_main, ydim1_d2fdx2_kernel_main, fderiv_p + (n_x * 1) + (n_y * xdim1_d2fdx2_kernel_main * 1) + (n_z * xdim1_d2fdx2_kernel_main * ydim1_d2fdx2_kernel_main * 1));

                const int *nxglbl = &nxglbl_val;
                const int *nendxl = &nendxl_val;
                const int *nendxr = &nendxr_val;
                const int *nbound = &nbound_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double fdifap;
    double fdifbp;
    double fdifcp;
    double fdifdp;
    double fdifep;
    double fdifam;
    double fdifbm;
    double fdifcm;
    double fdifdm;
    double fdifem;
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
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(-3, 0, 0);
        fdifdp = functn(4, 0, 0) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(-4, 0, 0);
        fdifep = functn(5, 0, 0) - functn(0, 0, 0);
        fdifem = functn(0, 0, 0) - functn(-5, 0, 0);
        fderiv(0, 0, 0) = acofsx_sycl[0] * (fdifap - fdifam) + bcofsx_sycl[0] * (fdifbp - fdifbm) + ccofsx_sycl[0] * (fdifcp - fdifcm) + dcofsx_sycl[0] * (fdifdp - fdifdm) + ecofsx_sycl[0] * (fdifep - fdifem);
    } else if (ic == 1) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifdp = functn(4, 0, 0) - functn(0, 0, 0);
        fdifep = functn(5, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs1x_sycl[0] * fdifap + bcfs1x_sycl[0] * fdifbp + ccfs1x_sycl[0] * fdifcp + dcfs1x_sycl[0] * fdifdp + ecfs1x_sycl[0] * fdifep;
    } else if (ic == 2) {
        fdifap = functn(-1, 0, 0) - functn(0, 0, 0);
        fdifbp = functn(1, 0, 0) - functn(0, 0, 0);
        fdifcp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifdp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifep = functn(4, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs2x_sycl[0] * fdifap + bcfs2x_sycl[0] * fdifbp + ccfs2x_sycl[0] * fdifcp + dcfs2x_sycl[0] * fdifdp + ecfs2x_sycl[0] * fdifep;
    } else if (ic == 3) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fderiv(0, 0, 0) = acfs3x_sycl[0] * (fdifap - fdifam) + bcfs3x_sycl[0] * (fdifbp - fdifbm);
    } else if (ic == 4) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(-3, 0, 0);
        fderiv(0, 0, 0) = acfs4x_sycl[0] * (fdifap - fdifam) + bcfs4x_sycl[0] * (fdifbp - fdifbm) + ccfs4x_sycl[0] * (fdifcp - fdifcm);
    } else if (ic == 5) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(-3, 0, 0);
        fdifdp = functn(4, 0, 0) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(-4, 0, 0);
        fderiv(0, 0, 0) = acfs5x_sycl[0] * (fdifap - fdifam) + bcfs5x_sycl[0] * (fdifbp - fdifbm) + ccfs5x_sycl[0] * (fdifcp - fdifcm) + dcfs5x_sycl[0] * (fdifdp - fdifdm);
    } else if (ic == nxglbl[0] - 4) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(-3, 0, 0);
        fdifdp = functn(4, 0, 0) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(-4, 0, 0);
        fderiv(0, 0, 0) = acfs5x_sycl[0] * (fdifap - fdifam) + bcfs5x_sycl[0] * (fdifbp - fdifbm) + ccfs5x_sycl[0] * (fdifcp - fdifcm) + dcfs5x_sycl[0] * (fdifdp - fdifdm);
    } else if (ic == nxglbl[0] - 3) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(-3, 0, 0);
        fderiv(0, 0, 0) = acfs4x_sycl[0] * (fdifap - fdifam) + bcfs4x_sycl[0] * (fdifbp - fdifbm) + ccfs4x_sycl[0] * (fdifcp - fdifcm);
    } else if (ic == nxglbl[0] - 2) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fderiv(0, 0, 0) = acfs3x_sycl[0] * (fdifap - fdifam) + bcfs3x_sycl[0] * (fdifbp - fdifbm);
    } else if (ic == nxglbl[0] - 1) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifbp = functn(-1, 0, 0) - functn(0, 0, 0);
        fdifcp = functn(-2, 0, 0) - functn(0, 0, 0);
        fdifdp = functn(-3, 0, 0) - functn(0, 0, 0);
        fdifep = functn(-4, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs2x_sycl[0] * fdifap + bcfs2x_sycl[0] * fdifbp + ccfs2x_sycl[0] * fdifcp + dcfs2x_sycl[0] * fdifdp + ecfs2x_sycl[0] * fdifep;
    } else if (ic == nxglbl[0]) {
        fdifap = functn(-1, 0, 0) - functn(0, 0, 0);
        fdifbp = functn(-2, 0, 0) - functn(0, 0, 0);
        fdifcp = functn(-3, 0, 0) - functn(0, 0, 0);
        fdifdp = functn(-4, 0, 0) - functn(0, 0, 0);
        fdifep = functn(-5, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs1x_sycl[0] * fdifap + bcfs1x_sycl[0] * fdifbp + ccfs1x_sycl[0] * fdifcp + dcfs1x_sycl[0] * fdifdp + ecfs1x_sycl[0] * fdifep;
    }
    fderiv(0, 0, 0) = fderiv(0, 0, 0) * ovdlx2_sycl[0];

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[7].time += __t1 - __t2;
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
        block->instance->OPS_kernels[7].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[7].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[7].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void d2fdx2_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("d2fdx2_kernel_main", args, 7, 7, dim, 1, range, block, d2fdx2_kernel_main_execute);
}
#endif

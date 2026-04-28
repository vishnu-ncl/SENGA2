// Auto-generated at 2026-04-28 18:44:34.954761 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void d2fdz2_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void d2fdz2_kernel_main_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 7, range, 15)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 15, "d2fdz2_kernel_main");
        block->instance->OPS_kernels[15].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "d2fdz2_kernel_main");
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
    int xdim0_d2fdz2_kernel_main = args[0].dat->size[0];
    int ydim0_d2fdz2_kernel_main = args[0].dat->size[1];
    int xdim1_d2fdz2_kernel_main = args[1].dat->size[0];
    int ydim1_d2fdz2_kernel_main = args[1].dat->size[1];

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
        block->instance->OPS_kernels[15].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto acofsz_sycl = (*acofsz_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcofsz_sycl = (*bcofsz_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccofsz_sycl = (*ccofsz_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcofsz_sycl = (*dcofsz_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ecofsz_sycl = (*ecofsz_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acfs1z_sycl = (*acfs1z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcfs1z_sycl = (*bcfs1z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccfs1z_sycl = (*ccfs1z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcfs1z_sycl = (*dcfs1z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ecfs1z_sycl = (*ecfs1z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acfs2z_sycl = (*acfs2z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcfs2z_sycl = (*bcfs2z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccfs2z_sycl = (*ccfs2z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcfs2z_sycl = (*dcfs2z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ecfs2z_sycl = (*ecfs2z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acfs3z_sycl = (*acfs3z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcfs3z_sycl = (*bcfs3z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acfs4z_sycl = (*acfs4z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcfs4z_sycl = (*bcfs4z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccfs4z_sycl = (*ccfs4z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto acfs5z_sycl = (*acfs5z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto bcfs5z_sycl = (*bcfs5z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ccfs5z_sycl = (*ccfs5z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto dcfs5z_sycl = (*dcfs5z_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ovdlz2_sycl = (*ovdlz2_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class d2fdz2_kernel_main_sycl>(
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

                const  ACC<double> functn(xdim0_d2fdz2_kernel_main, ydim0_d2fdz2_kernel_main, functn_p + (n_x * 1) + (n_y * xdim0_d2fdz2_kernel_main * 1) + (n_z * xdim0_d2fdz2_kernel_main * ydim0_d2fdz2_kernel_main * 1));
                 ACC<double> fderiv(xdim1_d2fdz2_kernel_main, ydim1_d2fdz2_kernel_main, fderiv_p + (n_x * 1) + (n_y * xdim1_d2fdz2_kernel_main * 1) + (n_z * xdim1_d2fdz2_kernel_main * ydim1_d2fdz2_kernel_main * 1));

                const int *nzglbl = &nzglbl_val;
                const int *nendzl = &nendzl_val;
                const int *nendzr = &nendzr_val;
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
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(0, 0, -3);
        fdifdp = functn(0, 0, 4) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(0, 0, -4);
        fdifep = functn(0, 0, 5) - functn(0, 0, 0);
        fdifem = functn(0, 0, 0) - functn(0, 0, -5);
        fderiv(0, 0, 0) = acofsz_sycl[0] * (fdifap - fdifam) + bcofsz_sycl[0] * (fdifbp - fdifbm) + ccofsz_sycl[0] * (fdifcp - fdifcm) + dcofsz_sycl[0] * (fdifdp - fdifdm) + ecofsz_sycl[0] * (fdifep - fdifem);
    } else if (kc == 1) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifdp = functn(0, 0, 4) - functn(0, 0, 0);
        fdifep = functn(0, 0, 5) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs1z_sycl[0] * fdifap + bcfs1z_sycl[0] * fdifbp + ccfs1z_sycl[0] * fdifcp + dcfs1z_sycl[0] * fdifdp + ecfs1z_sycl[0] * fdifep;
    } else if (kc == 2) {
        fdifap = functn(0, 0, -1) - functn(0, 0, 0);
        fdifbp = functn(0, 0, 1) - functn(0, 0, 0);
        fdifcp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifdp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifep = functn(0, 0, 4) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs2z_sycl[0] * fdifap + bcfs2z_sycl[0] * fdifbp + ccfs2z_sycl[0] * fdifcp + dcfs2z_sycl[0] * fdifdp + ecfs2z_sycl[0] * fdifep;
    } else if (kc == 3) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fderiv(0, 0, 0) = acfs3z_sycl[0] * (fdifap - fdifam) + bcfs3z_sycl[0] * (fdifbp - fdifbm);
    } else if (kc == 4) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(0, 0, -3);
        fderiv(0, 0, 0) = acfs4z_sycl[0] * (fdifap - fdifam) + bcfs4z_sycl[0] * (fdifbp - fdifbm) + ccfs4z_sycl[0] * (fdifcp - fdifcm);
    } else if (kc == 5) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(0, 0, -3);
        fdifdp = functn(0, 0, 4) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(0, 0, -4);
        fderiv(0, 0, 0) = acfs5z_sycl[0] * (fdifap - fdifam) + bcfs5z_sycl[0] * (fdifbp - fdifbm) + ccfs5z_sycl[0] * (fdifcp - fdifcm) + dcfs5z_sycl[0] * (fdifdp - fdifdm);
    } else if (kc == nzglbl[0] - 4) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(0, 0, -3);
        fdifdp = functn(0, 0, 4) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(0, 0, -4);
        fderiv(0, 0, 0) = acfs5z_sycl[0] * (fdifap - fdifam) + bcfs5z_sycl[0] * (fdifbp - fdifbm) + ccfs5z_sycl[0] * (fdifcp - fdifcm) + dcfs5z_sycl[0] * (fdifdp - fdifdm);
    } else if (kc == nzglbl[0] - 3) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fdifcp = functn(0, 0, 3) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(0, 0, -3);
        fderiv(0, 0, 0) = acfs4z_sycl[0] * (fdifap - fdifam) + bcfs4z_sycl[0] * (fdifbp - fdifbm) + ccfs4z_sycl[0] * (fdifcp - fdifcm);
    } else if (kc == nzglbl[0] - 2) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(0, 0, -1);
        fdifbp = functn(0, 0, 2) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(0, 0, -2);
        fderiv(0, 0, 0) = acfs3z_sycl[0] * (fdifap - fdifam) + bcfs3z_sycl[0] * (fdifbp - fdifbm);
    } else if (kc == nzglbl[0] - 1) {
        fdifap = functn(0, 0, 1) - functn(0, 0, 0);
        fdifbp = functn(0, 0, -1) - functn(0, 0, 0);
        fdifcp = functn(0, 0, -2) - functn(0, 0, 0);
        fdifdp = functn(0, 0, -3) - functn(0, 0, 0);
        fdifep = functn(0, 0, -4) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs2z_sycl[0] * fdifap + bcfs2z_sycl[0] * fdifbp + ccfs2z_sycl[0] * fdifcp + dcfs2z_sycl[0] * fdifdp + ecfs2z_sycl[0] * fdifep;
    } else if (kc == nzglbl[0]) {
        fdifap = functn(0, 0, -1) - functn(0, 0, 0);
        fdifbp = functn(0, 0, -2) - functn(0, 0, 0);
        fdifcp = functn(0, 0, -3) - functn(0, 0, 0);
        fdifdp = functn(0, 0, -4) - functn(0, 0, 0);
        fdifep = functn(0, 0, -5) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs1z_sycl[0] * fdifap + bcfs1z_sycl[0] * fdifbp + ccfs1z_sycl[0] * fdifcp + dcfs1z_sycl[0] * fdifdp + ecfs1z_sycl[0] * fdifep;
    }
    fderiv(0, 0, 0) = fderiv(0, 0, 0) * ovdlz2_sycl[0];

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[15].time += __t1 - __t2;
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
        block->instance->OPS_kernels[15].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[15].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[15].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void d2fdz2_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("d2fdz2_kernel_main", args, 7, 15, dim, 1, range, block, d2fdz2_kernel_main_execute);
}
#endif

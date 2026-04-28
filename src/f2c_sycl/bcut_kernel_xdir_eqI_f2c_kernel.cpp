// Auto-generated at 2026-04-28 18:44:49.235328 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bcut_kernel_xdir_eqI_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bcut_kernel_xdir_eqI_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[18];
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
    args[12] = desc->args[12];
    args[13] = desc->args[13];
    args[14] = desc->args[14];
    args[15] = desc->args[15];
    args[16] = desc->args[16];
    args[17] = desc->args[17];
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
    ops_arg arg12 = args[12];
    ops_arg arg13 = args[13];
    ops_arg arg14 = args[14];
    ops_arg arg15 = args[15];
    ops_arg arg16 = args[16];
    ops_arg arg17 = args[17];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 18, range, 326)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 326, "bcut_kernel_xdir_eqI");
        block->instance->OPS_kernels[326].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "bcut_kernel_xdir_eqI");
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
    int xdim0_bcut_kernel_xdir_eqI = args[0].dat->size[0];
    int ydim0_bcut_kernel_xdir_eqI = args[0].dat->size[1];
    int xdim1_bcut_kernel_xdir_eqI = args[1].dat->size[0];
    int ydim1_bcut_kernel_xdir_eqI = args[1].dat->size[1];
    int xdim2_bcut_kernel_xdir_eqI = args[2].dat->size[0];
    int ydim2_bcut_kernel_xdir_eqI = args[2].dat->size[1];
    int xdim3_bcut_kernel_xdir_eqI = args[3].dat->size[0];
    int ydim3_bcut_kernel_xdir_eqI = args[3].dat->size[1];
    int xdim4_bcut_kernel_xdir_eqI = args[4].dat->size[0];
    int ydim4_bcut_kernel_xdir_eqI = args[4].dat->size[1];
    int xdim5_bcut_kernel_xdir_eqI = args[5].dat->size[0];
    int ydim5_bcut_kernel_xdir_eqI = args[5].dat->size[1];
    int xdim6_bcut_kernel_xdir_eqI = args[6].dat->size[0];
    int ydim6_bcut_kernel_xdir_eqI = args[6].dat->size[1];
    int xdim7_bcut_kernel_xdir_eqI = args[7].dat->size[0];
    int ydim7_bcut_kernel_xdir_eqI = args[7].dat->size[1];
    int xdim8_bcut_kernel_xdir_eqI = args[8].dat->size[0];
    int ydim8_bcut_kernel_xdir_eqI = args[8].dat->size[1];
    int xdim9_bcut_kernel_xdir_eqI = args[9].dat->size[0];
    int ydim9_bcut_kernel_xdir_eqI = args[9].dat->size[1];
    int xdim10_bcut_kernel_xdir_eqI = args[10].dat->size[0];
    int ydim10_bcut_kernel_xdir_eqI = args[10].dat->size[1];
    int xdim11_bcut_kernel_xdir_eqI = args[11].dat->size[0];
    int ydim11_bcut_kernel_xdir_eqI = args[11].dat->size[1];
    int xdim12_bcut_kernel_xdir_eqI = args[12].dat->size[0];
    int ydim12_bcut_kernel_xdir_eqI = args[12].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ strux_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ strvx_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ strwx_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ dudtx_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ dvdtx_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ dwdtx_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ ustead_p = (double *)(args[6].data_d) + base6 - 1; // Subtracting 1 to convert to C-style

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ uinf1_p = (double *)(args[7].data_d) + base7 - 1; // Subtracting 1 to convert to C-style

    int base8 = getDatBaseFromOpsArg3D(&args[8], start_indx, 1);
    double * __restrict__ vinf1_p = (double *)(args[8].data_d) + base8 - 1; // Subtracting 1 to convert to C-style

    int base9 = getDatBaseFromOpsArg3D(&args[9], start_indx, 1);
    double * __restrict__ winf1_p = (double *)(args[9].data_d) + base9 - 1; // Subtracting 1 to convert to C-style

    int base10 = getDatBaseFromOpsArg3D(&args[10], start_indx, 1);
    double * __restrict__ uinf2_p = (double *)(args[10].data_d) + base10 - 1; // Subtracting 1 to convert to C-style

    int base11 = getDatBaseFromOpsArg3D(&args[11], start_indx, 1);
    double * __restrict__ vinf2_p = (double *)(args[11].data_d) + base11 - 1; // Subtracting 1 to convert to C-style

    int base12 = getDatBaseFromOpsArg3D(&args[12], start_indx, 1);
    double * __restrict__ winf2_p = (double *)(args[12].data_d) + base12 - 1; // Subtracting 1 to convert to C-style

    double rxlprm_1_val = *(double *)args[13].data;

    double vfac_val = *(double *)args[14].data;

    double coflow_val = *(double *)args[15].data;

    double tstep_val = *(double *)args[16].data;

    double dvfdt_val = *(double *)args[17].data;

//  Subtracting 1 here as start_indx and end_indx is in Fortran style - converting it to c-style range
    int start_0 = start_indx[0]-1;
    int end_0 = end_indx[0];
    int start_1 = start_indx[1]-1;
    int end_1 = end_indx[1];
    int start_2 = start_indx[2]-1;
    int end_2 = end_indx[2];

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 18);
    ops_halo_exchanges(args, 18, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[326].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            cgh.parallel_for<class bcut_kernel_xdir_eqI_sycl>(
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

                 ACC<double> strux(xdim0_bcut_kernel_xdir_eqI, ydim0_bcut_kernel_xdir_eqI, strux_p + (n_x * 0) + (n_y * xdim0_bcut_kernel_xdir_eqI * 1) + (n_z * xdim0_bcut_kernel_xdir_eqI * ydim0_bcut_kernel_xdir_eqI * 1));
                 ACC<double> strvx(xdim1_bcut_kernel_xdir_eqI, ydim1_bcut_kernel_xdir_eqI, strvx_p + (n_x * 0) + (n_y * xdim1_bcut_kernel_xdir_eqI * 1) + (n_z * xdim1_bcut_kernel_xdir_eqI * ydim1_bcut_kernel_xdir_eqI * 1));
                 ACC<double> strwx(xdim2_bcut_kernel_xdir_eqI, ydim2_bcut_kernel_xdir_eqI, strwx_p + (n_x * 0) + (n_y * xdim2_bcut_kernel_xdir_eqI * 1) + (n_z * xdim2_bcut_kernel_xdir_eqI * ydim2_bcut_kernel_xdir_eqI * 1));
                 ACC<double> dudtx(xdim3_bcut_kernel_xdir_eqI, ydim3_bcut_kernel_xdir_eqI, dudtx_p + (n_x * 0) + (n_y * xdim3_bcut_kernel_xdir_eqI * 1) + (n_z * xdim3_bcut_kernel_xdir_eqI * ydim3_bcut_kernel_xdir_eqI * 1));
                 ACC<double> dvdtx(xdim4_bcut_kernel_xdir_eqI, ydim4_bcut_kernel_xdir_eqI, dvdtx_p + (n_x * 0) + (n_y * xdim4_bcut_kernel_xdir_eqI * 1) + (n_z * xdim4_bcut_kernel_xdir_eqI * ydim4_bcut_kernel_xdir_eqI * 1));
                 ACC<double> dwdtx(xdim5_bcut_kernel_xdir_eqI, ydim5_bcut_kernel_xdir_eqI, dwdtx_p + (n_x * 0) + (n_y * xdim5_bcut_kernel_xdir_eqI * 1) + (n_z * xdim5_bcut_kernel_xdir_eqI * ydim5_bcut_kernel_xdir_eqI * 1));
                const  ACC<double> ustead(xdim6_bcut_kernel_xdir_eqI, ydim6_bcut_kernel_xdir_eqI, ustead_p + (n_x * 0) + (n_y * xdim6_bcut_kernel_xdir_eqI * 1) + (n_z * xdim6_bcut_kernel_xdir_eqI * ydim6_bcut_kernel_xdir_eqI * 1));
                const  ACC<double> uinf1(xdim7_bcut_kernel_xdir_eqI, ydim7_bcut_kernel_xdir_eqI, uinf1_p + (n_x * 0) + (n_y * xdim7_bcut_kernel_xdir_eqI * 1) + (n_z * xdim7_bcut_kernel_xdir_eqI * ydim7_bcut_kernel_xdir_eqI * 1));
                const  ACC<double> vinf1(xdim8_bcut_kernel_xdir_eqI, ydim8_bcut_kernel_xdir_eqI, vinf1_p + (n_x * 0) + (n_y * xdim8_bcut_kernel_xdir_eqI * 1) + (n_z * xdim8_bcut_kernel_xdir_eqI * ydim8_bcut_kernel_xdir_eqI * 1));
                const  ACC<double> winf1(xdim9_bcut_kernel_xdir_eqI, ydim9_bcut_kernel_xdir_eqI, winf1_p + (n_x * 0) + (n_y * xdim9_bcut_kernel_xdir_eqI * 1) + (n_z * xdim9_bcut_kernel_xdir_eqI * ydim9_bcut_kernel_xdir_eqI * 1));
                const  ACC<double> uinf2(xdim10_bcut_kernel_xdir_eqI, ydim10_bcut_kernel_xdir_eqI, uinf2_p + (n_x * 0) + (n_y * xdim10_bcut_kernel_xdir_eqI * 1) + (n_z * xdim10_bcut_kernel_xdir_eqI * ydim10_bcut_kernel_xdir_eqI * 1));
                const  ACC<double> vinf2(xdim11_bcut_kernel_xdir_eqI, ydim11_bcut_kernel_xdir_eqI, vinf2_p + (n_x * 0) + (n_y * xdim11_bcut_kernel_xdir_eqI * 1) + (n_z * xdim11_bcut_kernel_xdir_eqI * ydim11_bcut_kernel_xdir_eqI * 1));
                const  ACC<double> winf2(xdim12_bcut_kernel_xdir_eqI, ydim12_bcut_kernel_xdir_eqI, winf2_p + (n_x * 0) + (n_y * xdim12_bcut_kernel_xdir_eqI * 1) + (n_z * xdim12_bcut_kernel_xdir_eqI * ydim12_bcut_kernel_xdir_eqI * 1));

                const double *rxlprm_1 = &rxlprm_1_val;
                const double *vfac = &vfac_val;
                const double *coflow = &coflow_val;
                const double *tstep = &tstep_val;
                const double *dvfdt = &dvfdt_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    strux(0, 0, 0) = (rxlprm_1[0] * ustead(0, 0, 0) + uinf2(0, 0, 0)) * vfac[0] + coflow[0];
    strvx(0, 0, 0) = vinf2(0, 0, 0) * vfac[0];
    strwx(0, 0, 0) = winf2(0, 0, 0) * vfac[0];
    dudtx(0, 0, 0) = (uinf2(0, 0, 0) - uinf1(0, 0, 0)) / tstep[0] * vfac[0] + (rxlprm_1[0] * ustead(0, 0, 0) + uinf2(0, 0, 0)) * dvfdt[0];
    dvdtx(0, 0, 0) = (vinf2(0, 0, 0) - vinf1(0, 0, 0)) / tstep[0] * vfac[0] + vinf2(0, 0, 0) * dvfdt[0];
    dwdtx(0, 0, 0) = (winf2(0, 0, 0) - winf1(0, 0, 0)) / tstep[0] * vfac[0] + winf2(0, 0, 0) * dvfdt[0];

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[326].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 18);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
    ops_set_halo_dirtybit3(&args[3], range);
    ops_set_halo_dirtybit3(&args[4], range);
    ops_set_halo_dirtybit3(&args[5], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[326].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg8);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg9);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg10);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg11);
        block->instance->OPS_kernels[326].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg12);
    }
}

#ifdef OPS_LAZY
extern "C"
void bcut_kernel_xdir_eqI_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bcut_kernel_xdir_eqI", args, 18, 326, dim, 1, range, block, bcut_kernel_xdir_eqI_execute);
}
#endif

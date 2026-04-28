// Auto-generated at 2026-04-28 18:44:52.631158 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bounds_kernel_eqP_yl_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bounds_kernel_eqP_yl_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 18, range, 398)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 398, "bounds_kernel_eqP_yl");
        block->instance->OPS_kernels[398].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "bounds_kernel_eqP_yl");
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
    int xdim0_bounds_kernel_eqP_yl = args[0].dat->size[0];
    int ydim0_bounds_kernel_eqP_yl = args[0].dat->size[1];
    int xdim1_bounds_kernel_eqP_yl = args[1].dat->size[0];
    int ydim1_bounds_kernel_eqP_yl = args[1].dat->size[1];
    int xdim2_bounds_kernel_eqP_yl = args[2].dat->size[0];
    int ydim2_bounds_kernel_eqP_yl = args[2].dat->size[1];
    int xdim3_bounds_kernel_eqP_yl = args[3].dat->size[0];
    int ydim3_bounds_kernel_eqP_yl = args[3].dat->size[1];
    int xdim4_bounds_kernel_eqP_yl = args[4].dat->size[0];
    int ydim4_bounds_kernel_eqP_yl = args[4].dat->size[1];
    int xdim5_bounds_kernel_eqP_yl = args[5].dat->size[0];
    int ydim5_bounds_kernel_eqP_yl = args[5].dat->size[1];
    int xdim6_bounds_kernel_eqP_yl = args[6].dat->size[0];
    int ydim6_bounds_kernel_eqP_yl = args[6].dat->size[1];
    int xdim7_bounds_kernel_eqP_yl = args[7].dat->size[0];
    int ydim7_bounds_kernel_eqP_yl = args[7].dat->size[1];
    int xdim8_bounds_kernel_eqP_yl = args[8].dat->size[0];
    int ydim8_bounds_kernel_eqP_yl = args[8].dat->size[1];
    int xdim9_bounds_kernel_eqP_yl = args[9].dat->size[0];
    int ydim9_bounds_kernel_eqP_yl = args[9].dat->size[1];
    int xdim10_bounds_kernel_eqP_yl = args[10].dat->size[0];
    int ydim10_bounds_kernel_eqP_yl = args[10].dat->size[1];
    int xdim11_bounds_kernel_eqP_yl = args[11].dat->size[0];
    int ydim11_bounds_kernel_eqP_yl = args[11].dat->size[1];
    int xdim12_bounds_kernel_eqP_yl = args[12].dat->size[0];
    int ydim12_bounds_kernel_eqP_yl = args[12].dat->size[1];
    int xdim13_bounds_kernel_eqP_yl = args[13].dat->size[0];
    int ydim13_bounds_kernel_eqP_yl = args[13].dat->size[1];
    int xdim14_bounds_kernel_eqP_yl = args[14].dat->size[0];
    int ydim14_bounds_kernel_eqP_yl = args[14].dat->size[1];
    int xdim15_bounds_kernel_eqP_yl = args[15].dat->size[0];
    int ydim15_bounds_kernel_eqP_yl = args[15].dat->size[1];
    int xdim16_bounds_kernel_eqP_yl = args[16].dat->size[0];
    int ydim16_bounds_kernel_eqP_yl = args[16].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ urhs_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ vrhs_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ wrhs_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ bcl2yl_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ bcl3yl_p = (double *)(args[6].data_d) + base6 - 1; // Subtracting 1 to convert to C-style

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ bcl4yl_p = (double *)(args[7].data_d) + base7 - 1; // Subtracting 1 to convert to C-style

    int base8 = getDatBaseFromOpsArg3D(&args[8], start_indx, 1);
    double * __restrict__ bcl5yl_p = (double *)(args[8].data_d) + base8 - 1; // Subtracting 1 to convert to C-style

    int base9 = getDatBaseFromOpsArg3D(&args[9], start_indx, 1);
    double * __restrict__ strdyl_p = (double *)(args[9].data_d) + base9 - 1; // Subtracting 1 to convert to C-style

    int base10 = getDatBaseFromOpsArg3D(&args[10], start_indx, 1);
    double * __restrict__ struyl_p = (double *)(args[10].data_d) + base10 - 1; // Subtracting 1 to convert to C-style

    int base11 = getDatBaseFromOpsArg3D(&args[11], start_indx, 1);
    double * __restrict__ strvyl_p = (double *)(args[11].data_d) + base11 - 1; // Subtracting 1 to convert to C-style

    int base12 = getDatBaseFromOpsArg3D(&args[12], start_indx, 1);
    double * __restrict__ strwyl_p = (double *)(args[12].data_d) + base12 - 1; // Subtracting 1 to convert to C-style

    int base13 = getDatBaseFromOpsArg3D(&args[13], start_indx, 1);
    double * __restrict__ streyl_p = (double *)(args[13].data_d) + base13 - 1; // Subtracting 1 to convert to C-style

    int base14 = getDatBaseFromOpsArg3D(&args[14], start_indx, 1);
    double * __restrict__ ova2yl_p = (double *)(args[14].data_d) + base14 - 1; // Subtracting 1 to convert to C-style

    int base15 = getDatBaseFromOpsArg3D(&args[15], start_indx, 1);
    double * __restrict__ ovgmyl_p = (double *)(args[15].data_d) + base15 - 1; // Subtracting 1 to convert to C-style

    int base16 = getDatBaseFromOpsArg3D(&args[16], start_indx, 1);
    double * __restrict__ acouyl_p = (double *)(args[16].data_d) + base16 - 1; // Subtracting 1 to convert to C-style

    int flag_pio_yl_val = *(int *)args[17].data;

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
        block->instance->OPS_kernels[398].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            cgh.parallel_for<class bounds_kernel_eqP_yl_sycl>(
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

                 ACC<double> drhs(xdim0_bounds_kernel_eqP_yl, ydim0_bounds_kernel_eqP_yl, drhs_p + (n_x * 1) + (n_y * xdim0_bounds_kernel_eqP_yl * 1) + (n_z * xdim0_bounds_kernel_eqP_yl * ydim0_bounds_kernel_eqP_yl * 1));
                 ACC<double> urhs(xdim1_bounds_kernel_eqP_yl, ydim1_bounds_kernel_eqP_yl, urhs_p + (n_x * 1) + (n_y * xdim1_bounds_kernel_eqP_yl * 1) + (n_z * xdim1_bounds_kernel_eqP_yl * ydim1_bounds_kernel_eqP_yl * 1));
                 ACC<double> vrhs(xdim2_bounds_kernel_eqP_yl, ydim2_bounds_kernel_eqP_yl, vrhs_p + (n_x * 1) + (n_y * xdim2_bounds_kernel_eqP_yl * 1) + (n_z * xdim2_bounds_kernel_eqP_yl * ydim2_bounds_kernel_eqP_yl * 1));
                 ACC<double> wrhs(xdim3_bounds_kernel_eqP_yl, ydim3_bounds_kernel_eqP_yl, wrhs_p + (n_x * 1) + (n_y * xdim3_bounds_kernel_eqP_yl * 1) + (n_z * xdim3_bounds_kernel_eqP_yl * ydim3_bounds_kernel_eqP_yl * 1));
                 ACC<double> erhs(xdim4_bounds_kernel_eqP_yl, ydim4_bounds_kernel_eqP_yl, erhs_p + (n_x * 1) + (n_y * xdim4_bounds_kernel_eqP_yl * 1) + (n_z * xdim4_bounds_kernel_eqP_yl * ydim4_bounds_kernel_eqP_yl * 1));
                const  ACC<double> bcl2yl(xdim5_bounds_kernel_eqP_yl, ydim5_bounds_kernel_eqP_yl, bcl2yl_p + (n_x * 1) + (n_y * xdim5_bounds_kernel_eqP_yl * 0) + (n_z * xdim5_bounds_kernel_eqP_yl * ydim5_bounds_kernel_eqP_yl * 1));
                const  ACC<double> bcl3yl(xdim6_bounds_kernel_eqP_yl, ydim6_bounds_kernel_eqP_yl, bcl3yl_p + (n_x * 1) + (n_y * xdim6_bounds_kernel_eqP_yl * 0) + (n_z * xdim6_bounds_kernel_eqP_yl * ydim6_bounds_kernel_eqP_yl * 1));
                const  ACC<double> bcl4yl(xdim7_bounds_kernel_eqP_yl, ydim7_bounds_kernel_eqP_yl, bcl4yl_p + (n_x * 1) + (n_y * xdim7_bounds_kernel_eqP_yl * 0) + (n_z * xdim7_bounds_kernel_eqP_yl * ydim7_bounds_kernel_eqP_yl * 1));
                const  ACC<double> bcl5yl(xdim8_bounds_kernel_eqP_yl, ydim8_bounds_kernel_eqP_yl, bcl5yl_p + (n_x * 1) + (n_y * xdim8_bounds_kernel_eqP_yl * 0) + (n_z * xdim8_bounds_kernel_eqP_yl * ydim8_bounds_kernel_eqP_yl * 1));
                const  ACC<double> strdyl(xdim9_bounds_kernel_eqP_yl, ydim9_bounds_kernel_eqP_yl, strdyl_p + (n_x * 1) + (n_y * xdim9_bounds_kernel_eqP_yl * 0) + (n_z * xdim9_bounds_kernel_eqP_yl * ydim9_bounds_kernel_eqP_yl * 1));
                const  ACC<double> struyl(xdim10_bounds_kernel_eqP_yl, ydim10_bounds_kernel_eqP_yl, struyl_p + (n_x * 1) + (n_y * xdim10_bounds_kernel_eqP_yl * 0) + (n_z * xdim10_bounds_kernel_eqP_yl * ydim10_bounds_kernel_eqP_yl * 1));
                const  ACC<double> strvyl(xdim11_bounds_kernel_eqP_yl, ydim11_bounds_kernel_eqP_yl, strvyl_p + (n_x * 1) + (n_y * xdim11_bounds_kernel_eqP_yl * 0) + (n_z * xdim11_bounds_kernel_eqP_yl * ydim11_bounds_kernel_eqP_yl * 1));
                const  ACC<double> strwyl(xdim12_bounds_kernel_eqP_yl, ydim12_bounds_kernel_eqP_yl, strwyl_p + (n_x * 1) + (n_y * xdim12_bounds_kernel_eqP_yl * 0) + (n_z * xdim12_bounds_kernel_eqP_yl * ydim12_bounds_kernel_eqP_yl * 1));
                const  ACC<double> streyl(xdim13_bounds_kernel_eqP_yl, ydim13_bounds_kernel_eqP_yl, streyl_p + (n_x * 1) + (n_y * xdim13_bounds_kernel_eqP_yl * 0) + (n_z * xdim13_bounds_kernel_eqP_yl * ydim13_bounds_kernel_eqP_yl * 1));
                const  ACC<double> ova2yl(xdim14_bounds_kernel_eqP_yl, ydim14_bounds_kernel_eqP_yl, ova2yl_p + (n_x * 1) + (n_y * xdim14_bounds_kernel_eqP_yl * 0) + (n_z * xdim14_bounds_kernel_eqP_yl * ydim14_bounds_kernel_eqP_yl * 1));
                const  ACC<double> ovgmyl(xdim15_bounds_kernel_eqP_yl, ydim15_bounds_kernel_eqP_yl, ovgmyl_p + (n_x * 1) + (n_y * xdim15_bounds_kernel_eqP_yl * 0) + (n_z * xdim15_bounds_kernel_eqP_yl * ydim15_bounds_kernel_eqP_yl * 1));
                const  ACC<double> acouyl(xdim16_bounds_kernel_eqP_yl, ydim16_bounds_kernel_eqP_yl, acouyl_p + (n_x * 1) + (n_y * xdim16_bounds_kernel_eqP_yl * 0) + (n_z * xdim16_bounds_kernel_eqP_yl * ydim16_bounds_kernel_eqP_yl * 1));

                const int *flag_pio_yl = &flag_pio_yl_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    if ((strvyl(0, 0, 0) > 0.0) && (flag_pio_yl[0] == 1)) {
        drhs(0, 0, 0) = drhs(0, 0, 0) - bcl2yl(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0);
        urhs(0, 0, 0) = urhs(0, 0, 0) - bcl2yl(0, 0, 0) * struyl(0, 0, 0) - bcl3yl(0, 0, 0) * strdyl(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0) * struyl(0, 0, 0);
        vrhs(0, 0, 0) = vrhs(0, 0, 0) - bcl2yl(0, 0, 0) * strvyl(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0) * (strvyl(0, 0, 0) + acouyl(0, 0, 0));
        wrhs(0, 0, 0) = wrhs(0, 0, 0) - bcl2yl(0, 0, 0) * strwyl(0, 0, 0) - bcl4yl(0, 0, 0) * strdyl(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0) * strwyl(0, 0, 0);
        erhs(0, 0, 0) = erhs(0, 0, 0) - bcl2yl(0, 0, 0) * streyl(0, 0, 0) - bcl3yl(0, 0, 0) * strdyl(0, 0, 0) * struyl(0, 0, 0) - bcl4yl(0, 0, 0) * strdyl(0, 0, 0) * strwyl(0, 0, 0) - bcl5yl(0, 0, 0) * (ova2yl(0, 0, 0) * streyl(0, 0, 0) + strvyl(0, 0, 0) / acouyl(0, 0, 0) + ovgmyl(0, 0, 0));
    } else {
        drhs(0, 0, 0) = drhs(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0);
        urhs(0, 0, 0) = urhs(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0) * struyl(0, 0, 0);
        vrhs(0, 0, 0) = vrhs(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0) * (strvyl(0, 0, 0) + acouyl(0, 0, 0));
        wrhs(0, 0, 0) = wrhs(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0) * strwyl(0, 0, 0);
        erhs(0, 0, 0) = erhs(0, 0, 0) - bcl5yl(0, 0, 0) * (ova2yl(0, 0, 0) * streyl(0, 0, 0) + strvyl(0, 0, 0) / acouyl(0, 0, 0) + ovgmyl(0, 0, 0));
    }

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[398].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 18);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
    ops_set_halo_dirtybit3(&args[3], range);
    ops_set_halo_dirtybit3(&args[4], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[398].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg8);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg9);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg10);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg11);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg12);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg13);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg14);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg15);
        block->instance->OPS_kernels[398].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg16);
    }
}

#ifdef OPS_LAZY
extern "C"
void bounds_kernel_eqP_yl_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bounds_kernel_eqP_yl", args, 18, 398, dim, 1, range, block, bounds_kernel_eqP_yl_execute);
}
#endif

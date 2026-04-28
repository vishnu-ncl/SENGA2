// Auto-generated at 2026-04-28 18:44:50.492771 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bounds_kernel_eqAB_xdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bounds_kernel_eqAB_xdir_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[27];
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
    args[18] = desc->args[18];
    args[19] = desc->args[19];
    args[20] = desc->args[20];
    args[21] = desc->args[21];
    args[22] = desc->args[22];
    args[23] = desc->args[23];
    args[24] = desc->args[24];
    args[25] = desc->args[25];
    args[26] = desc->args[26];
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
    ops_arg arg18 = args[18];
    ops_arg arg19 = args[19];
    ops_arg arg20 = args[20];
    ops_arg arg21 = args[21];
    ops_arg arg22 = args[22];
    ops_arg arg23 = args[23];
    ops_arg arg24 = args[24];
    ops_arg arg25 = args[25];
    ops_arg arg26 = args[26];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 27, range, 356)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 356, "bounds_kernel_eqAB_xdir");
        block->instance->OPS_kernels[356].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "bounds_kernel_eqAB_xdir");
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
    int xdim0_bounds_kernel_eqAB_xdir = args[0].dat->size[0];
    int ydim0_bounds_kernel_eqAB_xdir = args[0].dat->size[1];
    int xdim1_bounds_kernel_eqAB_xdir = args[1].dat->size[0];
    int ydim1_bounds_kernel_eqAB_xdir = args[1].dat->size[1];
    int xdim2_bounds_kernel_eqAB_xdir = args[2].dat->size[0];
    int ydim2_bounds_kernel_eqAB_xdir = args[2].dat->size[1];
    int xdim3_bounds_kernel_eqAB_xdir = args[3].dat->size[0];
    int ydim3_bounds_kernel_eqAB_xdir = args[3].dat->size[1];
    int xdim4_bounds_kernel_eqAB_xdir = args[4].dat->size[0];
    int ydim4_bounds_kernel_eqAB_xdir = args[4].dat->size[1];
    int xdim5_bounds_kernel_eqAB_xdir = args[5].dat->size[0];
    int ydim5_bounds_kernel_eqAB_xdir = args[5].dat->size[1];
    int xdim6_bounds_kernel_eqAB_xdir = args[6].dat->size[0];
    int ydim6_bounds_kernel_eqAB_xdir = args[6].dat->size[1];
    int xdim7_bounds_kernel_eqAB_xdir = args[7].dat->size[0];
    int ydim7_bounds_kernel_eqAB_xdir = args[7].dat->size[1];
    int xdim8_bounds_kernel_eqAB_xdir = args[8].dat->size[0];
    int ydim8_bounds_kernel_eqAB_xdir = args[8].dat->size[1];
    int xdim9_bounds_kernel_eqAB_xdir = args[9].dat->size[0];
    int ydim9_bounds_kernel_eqAB_xdir = args[9].dat->size[1];
    int xdim10_bounds_kernel_eqAB_xdir = args[10].dat->size[0];
    int ydim10_bounds_kernel_eqAB_xdir = args[10].dat->size[1];
    int xdim11_bounds_kernel_eqAB_xdir = args[11].dat->size[0];
    int ydim11_bounds_kernel_eqAB_xdir = args[11].dat->size[1];
    int xdim12_bounds_kernel_eqAB_xdir = args[12].dat->size[0];
    int ydim12_bounds_kernel_eqAB_xdir = args[12].dat->size[1];
    int xdim13_bounds_kernel_eqAB_xdir = args[13].dat->size[0];
    int ydim13_bounds_kernel_eqAB_xdir = args[13].dat->size[1];
    int xdim14_bounds_kernel_eqAB_xdir = args[14].dat->size[0];
    int ydim14_bounds_kernel_eqAB_xdir = args[14].dat->size[1];
    int xdim15_bounds_kernel_eqAB_xdir = args[15].dat->size[0];
    int ydim15_bounds_kernel_eqAB_xdir = args[15].dat->size[1];
    int xdim16_bounds_kernel_eqAB_xdir = args[16].dat->size[0];
    int ydim16_bounds_kernel_eqAB_xdir = args[16].dat->size[1];
    int xdim17_bounds_kernel_eqAB_xdir = args[17].dat->size[0];
    int ydim17_bounds_kernel_eqAB_xdir = args[17].dat->size[1];
    int xdim18_bounds_kernel_eqAB_xdir = args[18].dat->size[0];
    int ydim18_bounds_kernel_eqAB_xdir = args[18].dat->size[1];
    int xdim19_bounds_kernel_eqAB_xdir = args[19].dat->size[0];
    int ydim19_bounds_kernel_eqAB_xdir = args[19].dat->size[1];
    int xdim20_bounds_kernel_eqAB_xdir = args[20].dat->size[0];
    int ydim20_bounds_kernel_eqAB_xdir = args[20].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ bcl2x_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ bcl3x_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ bcl4x_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ bcl5x_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ bcl1x_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ strdx_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ acoux_p = (double *)(args[6].data_d) + base6 - 1; // Subtracting 1 to convert to C-style

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ ova2x_p = (double *)(args[7].data_d) + base7 - 1; // Subtracting 1 to convert to C-style

    int base8 = getDatBaseFromOpsArg3D(&args[8], start_indx, 1);
    double * __restrict__ tt2x_p = (double *)(args[8].data_d) + base8 - 1; // Subtracting 1 to convert to C-style

    int base9 = getDatBaseFromOpsArg3D(&args[9], start_indx, 1);
    double * __restrict__ tt3x_p = (double *)(args[9].data_d) + base9 - 1; // Subtracting 1 to convert to C-style

    int base10 = getDatBaseFromOpsArg3D(&args[10], start_indx, 1);
    double * __restrict__ tt4x_p = (double *)(args[10].data_d) + base10 - 1; // Subtracting 1 to convert to C-style

    int base11 = getDatBaseFromOpsArg3D(&args[11], start_indx, 1);
    double * __restrict__ tt5x_p = (double *)(args[11].data_d) + base11 - 1; // Subtracting 1 to convert to C-style

    int base12 = getDatBaseFromOpsArg3D(&args[12], start_indx, 1);
    double * __restrict__ strrx_p = (double *)(args[12].data_d) + base12 - 1; // Subtracting 1 to convert to C-style

    int base13 = getDatBaseFromOpsArg3D(&args[13], start_indx, 1);
    double * __restrict__ stlux_p = (double *)(args[13].data_d) + base13 - 1; // Subtracting 1 to convert to C-style

    int base14 = getDatBaseFromOpsArg3D(&args[14], start_indx, 1);
    double * __restrict__ strux_p = (double *)(args[14].data_d) + base14 - 1; // Subtracting 1 to convert to C-style

    int base15 = getDatBaseFromOpsArg3D(&args[15], start_indx, 1);
    double * __restrict__ stlvx_p = (double *)(args[15].data_d) + base15 - 1; // Subtracting 1 to convert to C-style

    int base16 = getDatBaseFromOpsArg3D(&args[16], start_indx, 1);
    double * __restrict__ strvx_p = (double *)(args[16].data_d) + base16 - 1; // Subtracting 1 to convert to C-style

    int base17 = getDatBaseFromOpsArg3D(&args[17], start_indx, 1);
    double * __restrict__ stlwx_p = (double *)(args[17].data_d) + base17 - 1; // Subtracting 1 to convert to C-style

    int base18 = getDatBaseFromOpsArg3D(&args[18], start_indx, 1);
    double * __restrict__ strwx_p = (double *)(args[18].data_d) + base18 - 1; // Subtracting 1 to convert to C-style

    int base19 = getDatBaseFromOpsArg3D(&args[19], start_indx, 1);
    double * __restrict__ stltx_p = (double *)(args[19].data_d) + base19 - 1; // Subtracting 1 to convert to C-style

    int base20 = getDatBaseFromOpsArg3D(&args[20], start_indx, 1);
    double * __restrict__ strtx_p = (double *)(args[20].data_d) + base20 - 1; // Subtracting 1 to convert to C-style

    double xgdlen_val = *(double *)args[21].data;

    double nrieta2_val = *(double *)args[22].data;

    double nrieta3_val = *(double *)args[23].data;

    double nrieta4_val = *(double *)args[24].data;

    double nrieta5_val = *(double *)args[25].data;

    double m2max_val = *(double *)args[26].data;

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
    ops_H_D_exchanges_device(args, 27);
    ops_halo_exchanges(args, 27, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[356].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            cgh.parallel_for<class bounds_kernel_eqAB_xdir_sycl>(
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

                 ACC<double> bcl2x(xdim0_bounds_kernel_eqAB_xdir, ydim0_bounds_kernel_eqAB_xdir, bcl2x_p + (n_x * 0) + (n_y * xdim0_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim0_bounds_kernel_eqAB_xdir * ydim0_bounds_kernel_eqAB_xdir * 1));
                 ACC<double> bcl3x(xdim1_bounds_kernel_eqAB_xdir, ydim1_bounds_kernel_eqAB_xdir, bcl3x_p + (n_x * 0) + (n_y * xdim1_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim1_bounds_kernel_eqAB_xdir * ydim1_bounds_kernel_eqAB_xdir * 1));
                 ACC<double> bcl4x(xdim2_bounds_kernel_eqAB_xdir, ydim2_bounds_kernel_eqAB_xdir, bcl4x_p + (n_x * 0) + (n_y * xdim2_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim2_bounds_kernel_eqAB_xdir * ydim2_bounds_kernel_eqAB_xdir * 1));
                 ACC<double> bcl5x(xdim3_bounds_kernel_eqAB_xdir, ydim3_bounds_kernel_eqAB_xdir, bcl5x_p + (n_x * 0) + (n_y * xdim3_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim3_bounds_kernel_eqAB_xdir * ydim3_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> bcl1x(xdim4_bounds_kernel_eqAB_xdir, ydim4_bounds_kernel_eqAB_xdir, bcl1x_p + (n_x * 0) + (n_y * xdim4_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim4_bounds_kernel_eqAB_xdir * ydim4_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> strdx(xdim5_bounds_kernel_eqAB_xdir, ydim5_bounds_kernel_eqAB_xdir, strdx_p + (n_x * 0) + (n_y * xdim5_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim5_bounds_kernel_eqAB_xdir * ydim5_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> acoux(xdim6_bounds_kernel_eqAB_xdir, ydim6_bounds_kernel_eqAB_xdir, acoux_p + (n_x * 0) + (n_y * xdim6_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim6_bounds_kernel_eqAB_xdir * ydim6_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> ova2x(xdim7_bounds_kernel_eqAB_xdir, ydim7_bounds_kernel_eqAB_xdir, ova2x_p + (n_x * 0) + (n_y * xdim7_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim7_bounds_kernel_eqAB_xdir * ydim7_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> tt2x(xdim8_bounds_kernel_eqAB_xdir, ydim8_bounds_kernel_eqAB_xdir, tt2x_p + (n_x * 0) + (n_y * xdim8_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim8_bounds_kernel_eqAB_xdir * ydim8_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> tt3x(xdim9_bounds_kernel_eqAB_xdir, ydim9_bounds_kernel_eqAB_xdir, tt3x_p + (n_x * 0) + (n_y * xdim9_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim9_bounds_kernel_eqAB_xdir * ydim9_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> tt4x(xdim10_bounds_kernel_eqAB_xdir, ydim10_bounds_kernel_eqAB_xdir, tt4x_p + (n_x * 0) + (n_y * xdim10_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim10_bounds_kernel_eqAB_xdir * ydim10_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> tt5x(xdim11_bounds_kernel_eqAB_xdir, ydim11_bounds_kernel_eqAB_xdir, tt5x_p + (n_x * 0) + (n_y * xdim11_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim11_bounds_kernel_eqAB_xdir * ydim11_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> strrx(xdim12_bounds_kernel_eqAB_xdir, ydim12_bounds_kernel_eqAB_xdir, strrx_p + (n_x * 0) + (n_y * xdim12_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim12_bounds_kernel_eqAB_xdir * ydim12_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> stlux(xdim13_bounds_kernel_eqAB_xdir, ydim13_bounds_kernel_eqAB_xdir, stlux_p + (n_x * 0) + (n_y * xdim13_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim13_bounds_kernel_eqAB_xdir * ydim13_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> strux(xdim14_bounds_kernel_eqAB_xdir, ydim14_bounds_kernel_eqAB_xdir, strux_p + (n_x * 0) + (n_y * xdim14_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim14_bounds_kernel_eqAB_xdir * ydim14_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> stlvx(xdim15_bounds_kernel_eqAB_xdir, ydim15_bounds_kernel_eqAB_xdir, stlvx_p + (n_x * 0) + (n_y * xdim15_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim15_bounds_kernel_eqAB_xdir * ydim15_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> strvx(xdim16_bounds_kernel_eqAB_xdir, ydim16_bounds_kernel_eqAB_xdir, strvx_p + (n_x * 0) + (n_y * xdim16_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim16_bounds_kernel_eqAB_xdir * ydim16_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> stlwx(xdim17_bounds_kernel_eqAB_xdir, ydim17_bounds_kernel_eqAB_xdir, stlwx_p + (n_x * 0) + (n_y * xdim17_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim17_bounds_kernel_eqAB_xdir * ydim17_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> strwx(xdim18_bounds_kernel_eqAB_xdir, ydim18_bounds_kernel_eqAB_xdir, strwx_p + (n_x * 0) + (n_y * xdim18_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim18_bounds_kernel_eqAB_xdir * ydim18_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> stltx(xdim19_bounds_kernel_eqAB_xdir, ydim19_bounds_kernel_eqAB_xdir, stltx_p + (n_x * 0) + (n_y * xdim19_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim19_bounds_kernel_eqAB_xdir * ydim19_bounds_kernel_eqAB_xdir * 1));
                const  ACC<double> strtx(xdim20_bounds_kernel_eqAB_xdir, ydim20_bounds_kernel_eqAB_xdir, strtx_p + (n_x * 0) + (n_y * xdim20_bounds_kernel_eqAB_xdir * 1) + (n_z * xdim20_bounds_kernel_eqAB_xdir * ydim20_bounds_kernel_eqAB_xdir * 1));

                const double *xgdlen = &xgdlen_val;
                const double *nrieta2 = &nrieta2_val;
                const double *nrieta3 = &nrieta3_val;
                const double *nrieta4 = &nrieta4_val;
                const double *nrieta5 = &nrieta5_val;
                const double *m2max = &m2max_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double fornow;

    fornow = strdx(0, 0, 0) * acoux(0, 0, 0) * bcl1x(0, 0, 0);
    bcl2x(0, 0, 0) = stlux(0, 0, 0) * (bcl2x(0, 0, 0) - bcl5x(0, 0, 0) * ova2x(0, 0, 0));
    bcl3x(0, 0, 0) = stlux(0, 0, 0) * bcl3x(0, 0, 0);
    bcl4x(0, 0, 0) = stlux(0, 0, 0) * bcl4x(0, 0, 0);
    bcl5x(0, 0, 0) = 0.5 * (stlux(0, 0, 0) + acoux(0, 0, 0)) * (bcl5x(0, 0, 0) + fornow);
    bcl2x(0, 0, 0) = nrieta2[0] * strdx(0, 0, 0) * strrx(0, 0, 0) / (xgdlen[0] * acoux(0, 0, 0)) * (stltx(0, 0, 0) - strtx(0, 0, 0)) + tt2x(0, 0, 0) * ova2x(0, 0, 0) - bcl2x(0, 0, 0);
    bcl3x(0, 0, 0) = nrieta3[0] * acoux(0, 0, 0) / xgdlen[0] * (stlvx(0, 0, 0) - strvx(0, 0, 0)) + tt3x(0, 0, 0) - bcl3x(0, 0, 0);
    bcl4x(0, 0, 0) = nrieta4[0] * acoux(0, 0, 0) / xgdlen[0] * (stlwx(0, 0, 0) - strwx(0, 0, 0)) + tt4x(0, 0, 0) - bcl4x(0, 0, 0);
    bcl5x(0, 0, 0) = nrieta5[0] * strdx(0, 0, 0) * acoux(0, 0, 0) * acoux(0, 0, 0) * (1.0 - m2max[0]) / (2.0 * xgdlen[0]) * (stlux(0, 0, 0) - strux(0, 0, 0)) + 0.5 * tt5x(0, 0, 0) - bcl5x(0, 0, 0);

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[356].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 27);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
    ops_set_halo_dirtybit3(&args[3], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[356].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg8);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg9);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg10);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg11);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg12);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg13);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg14);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg15);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg16);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg17);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg18);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg19);
        block->instance->OPS_kernels[356].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg20);
    }
}

#ifdef OPS_LAZY
extern "C"
void bounds_kernel_eqAB_xdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bounds_kernel_eqAB_xdir", args, 27, 356, dim, 1, range, block, bounds_kernel_eqAB_xdir_execute);
}
#endif

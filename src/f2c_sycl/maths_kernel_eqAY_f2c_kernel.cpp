// Auto-generated at 2026-04-28 18:44:58.743968 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void maths_kernel_eqAY_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void maths_kernel_eqAY_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[19];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 19, range, 540)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 540, "maths_kernel_eqAY");
        block->instance->OPS_kernels[540].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "maths_kernel_eqAY");
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
    int xdim0_maths_kernel_eqAY = args[0].dat->size[0];
    int ydim0_maths_kernel_eqAY = args[0].dat->size[1];
    int xdim1_maths_kernel_eqAY = args[1].dat->size[0];
    int ydim1_maths_kernel_eqAY = args[1].dat->size[1];
    int xdim2_maths_kernel_eqAY = args[2].dat->size[0];
    int ydim2_maths_kernel_eqAY = args[2].dat->size[1];
    int xdim3_maths_kernel_eqAY = args[3].dat->size[0];
    int ydim3_maths_kernel_eqAY = args[3].dat->size[1];
    int xdim4_maths_kernel_eqAY = args[4].dat->size[0];
    int ydim4_maths_kernel_eqAY = args[4].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ out_arr1_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ out_arr2_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ out_arr3_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ in_arr1_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ in_arr2_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    double racnst_val = *(double *)args[5].data;

    double rncnst_val = *(double *)args[6].data;

    double reovrr_val = *(double *)args[7].data;

    double talpha_val = *(double *)args[8].data;

    double ovtst1_val = *(double *)args[9].data;

    double tstar2_val = *(double *)args[10].data;

    double ovtst3_val = *(double *)args[11].data;

    double cfcst1_val = *(double *)args[12].data;

    double cfcst2_val = *(double *)args[13].data;

    double encst1_val = *(double *)args[14].data;

    double encst2_val = *(double *)args[15].data;

    double dtcnst_val = *(double *)args[16].data;

    double omalph_val = *(double *)args[17].data;

    double clnten_val = *(double *)args[18].data;

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
    ops_H_D_exchanges_device(args, 19);
    ops_halo_exchanges(args, 19, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[540].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            cgh.parallel_for<class maths_kernel_eqAY_sycl>(
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

                 ACC<double> out_arr1(xdim0_maths_kernel_eqAY, ydim0_maths_kernel_eqAY, out_arr1_p + (n_x * 1) + (n_y * xdim0_maths_kernel_eqAY * 1) + (n_z * xdim0_maths_kernel_eqAY * ydim0_maths_kernel_eqAY * 1));
                 ACC<double> out_arr2(xdim1_maths_kernel_eqAY, ydim1_maths_kernel_eqAY, out_arr2_p + (n_x * 1) + (n_y * xdim1_maths_kernel_eqAY * 1) + (n_z * xdim1_maths_kernel_eqAY * ydim1_maths_kernel_eqAY * 1));
                 ACC<double> out_arr3(xdim2_maths_kernel_eqAY, ydim2_maths_kernel_eqAY, out_arr3_p + (n_x * 1) + (n_y * xdim2_maths_kernel_eqAY * 1) + (n_z * xdim2_maths_kernel_eqAY * ydim2_maths_kernel_eqAY * 1));
                const  ACC<double> in_arr1(xdim3_maths_kernel_eqAY, ydim3_maths_kernel_eqAY, in_arr1_p + (n_x * 1) + (n_y * xdim3_maths_kernel_eqAY * 1) + (n_z * xdim3_maths_kernel_eqAY * ydim3_maths_kernel_eqAY * 1));
                const  ACC<double> in_arr2(xdim4_maths_kernel_eqAY, ydim4_maths_kernel_eqAY, in_arr2_p + (n_x * 1) + (n_y * xdim4_maths_kernel_eqAY * 1) + (n_z * xdim4_maths_kernel_eqAY * ydim4_maths_kernel_eqAY * 1));

                const double *racnst = &racnst_val;
                const double *rncnst = &rncnst_val;
                const double *reovrr = &reovrr_val;
                const double *talpha = &talpha_val;
                const double *ovtst1 = &ovtst1_val;
                const double *tstar2 = &tstar2_val;
                const double *ovtst3 = &ovtst3_val;
                const double *cfcst1 = &cfcst1_val;
                const double *cfcst2 = &cfcst2_val;
                const double *encst1 = &encst1_val;
                const double *encst2 = &encst2_val;
                const double *dtcnst = &dtcnst_val;
                const double *omalph = &omalph_val;
                const double *clnten = &clnten_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double preduc;
    double fornow;
    double trats1;
    double trats2;
    double trats3;
    double ftcent;
    double cfactr;
    double enfact;
    double fbroad;

    out_arr3(0, 0, 0) = racnst[0] + rncnst[0] * cl::sycl::log(in_arr1(0, 0, 0)) - reovrr[0] / in_arr1(0, 0, 0);
    out_arr3(0, 0, 0) = cl::sycl::exp(out_arr3(0, 0, 0));
    preduc = in_arr2(0, 0, 0) * out_arr3(0, 0, 0) / out_arr1(0, 0, 0);
    trats1 = in_arr1(0, 0, 0) * ovtst1[0];
    trats2 = -tstar2[0] / in_arr1(0, 0, 0);
    trats3 = in_arr1(0, 0, 0) * ovtst3[0];
    ftcent = omalph[0] * cl::sycl::exp(trats3) + talpha[0] * cl::sycl::exp(trats1) + cl::sycl::exp(trats2);
    ftcent = cl::sycl::log10(ftcent);
    cfactr = cfcst1[0] + cfcst2[0] * ftcent;
    enfact = encst1[0] + encst2[0] * ftcent;
    fornow = cl::sycl::log10(preduc) + cfactr;
    fornow = fornow / (enfact - dtcnst[0] * fornow);
    fornow = 1.0 + fornow * fornow;
    fbroad = ftcent / fornow;
    fbroad = cl::sycl::exp(fbroad * clnten[0]);
    fornow = fbroad * preduc / (1.0 + preduc);
    out_arr1(0, 0, 0) = out_arr1(0, 0, 0) * fornow;
    out_arr2(0, 0, 0) = cl::sycl::log(out_arr1(0, 0, 0));

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[540].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 19);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[540].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[540].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[540].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[540].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[540].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[540].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
    }
}

#ifdef OPS_LAZY
extern "C"
void maths_kernel_eqAY_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("maths_kernel_eqAY", args, 19, 540, dim, 1, range, block, maths_kernel_eqAY_execute);
}
#endif

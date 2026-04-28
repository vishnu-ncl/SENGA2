// Auto-generated at 2026-04-28 18:43:28.359932 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bcut_kernel_xdir_eqH_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bcut_kernel_xdir_eqH_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[16];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 16, range, 325)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 325, "bcut_kernel_xdir_eqH");
        block->instance->OPS_kernels[325].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "bcut_kernel_xdir_eqH");
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
    int xdim0_bcut_kernel_xdir_eqH = args[0].dat->size[0];
    int ydim0_bcut_kernel_xdir_eqH = args[0].dat->size[1];
    int xdim1_bcut_kernel_xdir_eqH = args[1].dat->size[0];
    int ydim1_bcut_kernel_xdir_eqH = args[1].dat->size[1];
    int xdim2_bcut_kernel_xdir_eqH = args[2].dat->size[0];
    int ydim2_bcut_kernel_xdir_eqH = args[2].dat->size[1];
    int xdim3_bcut_kernel_xdir_eqH = args[3].dat->size[0];
    int ydim3_bcut_kernel_xdir_eqH = args[3].dat->size[1];
    int xdim4_bcut_kernel_xdir_eqH = args[4].dat->size[0];
    int ydim4_bcut_kernel_xdir_eqH = args[4].dat->size[1];
    int xdim5_bcut_kernel_xdir_eqH = args[5].dat->size[0];
    int ydim5_bcut_kernel_xdir_eqH = args[5].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ strux_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ strvx_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ strwx_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ dudtx_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ dvdtx_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ dwdtx_p = (double *)(args[5].data) + base5 - 1; // Subtracting 1 to convert to C-style

    double * __restrict__  slope = (double *)args[6].data;
    double * __restrict__  widthp = (double *)args[7].data;
    double * __restrict__  ptly = (double *)args[8].data;
    double * __restrict__  rxlprm_1 = (double *)args[9].data;
    double * __restrict__  rxlprc_7 = (double *)args[10].data;
    double * __restrict__  fornow = (double *)args[11].data;
    double * __restrict__  argmnt = (double *)args[12].data;
    double * __restrict__  deltagy = (double *)args[13].data;
    double * __restrict__  ygdlen = (double *)args[14].data;

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 16);
    ops_halo_exchanges(args, 16, range);
    ops_H_D_exchanges_host(args, 16);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[325].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {
                int idx[] = {arg_idx[0] + n_x, arg_idx[1] + n_y, arg_idx[2] + n_z};

                 ACC<double> strux(xdim0_bcut_kernel_xdir_eqH, ydim0_bcut_kernel_xdir_eqH, strux_p + (n_x * 0) + (n_y * xdim0_bcut_kernel_xdir_eqH * 1) + (n_z * xdim0_bcut_kernel_xdir_eqH * ydim0_bcut_kernel_xdir_eqH * 1));

                 ACC<double> strvx(xdim1_bcut_kernel_xdir_eqH, ydim1_bcut_kernel_xdir_eqH, strvx_p + (n_x * 0) + (n_y * xdim1_bcut_kernel_xdir_eqH * 1) + (n_z * xdim1_bcut_kernel_xdir_eqH * ydim1_bcut_kernel_xdir_eqH * 1));

                 ACC<double> strwx(xdim2_bcut_kernel_xdir_eqH, ydim2_bcut_kernel_xdir_eqH, strwx_p + (n_x * 0) + (n_y * xdim2_bcut_kernel_xdir_eqH * 1) + (n_z * xdim2_bcut_kernel_xdir_eqH * ydim2_bcut_kernel_xdir_eqH * 1));

                 ACC<double> dudtx(xdim3_bcut_kernel_xdir_eqH, ydim3_bcut_kernel_xdir_eqH, dudtx_p + (n_x * 0) + (n_y * xdim3_bcut_kernel_xdir_eqH * 1) + (n_z * xdim3_bcut_kernel_xdir_eqH * ydim3_bcut_kernel_xdir_eqH * 1));

                 ACC<double> dvdtx(xdim4_bcut_kernel_xdir_eqH, ydim4_bcut_kernel_xdir_eqH, dvdtx_p + (n_x * 0) + (n_y * xdim4_bcut_kernel_xdir_eqH * 1) + (n_z * xdim4_bcut_kernel_xdir_eqH * ydim4_bcut_kernel_xdir_eqH * 1));

                 ACC<double> dwdtx(xdim5_bcut_kernel_xdir_eqH, ydim5_bcut_kernel_xdir_eqH, dwdtx_p + (n_x * 0) + (n_y * xdim5_bcut_kernel_xdir_eqH * 1) + (n_z * xdim5_bcut_kernel_xdir_eqH * ydim5_bcut_kernel_xdir_eqH * 1));

    int jc;
    double hyptan;

    jc = idx[1];
    hyptan = (0.5 * f2c::tanh(slope[0] * (f2c::dble(jc - 1) * deltagy[0] - (ygdlen[0] / 2.0 - widthp[0] * ptly[0]))) + 0.5) * (-0.5 * f2c::tanh(slope[0] * (f2c::dble(jc - 1) * deltagy[0] - (ygdlen[0] / 2.0 + widthp[0] * ptly[0]))) + 0.5);
    strux(0, 0, 0) = rxlprm_1[0] + hyptan * rxlprc_7[0] * f2c::sin(argmnt[0]);
    strvx(0, 0, 0) = 0.0;
    strwx(0, 0, 0) = 0.0;
    dudtx(0, 0, 0) = hyptan * fornow[0] * rxlprc_7[0] * f2c::cos(argmnt[0]);
    dvdtx(0, 0, 0) = 0.0;
    dwdtx(0, 0, 0) = 0.0;

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[325].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 16);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
    ops_set_halo_dirtybit3(&args[3], range);
    ops_set_halo_dirtybit3(&args[4], range);
    ops_set_halo_dirtybit3(&args[5], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[325].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
    }
}

#ifdef OPS_LAZY
extern "C"
void bcut_kernel_xdir_eqH_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bcut_kernel_xdir_eqH", args, 16, 325, dim, 0, range, block, bcut_kernel_xdir_eqH_execute);
}
#endif

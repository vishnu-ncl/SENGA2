// Auto-generated at 2026-04-28 18:43:30.695223 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bounds_kernel_eqU_xr_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bounds_kernel_eqU_xr_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[13];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 13, range, 380)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 380, "bounds_kernel_eqU_xr");
        block->instance->OPS_kernels[380].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "bounds_kernel_eqU_xr");
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
    int xdim0_bounds_kernel_eqU_xr = args[0].dat->size[0];
    int ydim0_bounds_kernel_eqU_xr = args[0].dat->size[1];
    int xdim1_bounds_kernel_eqU_xr = args[1].dat->size[0];
    int ydim1_bounds_kernel_eqU_xr = args[1].dat->size[1];
    int xdim2_bounds_kernel_eqU_xr = args[2].dat->size[0];
    int ydim2_bounds_kernel_eqU_xr = args[2].dat->size[1];
    int xdim3_bounds_kernel_eqU_xr = args[3].dat->size[0];
    int ydim3_bounds_kernel_eqU_xr = args[3].dat->size[1];
    int xdim4_bounds_kernel_eqU_xr = args[4].dat->size[0];
    int ydim4_bounds_kernel_eqU_xr = args[4].dat->size[1];
    int xdim5_bounds_kernel_eqU_xr = args[5].dat->size[0];
    int ydim5_bounds_kernel_eqU_xr = args[5].dat->size[1];
    int xdim6_bounds_kernel_eqU_xr = args[6].dat->size[0];
    int ydim6_bounds_kernel_eqU_xr = args[6].dat->size[1];
    int xdim7_bounds_kernel_eqU_xr = args[7].dat->size[0];
    int ydim7_bounds_kernel_eqU_xr = args[7].dat->size[1];
    int xdim8_bounds_kernel_eqU_xr = args[8].dat->size[0];
    int ydim8_bounds_kernel_eqU_xr = args[8].dat->size[1];
    int xdim9_bounds_kernel_eqU_xr = args[9].dat->size[0];
    int ydim9_bounds_kernel_eqU_xr = args[9].dat->size[1];
    int xdim10_bounds_kernel_eqU_xr = args[10].dat->size[0];
    int ydim10_bounds_kernel_eqU_xr = args[10].dat->size[1];
    int xdim11_bounds_kernel_eqU_xr = args[11].dat->size[0];
    int ydim11_bounds_kernel_eqU_xr = args[11].dat->size[1];
    int xdim12_bounds_kernel_eqU_xr = args[12].dat->size[0];
    int ydim12_bounds_kernel_eqU_xr = args[12].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ bcl1xr_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ bcl2xr_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ bcl3xr_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ bcl4xr_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ strdxr_p = (double *)(args[5].data) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ strexr_p = (double *)(args[6].data) + base6 - 1; // Subtracting 1 to convert to C-style

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ struxr_p = (double *)(args[7].data) + base7 - 1; // Subtracting 1 to convert to C-style

    int base8 = getDatBaseFromOpsArg3D(&args[8], start_indx, 1);
    double * __restrict__ strvxr_p = (double *)(args[8].data) + base8 - 1; // Subtracting 1 to convert to C-style

    int base9 = getDatBaseFromOpsArg3D(&args[9], start_indx, 1);
    double * __restrict__ strwxr_p = (double *)(args[9].data) + base9 - 1; // Subtracting 1 to convert to C-style

    int base10 = getDatBaseFromOpsArg3D(&args[10], start_indx, 1);
    double * __restrict__ ova2xr_p = (double *)(args[10].data) + base10 - 1; // Subtracting 1 to convert to C-style

    int base11 = getDatBaseFromOpsArg3D(&args[11], start_indx, 1);
    double * __restrict__ ovgmxr_p = (double *)(args[11].data) + base11 - 1; // Subtracting 1 to convert to C-style

    int base12 = getDatBaseFromOpsArg3D(&args[12], start_indx, 1);
    double * __restrict__ acouxr_p = (double *)(args[12].data) + base12 - 1; // Subtracting 1 to convert to C-style

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 13);
    ops_halo_exchanges(args, 13, range);
    ops_H_D_exchanges_host(args, 13);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[380].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                 ACC<double> erhs(xdim0_bounds_kernel_eqU_xr, ydim0_bounds_kernel_eqU_xr, erhs_p + (n_x * 1) + (n_y * xdim0_bounds_kernel_eqU_xr * 1) + (n_z * xdim0_bounds_kernel_eqU_xr * ydim0_bounds_kernel_eqU_xr * 1));

                const  ACC<double> bcl1xr(xdim1_bounds_kernel_eqU_xr, ydim1_bounds_kernel_eqU_xr, bcl1xr_p + (n_x * 0) + (n_y * xdim1_bounds_kernel_eqU_xr * 1) + (n_z * xdim1_bounds_kernel_eqU_xr * ydim1_bounds_kernel_eqU_xr * 1));

                const  ACC<double> bcl2xr(xdim2_bounds_kernel_eqU_xr, ydim2_bounds_kernel_eqU_xr, bcl2xr_p + (n_x * 0) + (n_y * xdim2_bounds_kernel_eqU_xr * 1) + (n_z * xdim2_bounds_kernel_eqU_xr * ydim2_bounds_kernel_eqU_xr * 1));

                const  ACC<double> bcl3xr(xdim3_bounds_kernel_eqU_xr, ydim3_bounds_kernel_eqU_xr, bcl3xr_p + (n_x * 0) + (n_y * xdim3_bounds_kernel_eqU_xr * 1) + (n_z * xdim3_bounds_kernel_eqU_xr * ydim3_bounds_kernel_eqU_xr * 1));

                const  ACC<double> bcl4xr(xdim4_bounds_kernel_eqU_xr, ydim4_bounds_kernel_eqU_xr, bcl4xr_p + (n_x * 0) + (n_y * xdim4_bounds_kernel_eqU_xr * 1) + (n_z * xdim4_bounds_kernel_eqU_xr * ydim4_bounds_kernel_eqU_xr * 1));

                const  ACC<double> strdxr(xdim5_bounds_kernel_eqU_xr, ydim5_bounds_kernel_eqU_xr, strdxr_p + (n_x * 0) + (n_y * xdim5_bounds_kernel_eqU_xr * 1) + (n_z * xdim5_bounds_kernel_eqU_xr * ydim5_bounds_kernel_eqU_xr * 1));

                const  ACC<double> strexr(xdim6_bounds_kernel_eqU_xr, ydim6_bounds_kernel_eqU_xr, strexr_p + (n_x * 0) + (n_y * xdim6_bounds_kernel_eqU_xr * 1) + (n_z * xdim6_bounds_kernel_eqU_xr * ydim6_bounds_kernel_eqU_xr * 1));

                const  ACC<double> struxr(xdim7_bounds_kernel_eqU_xr, ydim7_bounds_kernel_eqU_xr, struxr_p + (n_x * 0) + (n_y * xdim7_bounds_kernel_eqU_xr * 1) + (n_z * xdim7_bounds_kernel_eqU_xr * ydim7_bounds_kernel_eqU_xr * 1));

                const  ACC<double> strvxr(xdim8_bounds_kernel_eqU_xr, ydim8_bounds_kernel_eqU_xr, strvxr_p + (n_x * 0) + (n_y * xdim8_bounds_kernel_eqU_xr * 1) + (n_z * xdim8_bounds_kernel_eqU_xr * ydim8_bounds_kernel_eqU_xr * 1));

                const  ACC<double> strwxr(xdim9_bounds_kernel_eqU_xr, ydim9_bounds_kernel_eqU_xr, strwxr_p + (n_x * 0) + (n_y * xdim9_bounds_kernel_eqU_xr * 1) + (n_z * xdim9_bounds_kernel_eqU_xr * ydim9_bounds_kernel_eqU_xr * 1));

                const  ACC<double> ova2xr(xdim10_bounds_kernel_eqU_xr, ydim10_bounds_kernel_eqU_xr, ova2xr_p + (n_x * 0) + (n_y * xdim10_bounds_kernel_eqU_xr * 1) + (n_z * xdim10_bounds_kernel_eqU_xr * ydim10_bounds_kernel_eqU_xr * 1));

                const  ACC<double> ovgmxr(xdim11_bounds_kernel_eqU_xr, ydim11_bounds_kernel_eqU_xr, ovgmxr_p + (n_x * 0) + (n_y * xdim11_bounds_kernel_eqU_xr * 1) + (n_z * xdim11_bounds_kernel_eqU_xr * ydim11_bounds_kernel_eqU_xr * 1));

                const  ACC<double> acouxr(xdim12_bounds_kernel_eqU_xr, ydim12_bounds_kernel_eqU_xr, acouxr_p + (n_x * 0) + (n_y * xdim12_bounds_kernel_eqU_xr * 1) + (n_z * xdim12_bounds_kernel_eqU_xr * ydim12_bounds_kernel_eqU_xr * 1));

    erhs(0, 0, 0) = erhs(0, 0, 0) - bcl1xr(0, 0, 0) * (ova2xr(0, 0, 0) * strexr(0, 0, 0) - struxr(0, 0, 0) / acouxr(0, 0, 0) + ovgmxr(0, 0, 0)) - bcl2xr(0, 0, 0) * strexr(0, 0, 0) - bcl3xr(0, 0, 0) * strdxr(0, 0, 0) * strvxr(0, 0, 0) - bcl4xr(0, 0, 0) * strdxr(0, 0, 0) * strwxr(0, 0, 0);

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[380].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 13);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[380].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg8);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg9);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg10);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg11);
        block->instance->OPS_kernels[380].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg12);
    }
}

#ifdef OPS_LAZY
extern "C"
void bounds_kernel_eqU_xr_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bounds_kernel_eqU_xr", args, 13, 380, dim, 0, range, block, bounds_kernel_eqU_xr_execute);
}
#endif

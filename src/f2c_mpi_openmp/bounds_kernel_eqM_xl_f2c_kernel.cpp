// Auto-generated at 2026-04-28 18:43:29.882626 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bounds_kernel_eqM_xl_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bounds_kernel_eqM_xl_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[9];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
    args[4] = desc->args[4];
    args[5] = desc->args[5];
    args[6] = desc->args[6];
    args[7] = desc->args[7];
    args[8] = desc->args[8];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 9, range, 363)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 363, "bounds_kernel_eqM_xl");
        block->instance->OPS_kernels[363].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "bounds_kernel_eqM_xl");
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
    int xdim0_bounds_kernel_eqM_xl = args[0].dat->size[0];
    int ydim0_bounds_kernel_eqM_xl = args[0].dat->size[1];
    int xdim1_bounds_kernel_eqM_xl = args[1].dat->size[0];
    int ydim1_bounds_kernel_eqM_xl = args[1].dat->size[1];
    int xdim2_bounds_kernel_eqM_xl = args[2].dat->size[0];
    int ydim2_bounds_kernel_eqM_xl = args[2].dat->size[1];
    int xdim3_bounds_kernel_eqM_xl = args[3].dat->size[0];
    int ydim3_bounds_kernel_eqM_xl = args[3].dat->size[1];
    int xdim4_bounds_kernel_eqM_xl = args[4].dat->size[0];
    int ydim4_bounds_kernel_eqM_xl = args[4].dat->size[1];
    int xdim5_bounds_kernel_eqM_xl = args[5].dat->size[0];
    int ydim5_bounds_kernel_eqM_xl = args[5].dat->size[1];
    int xdim6_bounds_kernel_eqM_xl = args[6].dat->size[0];
    int ydim6_bounds_kernel_eqM_xl = args[6].dat->size[1];
    int xdim7_bounds_kernel_eqM_xl = args[7].dat->size[0];
    int ydim7_bounds_kernel_eqM_xl = args[7].dat->size[1];
    int xdim8_bounds_kernel_eqM_xl = args[8].dat->size[0];
    int ydim8_bounds_kernel_eqM_xl = args[8].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ bcl1xl_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ bcl2xl_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ bcl3xl_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ bcl4xl_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ bcl5xl_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ strdxl_p = (double *)(args[5].data) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ acouxl_p = (double *)(args[6].data) + base6 - 1; // Subtracting 1 to convert to C-style

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ struxl_p = (double *)(args[7].data) + base7 - 1; // Subtracting 1 to convert to C-style

    int base8 = getDatBaseFromOpsArg3D(&args[8], start_indx, 1);
    double * __restrict__ ova2xl_p = (double *)(args[8].data) + base8 - 1; // Subtracting 1 to convert to C-style

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 9);
    ops_halo_exchanges(args, 9, range);
    ops_H_D_exchanges_host(args, 9);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[363].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                 ACC<double> bcl1xl(xdim0_bounds_kernel_eqM_xl, ydim0_bounds_kernel_eqM_xl, bcl1xl_p + (n_x * 0) + (n_y * xdim0_bounds_kernel_eqM_xl * 1) + (n_z * xdim0_bounds_kernel_eqM_xl * ydim0_bounds_kernel_eqM_xl * 1));

                 ACC<double> bcl2xl(xdim1_bounds_kernel_eqM_xl, ydim1_bounds_kernel_eqM_xl, bcl2xl_p + (n_x * 0) + (n_y * xdim1_bounds_kernel_eqM_xl * 1) + (n_z * xdim1_bounds_kernel_eqM_xl * ydim1_bounds_kernel_eqM_xl * 1));

                 ACC<double> bcl3xl(xdim2_bounds_kernel_eqM_xl, ydim2_bounds_kernel_eqM_xl, bcl3xl_p + (n_x * 0) + (n_y * xdim2_bounds_kernel_eqM_xl * 1) + (n_z * xdim2_bounds_kernel_eqM_xl * ydim2_bounds_kernel_eqM_xl * 1));

                 ACC<double> bcl4xl(xdim3_bounds_kernel_eqM_xl, ydim3_bounds_kernel_eqM_xl, bcl4xl_p + (n_x * 0) + (n_y * xdim3_bounds_kernel_eqM_xl * 1) + (n_z * xdim3_bounds_kernel_eqM_xl * ydim3_bounds_kernel_eqM_xl * 1));

                 ACC<double> bcl5xl(xdim4_bounds_kernel_eqM_xl, ydim4_bounds_kernel_eqM_xl, bcl5xl_p + (n_x * 0) + (n_y * xdim4_bounds_kernel_eqM_xl * 1) + (n_z * xdim4_bounds_kernel_eqM_xl * ydim4_bounds_kernel_eqM_xl * 1));

                const  ACC<double> strdxl(xdim5_bounds_kernel_eqM_xl, ydim5_bounds_kernel_eqM_xl, strdxl_p + (n_x * 0) + (n_y * xdim5_bounds_kernel_eqM_xl * 1) + (n_z * xdim5_bounds_kernel_eqM_xl * ydim5_bounds_kernel_eqM_xl * 1));

                const  ACC<double> acouxl(xdim6_bounds_kernel_eqM_xl, ydim6_bounds_kernel_eqM_xl, acouxl_p + (n_x * 0) + (n_y * xdim6_bounds_kernel_eqM_xl * 1) + (n_z * xdim6_bounds_kernel_eqM_xl * ydim6_bounds_kernel_eqM_xl * 1));

                const  ACC<double> struxl(xdim7_bounds_kernel_eqM_xl, ydim7_bounds_kernel_eqM_xl, struxl_p + (n_x * 0) + (n_y * xdim7_bounds_kernel_eqM_xl * 1) + (n_z * xdim7_bounds_kernel_eqM_xl * ydim7_bounds_kernel_eqM_xl * 1));

                const  ACC<double> ova2xl(xdim8_bounds_kernel_eqM_xl, ydim8_bounds_kernel_eqM_xl, ova2xl_p + (n_x * 0) + (n_y * xdim8_bounds_kernel_eqM_xl * 1) + (n_z * xdim8_bounds_kernel_eqM_xl * ydim8_bounds_kernel_eqM_xl * 1));

    double fornow;

    fornow = strdxl(0, 0, 0) * acouxl(0, 0, 0) * bcl1xl(0, 0, 0);
    bcl1xl(0, 0, 0) = 0.5 * (struxl(0, 0, 0) - acouxl(0, 0, 0)) * (bcl5xl(0, 0, 0) - fornow);
    bcl2xl(0, 0, 0) = struxl(0, 0, 0) * (bcl2xl(0, 0, 0) - bcl5xl(0, 0, 0) * ova2xl(0, 0, 0));
    bcl3xl(0, 0, 0) = struxl(0, 0, 0) * bcl3xl(0, 0, 0);
    bcl4xl(0, 0, 0) = struxl(0, 0, 0) * bcl4xl(0, 0, 0);
    bcl5xl(0, 0, 0) = 0.5 * (struxl(0, 0, 0) + acouxl(0, 0, 0)) * (bcl5xl(0, 0, 0) + fornow);
    bcl3xl(0, 0, 0) = 0.0;
    bcl4xl(0, 0, 0) = 0.0;
    bcl5xl(0, 0, 0) = bcl1xl(0, 0, 0) - bcl5xl(0, 0, 0);
    bcl2xl(0, 0, 0) = 0.0;

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[363].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 9);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
    ops_set_halo_dirtybit3(&args[3], range);
    ops_set_halo_dirtybit3(&args[4], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[363].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[363].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[363].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[363].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[363].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[363].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[363].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[363].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[363].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
        block->instance->OPS_kernels[363].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg8);
    }
}

#ifdef OPS_LAZY
extern "C"
void bounds_kernel_eqM_xl_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bounds_kernel_eqM_xl", args, 9, 363, dim, 0, range, block, bounds_kernel_eqM_xl_execute);
}
#endif

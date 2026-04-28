// Auto-generated at 2026-04-28 18:43:37.660129 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void temper_kernel_eqA_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void temper_kernel_eqA_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[6];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
    args[4] = desc->args[4];
    args[5] = desc->args[5];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 6, range, 558)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 558, "temper_kernel_eqA");
        block->instance->OPS_kernels[558].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "temper_kernel_eqA");
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
    int xdim0_temper_kernel_eqA = args[0].dat->size[0];
    int ydim0_temper_kernel_eqA = args[0].dat->size[1];
    int zdim0_temper_kernel_eqA = args[0].dat->size[2];
    int xdim1_temper_kernel_eqA = args[1].dat->size[0];
    int ydim1_temper_kernel_eqA = args[1].dat->size[1];
    int xdim2_temper_kernel_eqA = args[2].dat->size[0];
    int ydim2_temper_kernel_eqA = args[2].dat->size[1];
    int xdim3_temper_kernel_eqA = args[3].dat->size[0];
    int ydim3_temper_kernel_eqA = args[3].dat->size[1];
    int xdim4_temper_kernel_eqA = args[4].dat->size[0];
    int ydim4_temper_kernel_eqA = args[4].dat->size[1];
    int xdim5_temper_kernel_eqA = args[5].dat->size[0];
    int ydim5_temper_kernel_eqA = args[5].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int multi_d0 = getDatDimFromOpsArg(&args[0]);
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, multi_d0);
    double * __restrict__ tcoeff_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ urhs_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ vrhs_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ wrhs_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[5].data) + base5 - 1; // Subtracting 1 to convert to C-style

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 6);
    ops_halo_exchanges(args, 6, range);
    ops_H_D_exchanges_host(args, 6);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[558].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

#ifdef OPS_SOA
                 ACC<double> tcoeff(6,xdim0_temper_kernel_eqA, ydim0_temper_kernel_eqA, zdim0_temper_kernel_eqA, tcoeff_p + (n_x * 1) + (n_y * xdim0_temper_kernel_eqA * 1) + (n_z * xdim0_temper_kernel_eqA * ydim0_temper_kernel_eqA * 1));
#else
                 ACC<double> tcoeff(6,xdim0_temper_kernel_eqA, ydim0_temper_kernel_eqA, zdim0_temper_kernel_eqA, tcoeff_p + 6 * ((n_x * 1) + (n_y * xdim0_temper_kernel_eqA * 1) + (n_z * xdim0_temper_kernel_eqA * ydim0_temper_kernel_eqA * 1)));
#endif

                const  ACC<double> drhs(xdim1_temper_kernel_eqA, ydim1_temper_kernel_eqA, drhs_p + (n_x * 1) + (n_y * xdim1_temper_kernel_eqA * 1) + (n_z * xdim1_temper_kernel_eqA * ydim1_temper_kernel_eqA * 1));

                const  ACC<double> urhs(xdim2_temper_kernel_eqA, ydim2_temper_kernel_eqA, urhs_p + (n_x * 1) + (n_y * xdim2_temper_kernel_eqA * 1) + (n_z * xdim2_temper_kernel_eqA * ydim2_temper_kernel_eqA * 1));

                const  ACC<double> vrhs(xdim3_temper_kernel_eqA, ydim3_temper_kernel_eqA, vrhs_p + (n_x * 1) + (n_y * xdim3_temper_kernel_eqA * 1) + (n_z * xdim3_temper_kernel_eqA * ydim3_temper_kernel_eqA * 1));

                const  ACC<double> wrhs(xdim4_temper_kernel_eqA, ydim4_temper_kernel_eqA, wrhs_p + (n_x * 1) + (n_y * xdim4_temper_kernel_eqA * 1) + (n_z * xdim4_temper_kernel_eqA * ydim4_temper_kernel_eqA * 1));

                const  ACC<double> erhs(xdim5_temper_kernel_eqA, ydim5_temper_kernel_eqA, erhs_p + (n_x * 1) + (n_y * xdim5_temper_kernel_eqA * 1) + (n_z * xdim5_temper_kernel_eqA * ydim5_temper_kernel_eqA * 1));

    double ukuk;
    int icp;

    ukuk = (urhs(0, 0, 0) * urhs(0, 0, 0) + vrhs(0, 0, 0) * vrhs(0, 0, 0) + wrhs(0, 0, 0) * wrhs(0, 0, 0)) / drhs(0, 0, 0);
    tcoeff(1-1, 0, 0, 0) = 0.5 * ukuk - erhs(0, 0, 0);
    for (icp = 1; icp <= nctmax; ++icp) {
        tcoeff(icp + 1-1, 0, 0, 0) = 0.0;
    }

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[558].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 6);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[558].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[558].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[558].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[558].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[558].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[558].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[558].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
    }
}

#ifdef OPS_LAZY
extern "C"
void temper_kernel_eqA_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("temper_kernel_eqA", args, 6, 558, dim, 0, range, block, temper_kernel_eqA_execute);
}
#endif

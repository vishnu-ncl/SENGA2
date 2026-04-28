// Auto-generated at 2026-04-28 18:43:37.588746 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void turbin_kernel_eqA_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void turbin_kernel_eqA_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 7, range, 555)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 555, "turbin_kernel_eqA");
        block->instance->OPS_kernels[555].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "turbin_kernel_eqA");
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
    int xdim0_turbin_kernel_eqA = args[0].dat->size[0];
    int ydim0_turbin_kernel_eqA = args[0].dat->size[1];
    int xdim1_turbin_kernel_eqA = args[1].dat->size[0];
    int ydim1_turbin_kernel_eqA = args[1].dat->size[1];
    int xdim2_turbin_kernel_eqA = args[2].dat->size[0];
    int ydim2_turbin_kernel_eqA = args[2].dat->size[1];
    int xdim3_turbin_kernel_eqA = args[3].dat->size[0];
    int ydim3_turbin_kernel_eqA = args[3].dat->size[1];
    int xdim4_turbin_kernel_eqA = args[4].dat->size[0];
    int ydim4_turbin_kernel_eqA = args[4].dat->size[1];
    int xdim5_turbin_kernel_eqA = args[5].dat->size[0];
    int ydim5_turbin_kernel_eqA = args[5].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ urun_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ utmp_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ vrun_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ vtmp_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ wrun_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ wtmp_p = (double *)(args[5].data) + base5 - 1; // Subtracting 1 to convert to C-style

#ifdef OPS_MPI
    double * __restrict__ p_a6 = (double *)(((ops_reduction)args[6].data)->data + ((ops_reduction)args[6].data)->size * block->index);
#else //OPS_MPI
    double * __restrict__ p_a6 = (double *)((ops_reduction)args[6].data)->data;
#endif //OPS_MPI

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 7);
    ops_halo_exchanges(args, 7, range);
    ops_H_D_exchanges_host(args, 7);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[555].mpi_time += __t2 - __t1;
    }

    double p_a6_0 = p_a6[0];

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                const  ACC<double> urun(xdim0_turbin_kernel_eqA, ydim0_turbin_kernel_eqA, urun_p + (n_x * 1) + (n_y * xdim0_turbin_kernel_eqA * 1) + (n_z * xdim0_turbin_kernel_eqA * ydim0_turbin_kernel_eqA * 1));

                const  ACC<double> utmp(xdim1_turbin_kernel_eqA, ydim1_turbin_kernel_eqA, utmp_p + (n_x * 1) + (n_y * xdim1_turbin_kernel_eqA * 1) + (n_z * xdim1_turbin_kernel_eqA * ydim1_turbin_kernel_eqA * 1));

                const  ACC<double> vrun(xdim2_turbin_kernel_eqA, ydim2_turbin_kernel_eqA, vrun_p + (n_x * 1) + (n_y * xdim2_turbin_kernel_eqA * 1) + (n_z * xdim2_turbin_kernel_eqA * ydim2_turbin_kernel_eqA * 1));

                const  ACC<double> vtmp(xdim3_turbin_kernel_eqA, ydim3_turbin_kernel_eqA, vtmp_p + (n_x * 1) + (n_y * xdim3_turbin_kernel_eqA * 1) + (n_z * xdim3_turbin_kernel_eqA * ydim3_turbin_kernel_eqA * 1));

                const  ACC<double> wrun(xdim4_turbin_kernel_eqA, ydim4_turbin_kernel_eqA, wrun_p + (n_x * 1) + (n_y * xdim4_turbin_kernel_eqA * 1) + (n_z * xdim4_turbin_kernel_eqA * ydim4_turbin_kernel_eqA * 1));

                const  ACC<double> wtmp(xdim5_turbin_kernel_eqA, ydim5_turbin_kernel_eqA, wtmp_p + (n_x * 1) + (n_y * xdim5_turbin_kernel_eqA * 1) + (n_z * xdim5_turbin_kernel_eqA * ydim5_turbin_kernel_eqA * 1));

                double tket[1];
                tket[0] = ZERO_double;

    *tket = *tket + urun(0, 0, 0) * urun(0, 0, 0) + utmp(0, 0, 0) * utmp(0, 0, 0) + vrun(0, 0, 0) * vrun(0, 0, 0) + vtmp(0, 0, 0) * vtmp(0, 0, 0) + wrun(0, 0, 0) * wrun(0, 0, 0) + wtmp(0, 0, 0) * wtmp(0, 0, 0);

                            p_a6_0 += tket[0];

            }
      }

    }

    p_a6[0] = p_a6_0;    

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[555].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 7);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[555].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[555].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
    }
}

#ifdef OPS_LAZY
extern "C"
void turbin_kernel_eqA_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("turbin_kernel_eqA", args, 7, 555, dim, 0, range, block, turbin_kernel_eqA_execute);
}
#endif

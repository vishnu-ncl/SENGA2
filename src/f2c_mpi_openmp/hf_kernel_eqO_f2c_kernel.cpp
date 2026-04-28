// Auto-generated at 2026-04-28 18:43:26.545542 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void hf_kernel_eqO_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void hf_kernel_eqO_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[10];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 10, range, 259)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 259, "hf_kernel_eqO");
        block->instance->OPS_kernels[259].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "hf_kernel_eqO");
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
    int xdim0_hf_kernel_eqO = args[0].dat->size[0];
    int ydim0_hf_kernel_eqO = args[0].dat->size[1];
    int xdim1_hf_kernel_eqO = args[1].dat->size[0];
    int ydim1_hf_kernel_eqO = args[1].dat->size[1];
    int xdim2_hf_kernel_eqO = args[2].dat->size[0];
    int ydim2_hf_kernel_eqO = args[2].dat->size[1];
    int xdim3_hf_kernel_eqO = args[3].dat->size[0];
    int ydim3_hf_kernel_eqO = args[3].dat->size[1];
    int xdim4_hf_kernel_eqO = args[4].dat->size[0];
    int ydim4_hf_kernel_eqO = args[4].dat->size[1];
    int xdim5_hf_kernel_eqO = args[5].dat->size[0];
    int ydim5_hf_kernel_eqO = args[5].dat->size[1];
    int xdim6_hf_kernel_eqO = args[6].dat->size[0];
    int ydim6_hf_kernel_eqO = args[6].dat->size[1];
    int xdim7_hf_kernel_eqO = args[7].dat->size[0];
    int ydim7_hf_kernel_eqO = args[7].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ rate_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ vtmp_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ trun_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ tdrmix_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ utmp_p = (double *)(args[5].data) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ store7_p = (double *)(args[6].data) + base6 - 1; // Subtracting 1 to convert to C-style

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ td1y_p = (double *)(args[7].data) + base7 - 1; // Subtracting 1 to convert to C-style

    double * __restrict__  deltay = (double *)args[8].data;
    double * __restrict__  rgspec_ispec = (double *)args[9].data;

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 10);
    ops_halo_exchanges(args, 10, range);
    ops_H_D_exchanges_host(args, 10);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[259].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                 ACC<double> rate(xdim0_hf_kernel_eqO, ydim0_hf_kernel_eqO, rate_p + (n_x * 1) + (n_y * xdim0_hf_kernel_eqO * 1) + (n_z * xdim0_hf_kernel_eqO * ydim0_hf_kernel_eqO * 1));

                 ACC<double> erhs(xdim1_hf_kernel_eqO, ydim1_hf_kernel_eqO, erhs_p + (n_x * 1) + (n_y * xdim1_hf_kernel_eqO * 1) + (n_z * xdim1_hf_kernel_eqO * ydim1_hf_kernel_eqO * 1));

                 ACC<double> vtmp(xdim2_hf_kernel_eqO, ydim2_hf_kernel_eqO, vtmp_p + (n_x * 1) + (n_y * xdim2_hf_kernel_eqO * 1) + (n_z * xdim2_hf_kernel_eqO * ydim2_hf_kernel_eqO * 1));

                const  ACC<double> trun(xdim3_hf_kernel_eqO, ydim3_hf_kernel_eqO, trun_p + (n_x * 1) + (n_y * xdim3_hf_kernel_eqO * 1) + (n_z * xdim3_hf_kernel_eqO * ydim3_hf_kernel_eqO * 1));

                const  ACC<double> tdrmix(xdim4_hf_kernel_eqO, ydim4_hf_kernel_eqO, tdrmix_p + (n_x * 1) + (n_y * xdim4_hf_kernel_eqO * 1) + (n_z * xdim4_hf_kernel_eqO * ydim4_hf_kernel_eqO * 1));

                const  ACC<double> utmp(xdim5_hf_kernel_eqO, ydim5_hf_kernel_eqO, utmp_p + (n_x * 1) + (n_y * xdim5_hf_kernel_eqO * 1) + (n_z * xdim5_hf_kernel_eqO * ydim5_hf_kernel_eqO * 1));

                const  ACC<double> store7(xdim6_hf_kernel_eqO, ydim6_hf_kernel_eqO, store7_p + (n_x * 1) + (n_y * xdim6_hf_kernel_eqO * 1) + (n_z * xdim6_hf_kernel_eqO * ydim6_hf_kernel_eqO * 1));

                const  ACC<double> td1y(xdim7_hf_kernel_eqO, ydim7_hf_kernel_eqO, td1y_p + (n_x * 1) + (n_y * xdim7_hf_kernel_eqO * 1) + (n_z * xdim7_hf_kernel_eqO * ydim7_hf_kernel_eqO * 1));

    double fornow;
    double fornow0;
    double fornow1;
    double fornow2;
    double fornow3;
    double fornow4;

    fornow0 = store7(0, 0, 0) * td1y(0, 0, 0);
    fornow1 = store7(0, 1, 0) * td1y(0, 1, 0);
    fornow2 = store7(0, 2, 0) * td1y(0, 2, 0);
    fornow3 = store7(0, 3, 0) * td1y(0, 3, 0);
    fornow4 = store7(0, 4, 0) * td1y(0, 4, 0);
    fornow = 0.0;
    fornow = ((-25.0 / 12.0) * fornow0 + (4.0) * fornow1 - (3.0) * fornow2 + (4.0 / 3.0) * fornow3 - (1.0 / 4.0) * fornow4) / deltay[0];
    rate(0, 0, 0) = rate(0, 0, 0) + fornow;
    vtmp(0, 0, 0) = vtmp(0, 0, 0) + fornow;
    erhs(0, 0, 0) = erhs(0, 0, 0) + fornow * (utmp(0, 0, 0) + rgspec_ispec[0] * trun(0, 0, 0) * tdrmix(0, 0, 0));

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[259].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 10);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[259].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[259].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[259].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[259].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[259].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[259].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[259].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[259].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[259].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
    }
}

#ifdef OPS_LAZY
extern "C"
void hf_kernel_eqO_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("hf_kernel_eqO", args, 10, 259, dim, 0, range, block, hf_kernel_eqO_execute);
}
#endif

// Auto-generated at 2026-04-28 18:43:26.736779 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void hf_kernel_eqG_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void hf_kernel_eqG_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 7, range, 263)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 263, "hf_kernel_eqG");
        block->instance->OPS_kernels[263].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "hf_kernel_eqG");
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
    int xdim0_hf_kernel_eqG = args[0].dat->size[0];
    int ydim0_hf_kernel_eqG = args[0].dat->size[1];
    int xdim1_hf_kernel_eqG = args[1].dat->size[0];
    int ydim1_hf_kernel_eqG = args[1].dat->size[1];
    int xdim2_hf_kernel_eqG = args[2].dat->size[0];
    int ydim2_hf_kernel_eqG = args[2].dat->size[1];
    int xdim3_hf_kernel_eqG = args[3].dat->size[0];
    int ydim3_hf_kernel_eqG = args[3].dat->size[1];
    int xdim4_hf_kernel_eqG = args[4].dat->size[0];
    int ydim4_hf_kernel_eqG = args[4].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ trun_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ tdrmix_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ store7_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ td1x_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

    double * __restrict__  rgspec_ispec = (double *)args[5].data;
    double * __restrict__  deltax = (double *)args[6].data;

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
        block->instance->OPS_kernels[263].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                 ACC<double> erhs(xdim0_hf_kernel_eqG, ydim0_hf_kernel_eqG, erhs_p + (n_x * 1) + (n_y * xdim0_hf_kernel_eqG * 1) + (n_z * xdim0_hf_kernel_eqG * ydim0_hf_kernel_eqG * 1));

                const  ACC<double> trun(xdim1_hf_kernel_eqG, ydim1_hf_kernel_eqG, trun_p + (n_x * 1) + (n_y * xdim1_hf_kernel_eqG * 1) + (n_z * xdim1_hf_kernel_eqG * ydim1_hf_kernel_eqG * 1));

                const  ACC<double> tdrmix(xdim2_hf_kernel_eqG, ydim2_hf_kernel_eqG, tdrmix_p + (n_x * 1) + (n_y * xdim2_hf_kernel_eqG * 1) + (n_z * xdim2_hf_kernel_eqG * ydim2_hf_kernel_eqG * 1));

                const  ACC<double> store7(xdim3_hf_kernel_eqG, ydim3_hf_kernel_eqG, store7_p + (n_x * 1) + (n_y * xdim3_hf_kernel_eqG * 1) + (n_z * xdim3_hf_kernel_eqG * ydim3_hf_kernel_eqG * 1));

                const  ACC<double> td1x(xdim4_hf_kernel_eqG, ydim4_hf_kernel_eqG, td1x_p + (n_x * 1) + (n_y * xdim4_hf_kernel_eqG * 1) + (n_z * xdim4_hf_kernel_eqG * ydim4_hf_kernel_eqG * 1));

    double combo1;
    double fornow;
    double fornow0;
    double fornow1;
    double fornow2;
    double fornow3;
    double fornow4;

    fornow = 0.0;
    combo1 = trun(0, 0, 0) * tdrmix(0, 0, 0);
    fornow0 = combo1 * store7(0, 0, 0) * td1x(0, 0, 0);
    combo1 = trun(1, 0, 0) * tdrmix(1, 0, 0);
    fornow1 = combo1 * store7(1, 0, 0) * td1x(1, 0, 0);
    combo1 = trun(2, 0, 0) * tdrmix(2, 0, 0);
    fornow2 = combo1 * store7(2, 0, 0) * td1x(2, 0, 0);
    combo1 = trun(3, 0, 0) * tdrmix(3, 0, 0);
    fornow3 = combo1 * store7(3, 0, 0) * td1x(3, 0, 0);
    combo1 = trun(4, 0, 0) * tdrmix(4, 0, 0);
    fornow4 = combo1 * store7(4, 0, 0) * td1x(4, 0, 0);
    fornow = ((-25.0 / 12.0) * fornow0 + (4.0) * fornow1 - (3.0) * fornow2 + (4.0 / 3.0) * fornow3 - (1.0 / 4.0) * fornow4) / deltax[0];
    erhs(0, 0, 0) = erhs(0, 0, 0) + rgspec_ispec[0] * fornow;

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[263].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 7);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[263].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[263].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[263].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[263].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[263].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[263].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
    }
}

#ifdef OPS_LAZY
extern "C"
void hf_kernel_eqG_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("hf_kernel_eqG", args, 7, 263, dim, 0, range, block, hf_kernel_eqG_execute);
}
#endif

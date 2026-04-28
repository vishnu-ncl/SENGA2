// Auto-generated at 2026-04-28 18:43:37.174149 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void maths_kernel_eqAX_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void maths_kernel_eqAX_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 9, range, 539)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 539, "maths_kernel_eqAX");
        block->instance->OPS_kernels[539].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "maths_kernel_eqAX");
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
    int xdim0_maths_kernel_eqAX = args[0].dat->size[0];
    int ydim0_maths_kernel_eqAX = args[0].dat->size[1];
    int xdim1_maths_kernel_eqAX = args[1].dat->size[0];
    int ydim1_maths_kernel_eqAX = args[1].dat->size[1];
    int xdim2_maths_kernel_eqAX = args[2].dat->size[0];
    int ydim2_maths_kernel_eqAX = args[2].dat->size[1];
    int xdim3_maths_kernel_eqAX = args[3].dat->size[0];
    int ydim3_maths_kernel_eqAX = args[3].dat->size[1];
    int xdim4_maths_kernel_eqAX = args[4].dat->size[0];
    int ydim4_maths_kernel_eqAX = args[4].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ out_arr1_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ out_arr2_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ out_arr3_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ in_arr1_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ in_arr2_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

    double * __restrict__  racnst = (double *)args[5].data;
    double * __restrict__  rncnst = (double *)args[6].data;
    double * __restrict__  reovrr = (double *)args[7].data;
    double * __restrict__  flcnst = (double *)args[8].data;

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
        block->instance->OPS_kernels[539].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                 ACC<double> out_arr1(xdim0_maths_kernel_eqAX, ydim0_maths_kernel_eqAX, out_arr1_p + (n_x * 1) + (n_y * xdim0_maths_kernel_eqAX * 1) + (n_z * xdim0_maths_kernel_eqAX * ydim0_maths_kernel_eqAX * 1));

                 ACC<double> out_arr2(xdim1_maths_kernel_eqAX, ydim1_maths_kernel_eqAX, out_arr2_p + (n_x * 1) + (n_y * xdim1_maths_kernel_eqAX * 1) + (n_z * xdim1_maths_kernel_eqAX * ydim1_maths_kernel_eqAX * 1));

                 ACC<double> out_arr3(xdim2_maths_kernel_eqAX, ydim2_maths_kernel_eqAX, out_arr3_p + (n_x * 1) + (n_y * xdim2_maths_kernel_eqAX * 1) + (n_z * xdim2_maths_kernel_eqAX * ydim2_maths_kernel_eqAX * 1));

                const  ACC<double> in_arr1(xdim3_maths_kernel_eqAX, ydim3_maths_kernel_eqAX, in_arr1_p + (n_x * 1) + (n_y * xdim3_maths_kernel_eqAX * 1) + (n_z * xdim3_maths_kernel_eqAX * ydim3_maths_kernel_eqAX * 1));

                const  ACC<double> in_arr2(xdim4_maths_kernel_eqAX, ydim4_maths_kernel_eqAX, in_arr2_p + (n_x * 1) + (n_y * xdim4_maths_kernel_eqAX * 1) + (n_z * xdim4_maths_kernel_eqAX * ydim4_maths_kernel_eqAX * 1));

    double preduc;
    double fornow;

    out_arr3(0, 0, 0) = racnst[0] + rncnst[0] * f2c::log(in_arr1(0, 0, 0)) - reovrr[0] / in_arr1(0, 0, 0);
    out_arr3(0, 0, 0) = f2c::exp(out_arr3(0, 0, 0));
    preduc = in_arr2(0, 0, 0) * out_arr3(0, 0, 0) / out_arr1(0, 0, 0);
    fornow = flcnst[0] * preduc / (1.0 + preduc);
    out_arr1(0, 0, 0) = out_arr1(0, 0, 0) * fornow;
    out_arr2(0, 0, 0) = f2c::log(out_arr1(0, 0, 0));

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[539].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 9);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[539].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[539].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[539].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[539].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[539].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[539].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
    }
}

#ifdef OPS_LAZY
extern "C"
void maths_kernel_eqAX_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("maths_kernel_eqAX", args, 9, 539, dim, 0, range, block, maths_kernel_eqAX_execute);
}
#endif

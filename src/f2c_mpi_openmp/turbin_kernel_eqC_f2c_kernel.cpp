// Auto-generated at 2026-04-28 18:43:37.625020 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void turbin_kernel_eqC_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void turbin_kernel_eqC_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[3];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
#endif

//  ======
//  Timing
//  ======
    double __t1, __t2, __c1, __c2;

    ops_arg arg0 = args[0];
    ops_arg arg1 = args[1];
    ops_arg arg2 = args[2];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 3, range, 557)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 557, "turbin_kernel_eqC");
        block->instance->OPS_kernels[557].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "turbin_kernel_eqC");
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
    int xdim0_turbin_kernel_eqC = args[0].dat->size[0];
    int ydim0_turbin_kernel_eqC = args[0].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ in_arr_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    double * __restrict__  val1 = (double *)args[1].data;
#ifdef OPS_MPI
    double * __restrict__ p_a2 = (double *)(((ops_reduction)args[2].data)->data + ((ops_reduction)args[2].data)->size * block->index);
#else //OPS_MPI
    double * __restrict__ p_a2 = (double *)((ops_reduction)args[2].data)->data;
#endif //OPS_MPI

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 3);
    ops_halo_exchanges(args, 3, range);
    ops_H_D_exchanges_host(args, 3);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[557].mpi_time += __t2 - __t1;
    }

    double p_a2_0 = p_a2[0];

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                const  ACC<double> in_arr(xdim0_turbin_kernel_eqC, ydim0_turbin_kernel_eqC, in_arr_p + (n_x * 1) + (n_y * xdim0_turbin_kernel_eqC * 1) + (n_z * xdim0_turbin_kernel_eqC * ydim0_turbin_kernel_eqC * 1));

                double reduction_var[1];
                reduction_var[0] = ZERO_double;

    double fornow;

    fornow = in_arr(0, 0, 0) - val1[0];
    *reduction_var = *reduction_var + fornow * fornow;

                            p_a2_0 += reduction_var[0];

            }
      }

    }

    p_a2[0] = p_a2_0;    

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[557].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 3);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[557].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[557].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
    }
}

#ifdef OPS_LAZY
extern "C"
void turbin_kernel_eqC_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("turbin_kernel_eqC", args, 3, 557, dim, 0, range, block, turbin_kernel_eqC_execute);
}
#endif

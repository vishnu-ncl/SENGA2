// Auto-generated at 2026-04-28 18:43:37.033315 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void radcal_kernel_meancoef_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void radcal_kernel_meancoef_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 6, range, 534)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 534, "radcal_kernel_meancoef");
        block->instance->OPS_kernels[534].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "radcal_kernel_meancoef");
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
    int xdim0_radcal_kernel_meancoef = args[0].dat->size[0];
    int ydim0_radcal_kernel_meancoef = args[0].dat->size[1];
    int xdim1_radcal_kernel_meancoef = args[1].dat->size[0];
    int ydim1_radcal_kernel_meancoef = args[1].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ store2_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ trun_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    double * __restrict__  akprad = (double *)args[2].data;
    int * __restrict__  nkprad = (int *)args[3].data;
    int * __restrict__  nkprm1 = (int *)args[4].data;
    int * __restrict__  jspec = (int *)args[5].data;

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
        block->instance->OPS_kernels[534].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                 ACC<double> store2(xdim0_radcal_kernel_meancoef, ydim0_radcal_kernel_meancoef, store2_p + (n_x * 1) + (n_y * xdim0_radcal_kernel_meancoef * 1) + (n_z * xdim0_radcal_kernel_meancoef * ydim0_radcal_kernel_meancoef * 1));

                const  ACC<double> trun(xdim1_radcal_kernel_meancoef, ydim1_radcal_kernel_meancoef, trun_p + (n_x * 1) + (n_y * xdim1_radcal_kernel_meancoef * 1) + (n_z * xdim1_radcal_kernel_meancoef * ydim1_radcal_kernel_meancoef * 1));

    double plspec;
    double fornow;
    int icp;

    fornow = trun(0, 0, 0);
    plspec = akprad[(nkprad[jspec[0]-1]-1)+(jspec[0]-1)*ncfrmx];
    for (icp = nkprm1[jspec[0]-1]; icp >= 1; icp -= 1) {
        plspec = plspec * fornow + akprad[(icp-1)+(jspec[0]-1)*ncfrmx];
    }
    store2(0, 0, 0) = plspec;

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[534].time += __t1 - __t2;
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
        block->instance->OPS_kernels[534].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[534].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[534].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void radcal_kernel_meancoef_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("radcal_kernel_meancoef", args, 6, 534, dim, 0, range, block, radcal_kernel_meancoef_execute);
}
#endif

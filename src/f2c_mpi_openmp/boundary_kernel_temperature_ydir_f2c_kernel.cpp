// Auto-generated at 2026-04-28 18:43:25.060329 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void boundary_kernel_temperature_ydir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void boundary_kernel_temperature_ydir_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[4];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
#endif

//  ======
//  Timing
//  ======
    double __t1, __t2, __c1, __c2;

    ops_arg arg0 = args[0];
    ops_arg arg1 = args[1];
    ops_arg arg2 = args[2];
    ops_arg arg3 = args[3];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 4, range, 208)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 208, "boundary_kernel_temperature_ydir");
        block->instance->OPS_kernels[208].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "boundary_kernel_temperature_ydir");
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
    int xdim0_boundary_kernel_temperature_ydir = args[0].dat->size[0];
    int ydim0_boundary_kernel_temperature_ydir = args[0].dat->size[1];
    int xdim1_boundary_kernel_temperature_ydir = args[1].dat->size[0];
    int ydim1_boundary_kernel_temperature_ydir = args[1].dat->size[1];
    int xdim2_boundary_kernel_temperature_ydir = args[2].dat->size[0];
    int ydim2_boundary_kernel_temperature_ydir = args[2].dat->size[1];
    int xdim3_boundary_kernel_temperature_ydir = args[3].dat->size[0];
    int ydim3_boundary_kernel_temperature_ydir = args[3].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ trun_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ store1_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ strty_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ bclty_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 4);
    ops_halo_exchanges(args, 4, range);
    ops_H_D_exchanges_host(args, 4);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[208].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                const  ACC<double> trun(xdim0_boundary_kernel_temperature_ydir, ydim0_boundary_kernel_temperature_ydir, trun_p + (n_x * 1) + (n_y * xdim0_boundary_kernel_temperature_ydir * 1) + (n_z * xdim0_boundary_kernel_temperature_ydir * ydim0_boundary_kernel_temperature_ydir * 1));

                const  ACC<double> store1(xdim1_boundary_kernel_temperature_ydir, ydim1_boundary_kernel_temperature_ydir, store1_p + (n_x * 1) + (n_y * xdim1_boundary_kernel_temperature_ydir * 1) + (n_z * xdim1_boundary_kernel_temperature_ydir * ydim1_boundary_kernel_temperature_ydir * 1));

                 ACC<double> strty(xdim2_boundary_kernel_temperature_ydir, ydim2_boundary_kernel_temperature_ydir, strty_p + (n_x * 1) + (n_y * xdim2_boundary_kernel_temperature_ydir * 0) + (n_z * xdim2_boundary_kernel_temperature_ydir * ydim2_boundary_kernel_temperature_ydir * 1));

                 ACC<double> bclty(xdim3_boundary_kernel_temperature_ydir, ydim3_boundary_kernel_temperature_ydir, bclty_p + (n_x * 1) + (n_y * xdim3_boundary_kernel_temperature_ydir * 0) + (n_z * xdim3_boundary_kernel_temperature_ydir * ydim3_boundary_kernel_temperature_ydir * 1));

    strty(0, 0, 0) = trun(0, 0, 0);
    bclty(0, 0, 0) = store1(0, 0, 0);

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[208].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 4);
    ops_set_halo_dirtybit3(&args[2], range);
    ops_set_halo_dirtybit3(&args[3], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[208].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[208].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[208].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[208].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[208].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
    }
}

#ifdef OPS_LAZY
extern "C"
void boundary_kernel_temperature_ydir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("boundary_kernel_temperature_ydir", args, 4, 208, dim, 0, range, block, boundary_kernel_temperature_ydir_execute);
}
#endif

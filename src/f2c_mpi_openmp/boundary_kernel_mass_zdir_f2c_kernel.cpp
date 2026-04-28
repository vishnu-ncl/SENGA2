// Auto-generated at 2026-04-28 18:43:25.636079 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void boundary_kernel_mass_zdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void boundary_kernel_mass_zdir_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 10, range, 231)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 231, "boundary_kernel_mass_zdir");
        block->instance->OPS_kernels[231].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "boundary_kernel_mass_zdir");
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
    int xdim0_boundary_kernel_mass_zdir = args[0].dat->size[0];
    int ydim0_boundary_kernel_mass_zdir = args[0].dat->size[1];
    int xdim1_boundary_kernel_mass_zdir = args[1].dat->size[0];
    int ydim1_boundary_kernel_mass_zdir = args[1].dat->size[1];
    int xdim2_boundary_kernel_mass_zdir = args[2].dat->size[0];
    int ydim2_boundary_kernel_mass_zdir = args[2].dat->size[1];
    int xdim3_boundary_kernel_mass_zdir = args[3].dat->size[0];
    int ydim3_boundary_kernel_mass_zdir = args[3].dat->size[1];
    int xdim4_boundary_kernel_mass_zdir = args[4].dat->size[0];
    int ydim4_boundary_kernel_mass_zdir = args[4].dat->size[1];
    int xdim5_boundary_kernel_mass_zdir = args[5].dat->size[0];
    int ydim5_boundary_kernel_mass_zdir = args[5].dat->size[1];
    int xdim6_boundary_kernel_mass_zdir = args[6].dat->size[0];
    int ydim6_boundary_kernel_mass_zdir = args[6].dat->size[1];
    int xdim7_boundary_kernel_mass_zdir = args[7].dat->size[0];
    int ydim7_boundary_kernel_mass_zdir = args[7].dat->size[1];
    int xdim8_boundary_kernel_mass_zdir = args[8].dat->size[0];
    int ydim8_boundary_kernel_mass_zdir = args[8].dat->size[1];
    int xdim9_boundary_kernel_mass_zdir = args[9].dat->size[0];
    int ydim9_boundary_kernel_mass_zdir = args[9].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ yrhs_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ store1_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ store2_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ store3_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ urhs_p = (double *)(args[5].data) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ vrhs_p = (double *)(args[6].data) + base6 - 1; // Subtracting 1 to convert to C-style

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ stryz_p = (double *)(args[7].data) + base7 - 1; // Subtracting 1 to convert to C-style

    int base8 = getDatBaseFromOpsArg3D(&args[8], start_indx, 1);
    double * __restrict__ bclyz_p = (double *)(args[8].data) + base8 - 1; // Subtracting 1 to convert to C-style

    int base9 = getDatBaseFromOpsArg3D(&args[9], start_indx, 1);
    double * __restrict__ t6bz_p = (double *)(args[9].data) + base9 - 1; // Subtracting 1 to convert to C-style

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
        block->instance->OPS_kernels[231].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                const  ACC<double> yrhs(xdim0_boundary_kernel_mass_zdir, ydim0_boundary_kernel_mass_zdir, yrhs_p + (n_x * 1) + (n_y * xdim0_boundary_kernel_mass_zdir * 1) + (n_z * xdim0_boundary_kernel_mass_zdir * ydim0_boundary_kernel_mass_zdir * 1));

                const  ACC<double> store1(xdim1_boundary_kernel_mass_zdir, ydim1_boundary_kernel_mass_zdir, store1_p + (n_x * 1) + (n_y * xdim1_boundary_kernel_mass_zdir * 1) + (n_z * xdim1_boundary_kernel_mass_zdir * ydim1_boundary_kernel_mass_zdir * 1));

                const  ACC<double> store2(xdim2_boundary_kernel_mass_zdir, ydim2_boundary_kernel_mass_zdir, store2_p + (n_x * 1) + (n_y * xdim2_boundary_kernel_mass_zdir * 1) + (n_z * xdim2_boundary_kernel_mass_zdir * ydim2_boundary_kernel_mass_zdir * 1));

                const  ACC<double> store3(xdim3_boundary_kernel_mass_zdir, ydim3_boundary_kernel_mass_zdir, store3_p + (n_x * 1) + (n_y * xdim3_boundary_kernel_mass_zdir * 1) + (n_z * xdim3_boundary_kernel_mass_zdir * ydim3_boundary_kernel_mass_zdir * 1));

                const  ACC<double> drhs(xdim4_boundary_kernel_mass_zdir, ydim4_boundary_kernel_mass_zdir, drhs_p + (n_x * 1) + (n_y * xdim4_boundary_kernel_mass_zdir * 1) + (n_z * xdim4_boundary_kernel_mass_zdir * ydim4_boundary_kernel_mass_zdir * 1));

                const  ACC<double> urhs(xdim5_boundary_kernel_mass_zdir, ydim5_boundary_kernel_mass_zdir, urhs_p + (n_x * 1) + (n_y * xdim5_boundary_kernel_mass_zdir * 1) + (n_z * xdim5_boundary_kernel_mass_zdir * ydim5_boundary_kernel_mass_zdir * 1));

                const  ACC<double> vrhs(xdim6_boundary_kernel_mass_zdir, ydim6_boundary_kernel_mass_zdir, vrhs_p + (n_x * 1) + (n_y * xdim6_boundary_kernel_mass_zdir * 1) + (n_z * xdim6_boundary_kernel_mass_zdir * ydim6_boundary_kernel_mass_zdir * 1));

                 ACC<double> stryz(xdim7_boundary_kernel_mass_zdir, ydim7_boundary_kernel_mass_zdir, stryz_p + (n_x * 1) + (n_y * xdim7_boundary_kernel_mass_zdir * 1) + (n_z * xdim7_boundary_kernel_mass_zdir * ydim7_boundary_kernel_mass_zdir * 0));

                 ACC<double> bclyz(xdim8_boundary_kernel_mass_zdir, ydim8_boundary_kernel_mass_zdir, bclyz_p + (n_x * 1) + (n_y * xdim8_boundary_kernel_mass_zdir * 1) + (n_z * xdim8_boundary_kernel_mass_zdir * ydim8_boundary_kernel_mass_zdir * 0));

                 ACC<double> t6bz(xdim9_boundary_kernel_mass_zdir, ydim9_boundary_kernel_mass_zdir, t6bz_p + (n_x * 1) + (n_y * xdim9_boundary_kernel_mass_zdir * 1) + (n_z * xdim9_boundary_kernel_mass_zdir * ydim9_boundary_kernel_mass_zdir * 0));

    stryz(0, 0, 0) = yrhs(0, 0, 0);
    bclyz(0, 0, 0) = store3(0, 0, 0);
    t6bz(0, 0, 0) = -store1(0, 0, 0) * urhs(0, 0, 0) / drhs(0, 0, 0) - store2(0, 0, 0) * vrhs(0, 0, 0) / drhs(0, 0, 0);

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[231].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 10);
    ops_set_halo_dirtybit3(&args[7], range);
    ops_set_halo_dirtybit3(&args[8], range);
    ops_set_halo_dirtybit3(&args[9], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[231].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[231].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[231].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[231].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[231].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[231].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[231].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[231].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[231].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
        block->instance->OPS_kernels[231].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg8);
        block->instance->OPS_kernels[231].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg9);
    }
}

#ifdef OPS_LAZY
extern "C"
void boundary_kernel_mass_zdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("boundary_kernel_mass_zdir", args, 10, 231, dim, 0, range, block, boundary_kernel_mass_zdir_execute);
}
#endif

// Auto-generated at 2026-04-28 18:43:19.644758 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void d2fdxz_kernel_eqS_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void d2fdxz_kernel_eqS_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[2];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
#endif

//  ======
//  Timing
//  ======
    double __t1, __t2, __c1, __c2;

    ops_arg arg0 = args[0];
    ops_arg arg1 = args[1];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 2, range, 95)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 95, "d2fdxz_kernel_eqS");
        block->instance->OPS_kernels[95].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "d2fdxz_kernel_eqS");
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
    int xdim0_d2fdxz_kernel_eqS = args[0].dat->size[0];
    int ydim0_d2fdxz_kernel_eqS = args[0].dat->size[1];
    int xdim1_d2fdxz_kernel_eqS = args[1].dat->size[0];
    int ydim1_d2fdxz_kernel_eqS = args[1].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ functn_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ fderiv_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 2);
    ops_halo_exchanges(args, 2, range);
    ops_H_D_exchanges_host(args, 2);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[95].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                const  ACC<double> functn(xdim0_d2fdxz_kernel_eqS, ydim0_d2fdxz_kernel_eqS, functn_p + (n_x * 1) + (n_y * xdim0_d2fdxz_kernel_eqS * 1) + (n_z * xdim0_d2fdxz_kernel_eqS * ydim0_d2fdxz_kernel_eqS * 1));

                 ACC<double> fderiv(xdim1_d2fdxz_kernel_eqS, ydim1_d2fdxz_kernel_eqS, fderiv_p + (n_x * 1) + (n_y * xdim1_d2fdxz_kernel_eqS * 1) + (n_z * xdim1_d2fdxz_kernel_eqS * ydim1_d2fdxz_kernel_eqS * 1));

    double fdiffa;
    double fdiffb;
    double fdiffc;
    double fdiffd;

    fdiffa = acf2xz * (functn(-1, 0, 1) - functn(0, 0, 1) - functn(-1, 0, 0) + functn(0, 0, 0)) + bcf2xz * (functn(1, 0, 1) - functn(0, 0, 1) - functn(1, 0, 0) + functn(0, 0, 0)) + ccf2xz * (functn(2, 0, 1) - functn(0, 0, 1) - functn(2, 0, 0) + functn(0, 0, 0)) + dcf2xz * (functn(3, 0, 1) - functn(0, 0, 1) - functn(3, 0, 0) + functn(0, 0, 0));
    fdiffb = acf2xz * (functn(-1, 0, 2) - functn(0, 0, 2) - functn(-1, 0, 0) + functn(0, 0, 0)) + bcf2xz * (functn(1, 0, 2) - functn(0, 0, 2) - functn(1, 0, 0) + functn(0, 0, 0)) + ccf2xz * (functn(2, 0, 2) - functn(0, 0, 2) - functn(2, 0, 0) + functn(0, 0, 0)) + dcf2xz * (functn(3, 0, 2) - functn(0, 0, 2) - functn(3, 0, 0) + functn(0, 0, 0));
    fdiffc = acf2xz * (functn(-1, 0, 3) - functn(0, 0, 3) - functn(-1, 0, 0) + functn(0, 0, 0)) + bcf2xz * (functn(1, 0, 3) - functn(0, 0, 3) - functn(1, 0, 0) + functn(0, 0, 0)) + ccf2xz * (functn(2, 0, 3) - functn(0, 0, 3) - functn(2, 0, 0) + functn(0, 0, 0)) + dcf2xz * (functn(3, 0, 3) - functn(0, 0, 3) - functn(3, 0, 0) + functn(0, 0, 0));
    fdiffd = acf2xz * (functn(-1, 0, 4) - functn(0, 0, 4) - functn(-1, 0, 0) + functn(0, 0, 0)) + bcf2xz * (functn(1, 0, 4) - functn(0, 0, 4) - functn(1, 0, 0) + functn(0, 0, 0)) + ccf2xz * (functn(2, 0, 4) - functn(0, 0, 4) - functn(2, 0, 0) + functn(0, 0, 0)) + dcf2xz * (functn(3, 0, 4) - functn(0, 0, 4) - functn(3, 0, 0) + functn(0, 0, 0));
    fderiv(0, 0, 0) = acf1xz * fdiffa + bcf1xz * fdiffb + ccf1xz * fdiffc + dcf1xz * fdiffd;

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[95].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 2);
    ops_set_halo_dirtybit3(&args[1], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[95].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[95].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[95].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void d2fdxz_kernel_eqS_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("d2fdxz_kernel_eqS", args, 2, 95, dim, 0, range, block, d2fdxz_kernel_eqS_execute);
}
#endif

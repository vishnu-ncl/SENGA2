// Auto-generated at 2026-04-28 18:43:36.521814 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void boundt_kernel_eqG_zdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void boundt_kernel_eqG_zdir_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[12];
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
    args[10] = desc->args[10];
    args[11] = desc->args[11];
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
    ops_arg arg10 = args[10];
    ops_arg arg11 = args[11];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 12, range, 516)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 516, "boundt_kernel_eqG_zdir");
        block->instance->OPS_kernels[516].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "boundt_kernel_eqG_zdir");
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
    int xdim0_boundt_kernel_eqG_zdir = args[0].dat->size[0];
    int ydim0_boundt_kernel_eqG_zdir = args[0].dat->size[1];
    int xdim1_boundt_kernel_eqG_zdir = args[1].dat->size[0];
    int ydim1_boundt_kernel_eqG_zdir = args[1].dat->size[1];
    int xdim2_boundt_kernel_eqG_zdir = args[2].dat->size[0];
    int ydim2_boundt_kernel_eqG_zdir = args[2].dat->size[1];
    int xdim3_boundt_kernel_eqG_zdir = args[3].dat->size[0];
    int ydim3_boundt_kernel_eqG_zdir = args[3].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ yrhs_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    int * __restrict__ itndex_p = (int *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ strtz_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    double * __restrict__  amasch = (double *)args[4].data;
    double * __restrict__  rgspec = (double *)args[5].data;
    int * __restrict__  ncpoly = (int *)args[6].data;
    int * __restrict__  ncpom1 = (int *)args[7].data;
    int * __restrict__  ncenth = (int *)args[8].data;
    int * __restrict__  ispec = (int *)args[9].data;
    int * __restrict__  icoef1 = (int *)args[10].data;
    int * __restrict__  icoef2 = (int *)args[11].data;

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 12);
    ops_halo_exchanges(args, 12, range);
    ops_H_D_exchanges_host(args, 12);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[516].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                 ACC<double> erhs(xdim0_boundt_kernel_eqG_zdir, ydim0_boundt_kernel_eqG_zdir, erhs_p + (n_x * 1) + (n_y * xdim0_boundt_kernel_eqG_zdir * 1) + (n_z * xdim0_boundt_kernel_eqG_zdir * ydim0_boundt_kernel_eqG_zdir * 1));

                const  ACC<double> yrhs(xdim1_boundt_kernel_eqG_zdir, ydim1_boundt_kernel_eqG_zdir, yrhs_p + (n_x * 1) + (n_y * xdim1_boundt_kernel_eqG_zdir * 1) + (n_z * xdim1_boundt_kernel_eqG_zdir * ydim1_boundt_kernel_eqG_zdir * 1));

                const  ACC<int> itndex(xdim2_boundt_kernel_eqG_zdir, ydim2_boundt_kernel_eqG_zdir, itndex_p + (n_x * 1) + (n_y * xdim2_boundt_kernel_eqG_zdir * 1) + (n_z * xdim2_boundt_kernel_eqG_zdir * ydim2_boundt_kernel_eqG_zdir * 1));

                const  ACC<double> strtz(xdim3_boundt_kernel_eqG_zdir, ydim3_boundt_kernel_eqG_zdir, strtz_p + (n_x * 1) + (n_y * xdim3_boundt_kernel_eqG_zdir * 1) + (n_z * xdim3_boundt_kernel_eqG_zdir * ydim3_boundt_kernel_eqG_zdir * 0));

    double fornow;
    int itint;
    int icp;

    itint = 1 + f2c::mod(itndex(0, 0, 0), icoef1[0]) / icoef2[0];
    fornow = amasch[(ncpoly[(itint-1)+(ispec[0]-1)*ntinmx]-1)+(itint-1)*ncofmx+(ispec[0]-1)*ncofmx*ntinmx];
    for (icp = ncpom1[(itint-1)+(ispec[0]-1)*ntinmx]; icp >= 1; icp -= 1) {
        fornow = fornow * strtz(0, 0, 0) + amasch[(icp-1)+(itint-1)*ncofmx+(ispec[0]-1)*ncofmx*ntinmx];
    }
    fornow = amasch[(ncenth[(itint-1)+(ispec[0]-1)*ntinmx]-1)+(itint-1)*ncofmx+(ispec[0]-1)*ncofmx*ntinmx] + fornow * strtz(0, 0, 0);
    erhs(0, 0, 0) = erhs(0, 0, 0) + (fornow - rgspec[ispec[0]-1] * strtz(0, 0, 0)) * yrhs(0, 0, 0);

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[516].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 12);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[516].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[516].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[516].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[516].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[516].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
    }
}

#ifdef OPS_LAZY
extern "C"
void boundt_kernel_eqG_zdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("boundt_kernel_eqG_zdir", args, 12, 516, dim, 0, range, block, boundt_kernel_eqG_zdir_execute);
}
#endif

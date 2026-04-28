// Auto-generated at 2026-04-28 18:43:25.764168 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void maths_kernel_eqBFG_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void maths_kernel_eqBFG_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[14];
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
    args[12] = desc->args[12];
    args[13] = desc->args[13];
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
    ops_arg arg12 = args[12];
    ops_arg arg13 = args[13];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 14, range, 235)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 235, "maths_kernel_eqBFG");
        block->instance->OPS_kernels[235].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "maths_kernel_eqBFG");
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
    int xdim0_maths_kernel_eqBFG = args[0].dat->size[0];
    int ydim0_maths_kernel_eqBFG = args[0].dat->size[1];
    int xdim1_maths_kernel_eqBFG = args[1].dat->size[0];
    int ydim1_maths_kernel_eqBFG = args[1].dat->size[1];
    int xdim2_maths_kernel_eqBFG = args[2].dat->size[0];
    int ydim2_maths_kernel_eqBFG = args[2].dat->size[1];
    int xdim3_maths_kernel_eqBFG = args[3].dat->size[0];
    int ydim3_maths_kernel_eqBFG = args[3].dat->size[1];
    int xdim4_maths_kernel_eqBFG = args[4].dat->size[0];
    int ydim4_maths_kernel_eqBFG = args[4].dat->size[1];
    int zdim4_maths_kernel_eqBFG = args[4].dat->size[2];
    int xdim5_maths_kernel_eqBFG = args[5].dat->size[0];
    int ydim5_maths_kernel_eqBFG = args[5].dat->size[1];
    int xdim6_maths_kernel_eqBFG = args[6].dat->size[0];
    int ydim6_maths_kernel_eqBFG = args[6].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ difmix_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ store7_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ transp_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ prun_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int multi_d4 = getDatDimFromOpsArg(&args[4]);
    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, multi_d4);
    double * __restrict__ yrhs_md_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ wmomix_p = (double *)(args[5].data) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[6].data) + base6 - 1; // Subtracting 1 to convert to C-style

    double * __restrict__  diffco = (double *)args[7].data;
    double * __restrict__  ovwmol = (double *)args[8].data;
    double * __restrict__  pdifgb = (double *)args[9].data;
    double * __restrict__  dfctol = (double *)args[10].data;
    int * __restrict__  ncodif = (int *)args[11].data;
    int * __restrict__  ncodm1 = (int *)args[12].data;
    int * __restrict__  ispec = (int *)args[13].data;

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 14);
    ops_halo_exchanges(args, 14, range);
    ops_H_D_exchanges_host(args, 14);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[235].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {

                 ACC<double> difmix(xdim0_maths_kernel_eqBFG, ydim0_maths_kernel_eqBFG, difmix_p + (n_x * 1) + (n_y * xdim0_maths_kernel_eqBFG * 1) + (n_z * xdim0_maths_kernel_eqBFG * ydim0_maths_kernel_eqBFG * 1));

                 ACC<double> store7(xdim1_maths_kernel_eqBFG, ydim1_maths_kernel_eqBFG, store7_p + (n_x * 1) + (n_y * xdim1_maths_kernel_eqBFG * 1) + (n_z * xdim1_maths_kernel_eqBFG * ydim1_maths_kernel_eqBFG * 1));

                const  ACC<double> transp(xdim2_maths_kernel_eqBFG, ydim2_maths_kernel_eqBFG, transp_p + (n_x * 1) + (n_y * xdim2_maths_kernel_eqBFG * 1) + (n_z * xdim2_maths_kernel_eqBFG * ydim2_maths_kernel_eqBFG * 1));

                const  ACC<double> prun(xdim3_maths_kernel_eqBFG, ydim3_maths_kernel_eqBFG, prun_p + (n_x * 1) + (n_y * xdim3_maths_kernel_eqBFG * 1) + (n_z * xdim3_maths_kernel_eqBFG * ydim3_maths_kernel_eqBFG * 1));

#ifdef OPS_SOA
                const  ACC<double> yrhs_md(22,xdim4_maths_kernel_eqBFG, ydim4_maths_kernel_eqBFG, zdim4_maths_kernel_eqBFG, yrhs_md_p + (n_x * 1) + (n_y * xdim4_maths_kernel_eqBFG * 1) + (n_z * xdim4_maths_kernel_eqBFG * ydim4_maths_kernel_eqBFG * 1));
#else
                const  ACC<double> yrhs_md(22,xdim4_maths_kernel_eqBFG, ydim4_maths_kernel_eqBFG, zdim4_maths_kernel_eqBFG, yrhs_md_p + 22 * ((n_x * 1) + (n_y * xdim4_maths_kernel_eqBFG * 1) + (n_z * xdim4_maths_kernel_eqBFG * ydim4_maths_kernel_eqBFG * 1)));
#endif

                const  ACC<double> wmomix(xdim5_maths_kernel_eqBFG, ydim5_maths_kernel_eqBFG, wmomix_p + (n_x * 1) + (n_y * xdim5_maths_kernel_eqBFG * 1) + (n_z * xdim5_maths_kernel_eqBFG * ydim5_maths_kernel_eqBFG * 1));

                const  ACC<double> drhs(xdim6_maths_kernel_eqBFG, ydim6_maths_kernel_eqBFG, drhs_p + (n_x * 1) + (n_y * xdim6_maths_kernel_eqBFG * 1) + (n_z * xdim6_maths_kernel_eqBFG * ydim6_maths_kernel_eqBFG * 1));

    double fornow;
    double combo1;
    double combo2;
    double ctrans[(22)];
    int jspec;
    int icp;

    for (jspec = 1; jspec <= nspcmx; ++jspec) {
        fornow = diffco[(ncodif[0]-1)+(jspec-1)*ndcfmx+(ispec[0]-1)*ndcfmx*nspcmx];
        for (icp = ncodm1[0]; icp >= 1; icp -= 1) {
            fornow = fornow * transp(0, 0, 0) + diffco[(icp-1)+(jspec-1)*ndcfmx+(ispec[0]-1)*ndcfmx*nspcmx];
        }
        ctrans[jspec-1] = f2c::exp(fornow) * pdifgb[0] / prun(0, 0, 0);
    }
    combo1 = 0.0;
    combo2 = 0.0;
    for (jspec = 1; jspec <= nspcmx; ++jspec) {
        fornow = yrhs_md(jspec-1, 0, 0, 0) + dfctol[0];
        combo1 = combo1 + fornow;
        combo2 = combo2 + fornow * ovwmol[jspec-1] / ctrans[jspec-1];
    }
    fornow = yrhs_md(ispec[0]-1, 0, 0, 0) + dfctol[0];
    combo1 = combo1 - fornow;
    combo2 = combo2 - fornow * ovwmol[ispec[0]-1] / ctrans[ispec[0]-1];
    combo2 = combo2 * wmomix(0, 0, 0);
    difmix(0, 0, 0) = drhs(0, 0, 0) * combo1 / combo2;
    store7(0, 0, 0) = difmix(0, 0, 0);

            }
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[235].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 14);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[235].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[235].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[235].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[235].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[235].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[235].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[235].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[235].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
    }
}

#ifdef OPS_LAZY
extern "C"
void maths_kernel_eqBFG_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("maths_kernel_eqBFG", args, 14, 235, dim, 0, range, block, maths_kernel_eqBFG_execute);
}
#endif

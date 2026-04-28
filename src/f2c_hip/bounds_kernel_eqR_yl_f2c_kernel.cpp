// Auto-generated at 2026-04-28 18:44:22.184830 by ops-translator

__constant__ int dims_bounds_kernel_eqR_yl[17][3];
static int dims_bounds_kernel_eqR_yl_h[17][3] = {{0}};

//  =============
//  User function
//  =============
__device__ void bounds_kernel_eqR_yl_gpu(ACC<double> &drhs, ACC<double> &urhs, ACC<double> &vrhs, ACC<double> &wrhs, ACC<double> &erhs, const ACC<double> &bcl2yl, const ACC<double> &bcl3yl, const ACC<double> &bcl4yl, const ACC<double> &bcl5yl, const ACC<double> &struyl, const ACC<double> &strvyl, const ACC<double> &strwyl, const ACC<double> &strdyl, const ACC<double> &streyl, const ACC<double> &acouyl, const ACC<double> &ova2yl, const ACC<double> &ovgmyl) {

    drhs(0, 0, 0) = drhs(0, 0, 0) - bcl2yl(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0);
    urhs(0, 0, 0) = urhs(0, 0, 0) - bcl2yl(0, 0, 0) * struyl(0, 0, 0) - bcl3yl(0, 0, 0) * strdyl(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0) * struyl(0, 0, 0);
    vrhs(0, 0, 0) = vrhs(0, 0, 0) - bcl2yl(0, 0, 0) * strvyl(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0) * (strvyl(0, 0, 0) + acouyl(0, 0, 0));
    wrhs(0, 0, 0) = wrhs(0, 0, 0) - bcl2yl(0, 0, 0) * strwyl(0, 0, 0) - bcl4yl(0, 0, 0) * strdyl(0, 0, 0) - bcl5yl(0, 0, 0) * ova2yl(0, 0, 0) * strwyl(0, 0, 0);
    erhs(0, 0, 0) = erhs(0, 0, 0) - bcl2yl(0, 0, 0) * streyl(0, 0, 0) - bcl3yl(0, 0, 0) * strdyl(0, 0, 0) * struyl(0, 0, 0) - bcl4yl(0, 0, 0) * strdyl(0, 0, 0) * strwyl(0, 0, 0) - bcl5yl(0, 0, 0) * (ova2yl(0, 0, 0) * streyl(0, 0, 0) + strvyl(0, 0, 0) / acouyl(0, 0, 0) + ovgmyl(0, 0, 0));

}

//  ============================
//  Cuda kernel wrapper function
//  ============================
__global__ void ops_bounds_kernel_eqR_yl(double* __restrict arg0, int xstride_0, int ystride_0, int zstride_0, 
double* __restrict arg1, int xstride_1, int ystride_1, int zstride_1, 
double* __restrict arg2, int xstride_2, int ystride_2, int zstride_2, 
double* __restrict arg3, int xstride_3, int ystride_3, int zstride_3, 
double* __restrict arg4, int xstride_4, int ystride_4, int zstride_4, 
double* __restrict arg5, int xstride_5, int ystride_5, int zstride_5, 
double* __restrict arg6, int xstride_6, int ystride_6, int zstride_6, 
double* __restrict arg7, int xstride_7, int ystride_7, int zstride_7, 
double* __restrict arg8, int xstride_8, int ystride_8, int zstride_8, 
double* __restrict arg9, int xstride_9, int ystride_9, int zstride_9, 
double* __restrict arg10, int xstride_10, int ystride_10, int zstride_10, 
double* __restrict arg11, int xstride_11, int ystride_11, int zstride_11, 
double* __restrict arg12, int xstride_12, int ystride_12, int zstride_12, 
double* __restrict arg13, int xstride_13, int ystride_13, int zstride_13, 
double* __restrict arg14, int xstride_14, int ystride_14, int zstride_14, 
double* __restrict arg15, int xstride_15, int ystride_15, int zstride_15, 
double* __restrict arg16, int xstride_16, int ystride_16, int zstride_16, 
int size0, int size1, int size2) {

    int idx_z = blockDim.z * blockIdx.z + threadIdx.z;
    int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
    int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

    arg0 += idx_x * xstride_0 + idx_y * ystride_0 * dims_bounds_kernel_eqR_yl[0][0] + idx_z * zstride_0 * dims_bounds_kernel_eqR_yl[0][0] * dims_bounds_kernel_eqR_yl[0][1];
    arg1 += idx_x * xstride_1 + idx_y * ystride_1 * dims_bounds_kernel_eqR_yl[1][0] + idx_z * zstride_1 * dims_bounds_kernel_eqR_yl[1][0] * dims_bounds_kernel_eqR_yl[1][1];
    arg2 += idx_x * xstride_2 + idx_y * ystride_2 * dims_bounds_kernel_eqR_yl[2][0] + idx_z * zstride_2 * dims_bounds_kernel_eqR_yl[2][0] * dims_bounds_kernel_eqR_yl[2][1];
    arg3 += idx_x * xstride_3 + idx_y * ystride_3 * dims_bounds_kernel_eqR_yl[3][0] + idx_z * zstride_3 * dims_bounds_kernel_eqR_yl[3][0] * dims_bounds_kernel_eqR_yl[3][1];
    arg4 += idx_x * xstride_4 + idx_y * ystride_4 * dims_bounds_kernel_eqR_yl[4][0] + idx_z * zstride_4 * dims_bounds_kernel_eqR_yl[4][0] * dims_bounds_kernel_eqR_yl[4][1];
    arg5 += idx_x * xstride_5 + idx_y * ystride_5 * dims_bounds_kernel_eqR_yl[5][0] + idx_z * zstride_5 * dims_bounds_kernel_eqR_yl[5][0] * dims_bounds_kernel_eqR_yl[5][1];
    arg6 += idx_x * xstride_6 + idx_y * ystride_6 * dims_bounds_kernel_eqR_yl[6][0] + idx_z * zstride_6 * dims_bounds_kernel_eqR_yl[6][0] * dims_bounds_kernel_eqR_yl[6][1];
    arg7 += idx_x * xstride_7 + idx_y * ystride_7 * dims_bounds_kernel_eqR_yl[7][0] + idx_z * zstride_7 * dims_bounds_kernel_eqR_yl[7][0] * dims_bounds_kernel_eqR_yl[7][1];
    arg8 += idx_x * xstride_8 + idx_y * ystride_8 * dims_bounds_kernel_eqR_yl[8][0] + idx_z * zstride_8 * dims_bounds_kernel_eqR_yl[8][0] * dims_bounds_kernel_eqR_yl[8][1];
    arg9 += idx_x * xstride_9 + idx_y * ystride_9 * dims_bounds_kernel_eqR_yl[9][0] + idx_z * zstride_9 * dims_bounds_kernel_eqR_yl[9][0] * dims_bounds_kernel_eqR_yl[9][1];
    arg10 += idx_x * xstride_10 + idx_y * ystride_10 * dims_bounds_kernel_eqR_yl[10][0] + idx_z * zstride_10 * dims_bounds_kernel_eqR_yl[10][0] * dims_bounds_kernel_eqR_yl[10][1];
    arg11 += idx_x * xstride_11 + idx_y * ystride_11 * dims_bounds_kernel_eqR_yl[11][0] + idx_z * zstride_11 * dims_bounds_kernel_eqR_yl[11][0] * dims_bounds_kernel_eqR_yl[11][1];
    arg12 += idx_x * xstride_12 + idx_y * ystride_12 * dims_bounds_kernel_eqR_yl[12][0] + idx_z * zstride_12 * dims_bounds_kernel_eqR_yl[12][0] * dims_bounds_kernel_eqR_yl[12][1];
    arg13 += idx_x * xstride_13 + idx_y * ystride_13 * dims_bounds_kernel_eqR_yl[13][0] + idx_z * zstride_13 * dims_bounds_kernel_eqR_yl[13][0] * dims_bounds_kernel_eqR_yl[13][1];
    arg14 += idx_x * xstride_14 + idx_y * ystride_14 * dims_bounds_kernel_eqR_yl[14][0] + idx_z * zstride_14 * dims_bounds_kernel_eqR_yl[14][0] * dims_bounds_kernel_eqR_yl[14][1];
    arg15 += idx_x * xstride_15 + idx_y * ystride_15 * dims_bounds_kernel_eqR_yl[15][0] + idx_z * zstride_15 * dims_bounds_kernel_eqR_yl[15][0] * dims_bounds_kernel_eqR_yl[15][1];
    arg16 += idx_x * xstride_16 + idx_y * ystride_16 * dims_bounds_kernel_eqR_yl[16][0] + idx_z * zstride_16 * dims_bounds_kernel_eqR_yl[16][0] * dims_bounds_kernel_eqR_yl[16][1];

    if(idx_x < size0 && idx_y < size1 && idx_z < size2) {

        ACC<double> argp0(dims_bounds_kernel_eqR_yl[0][0], dims_bounds_kernel_eqR_yl[0][1], arg0);
        ACC<double> argp1(dims_bounds_kernel_eqR_yl[1][0], dims_bounds_kernel_eqR_yl[1][1], arg1);
        ACC<double> argp2(dims_bounds_kernel_eqR_yl[2][0], dims_bounds_kernel_eqR_yl[2][1], arg2);
        ACC<double> argp3(dims_bounds_kernel_eqR_yl[3][0], dims_bounds_kernel_eqR_yl[3][1], arg3);
        ACC<double> argp4(dims_bounds_kernel_eqR_yl[4][0], dims_bounds_kernel_eqR_yl[4][1], arg4);
        const ACC<double> argp5(dims_bounds_kernel_eqR_yl[5][0], dims_bounds_kernel_eqR_yl[5][1], arg5);
        const ACC<double> argp6(dims_bounds_kernel_eqR_yl[6][0], dims_bounds_kernel_eqR_yl[6][1], arg6);
        const ACC<double> argp7(dims_bounds_kernel_eqR_yl[7][0], dims_bounds_kernel_eqR_yl[7][1], arg7);
        const ACC<double> argp8(dims_bounds_kernel_eqR_yl[8][0], dims_bounds_kernel_eqR_yl[8][1], arg8);
        const ACC<double> argp9(dims_bounds_kernel_eqR_yl[9][0], dims_bounds_kernel_eqR_yl[9][1], arg9);
        const ACC<double> argp10(dims_bounds_kernel_eqR_yl[10][0], dims_bounds_kernel_eqR_yl[10][1], arg10);
        const ACC<double> argp11(dims_bounds_kernel_eqR_yl[11][0], dims_bounds_kernel_eqR_yl[11][1], arg11);
        const ACC<double> argp12(dims_bounds_kernel_eqR_yl[12][0], dims_bounds_kernel_eqR_yl[12][1], arg12);
        const ACC<double> argp13(dims_bounds_kernel_eqR_yl[13][0], dims_bounds_kernel_eqR_yl[13][1], arg13);
        const ACC<double> argp14(dims_bounds_kernel_eqR_yl[14][0], dims_bounds_kernel_eqR_yl[14][1], arg14);
        const ACC<double> argp15(dims_bounds_kernel_eqR_yl[15][0], dims_bounds_kernel_eqR_yl[15][1], arg15);
        const ACC<double> argp16(dims_bounds_kernel_eqR_yl[16][0], dims_bounds_kernel_eqR_yl[16][1], arg16);

        bounds_kernel_eqR_yl_gpu(argp0, argp1, argp2, argp3, argp4, argp5, argp6, argp7, argp8, argp9, argp10, argp11, argp12, argp13, argp14, argp15, argp16);

    }// End of cuda index in_range check

}// End of cuda kernel wrapper function

//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bounds_kernel_eqR_yl_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bounds_kernel_eqR_yl_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[17];
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
    args[14] = desc->args[14];
    args[15] = desc->args[15];
    args[16] = desc->args[16];
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
    ops_arg arg14 = args[14];
    ops_arg arg15 = args[15];
    ops_arg arg16 = args[16];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 17, range, 402)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 402, "bounds_kernel_eqR_yl");
        block->instance->OPS_kernels[402].count++;
        ops_timers_core(&__c1, &__t1);
    }

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
    int xdim0 = args[0].dat->size[0];
    int ydim0 = args[0].dat->size[1];
    int zdim0 = args[0].dat->size[2];
    int xdim1 = args[1].dat->size[0];
    int ydim1 = args[1].dat->size[1];
    int zdim1 = args[1].dat->size[2];
    int xdim2 = args[2].dat->size[0];
    int ydim2 = args[2].dat->size[1];
    int zdim2 = args[2].dat->size[2];
    int xdim3 = args[3].dat->size[0];
    int ydim3 = args[3].dat->size[1];
    int zdim3 = args[3].dat->size[2];
    int xdim4 = args[4].dat->size[0];
    int ydim4 = args[4].dat->size[1];
    int zdim4 = args[4].dat->size[2];
    int xdim5 = args[5].dat->size[0];
    int ydim5 = args[5].dat->size[1];
    int zdim5 = args[5].dat->size[2];
    int xdim6 = args[6].dat->size[0];
    int ydim6 = args[6].dat->size[1];
    int zdim6 = args[6].dat->size[2];
    int xdim7 = args[7].dat->size[0];
    int ydim7 = args[7].dat->size[1];
    int zdim7 = args[7].dat->size[2];
    int xdim8 = args[8].dat->size[0];
    int ydim8 = args[8].dat->size[1];
    int zdim8 = args[8].dat->size[2];
    int xdim9 = args[9].dat->size[0];
    int ydim9 = args[9].dat->size[1];
    int zdim9 = args[9].dat->size[2];
    int xdim10 = args[10].dat->size[0];
    int ydim10 = args[10].dat->size[1];
    int zdim10 = args[10].dat->size[2];
    int xdim11 = args[11].dat->size[0];
    int ydim11 = args[11].dat->size[1];
    int zdim11 = args[11].dat->size[2];
    int xdim12 = args[12].dat->size[0];
    int ydim12 = args[12].dat->size[1];
    int zdim12 = args[12].dat->size[2];
    int xdim13 = args[13].dat->size[0];
    int ydim13 = args[13].dat->size[1];
    int zdim13 = args[13].dat->size[2];
    int xdim14 = args[14].dat->size[0];
    int ydim14 = args[14].dat->size[1];
    int zdim14 = args[14].dat->size[2];
    int xdim15 = args[15].dat->size[0];
    int ydim15 = args[15].dat->size[1];
    int zdim15 = args[15].dat->size[2];
    int xdim16 = args[16].dat->size[0];
    int ydim16 = args[16].dat->size[1];
    int zdim16 = args[16].dat->size[2];

    if (xdim0 != dims_bounds_kernel_eqR_yl_h[0][0] || ydim0 != dims_bounds_kernel_eqR_yl_h[0][1] || zdim0 != dims_bounds_kernel_eqR_yl_h[0][2] || xdim1 != dims_bounds_kernel_eqR_yl_h[1][0] || ydim1 != dims_bounds_kernel_eqR_yl_h[1][1] || zdim1 != dims_bounds_kernel_eqR_yl_h[1][2] || xdim2 != dims_bounds_kernel_eqR_yl_h[2][0] || ydim2 != dims_bounds_kernel_eqR_yl_h[2][1] || zdim2 != dims_bounds_kernel_eqR_yl_h[2][2] || xdim3 != dims_bounds_kernel_eqR_yl_h[3][0] || ydim3 != dims_bounds_kernel_eqR_yl_h[3][1] || zdim3 != dims_bounds_kernel_eqR_yl_h[3][2] || xdim4 != dims_bounds_kernel_eqR_yl_h[4][0] || ydim4 != dims_bounds_kernel_eqR_yl_h[4][1] || zdim4 != dims_bounds_kernel_eqR_yl_h[4][2] || xdim5 != dims_bounds_kernel_eqR_yl_h[5][0] || ydim5 != dims_bounds_kernel_eqR_yl_h[5][1] || zdim5 != dims_bounds_kernel_eqR_yl_h[5][2] || xdim6 != dims_bounds_kernel_eqR_yl_h[6][0] || ydim6 != dims_bounds_kernel_eqR_yl_h[6][1] || zdim6 != dims_bounds_kernel_eqR_yl_h[6][2] || xdim7 != dims_bounds_kernel_eqR_yl_h[7][0] || ydim7 != dims_bounds_kernel_eqR_yl_h[7][1] || zdim7 != dims_bounds_kernel_eqR_yl_h[7][2] || xdim8 != dims_bounds_kernel_eqR_yl_h[8][0] || ydim8 != dims_bounds_kernel_eqR_yl_h[8][1] || zdim8 != dims_bounds_kernel_eqR_yl_h[8][2] || xdim9 != dims_bounds_kernel_eqR_yl_h[9][0] || ydim9 != dims_bounds_kernel_eqR_yl_h[9][1] || zdim9 != dims_bounds_kernel_eqR_yl_h[9][2] || xdim10 != dims_bounds_kernel_eqR_yl_h[10][0] || ydim10 != dims_bounds_kernel_eqR_yl_h[10][1] || zdim10 != dims_bounds_kernel_eqR_yl_h[10][2] || xdim11 != dims_bounds_kernel_eqR_yl_h[11][0] || ydim11 != dims_bounds_kernel_eqR_yl_h[11][1] || zdim11 != dims_bounds_kernel_eqR_yl_h[11][2] || xdim12 != dims_bounds_kernel_eqR_yl_h[12][0] || ydim12 != dims_bounds_kernel_eqR_yl_h[12][1] || zdim12 != dims_bounds_kernel_eqR_yl_h[12][2] || xdim13 != dims_bounds_kernel_eqR_yl_h[13][0] || ydim13 != dims_bounds_kernel_eqR_yl_h[13][1] || zdim13 != dims_bounds_kernel_eqR_yl_h[13][2] || xdim14 != dims_bounds_kernel_eqR_yl_h[14][0] || ydim14 != dims_bounds_kernel_eqR_yl_h[14][1] || zdim14 != dims_bounds_kernel_eqR_yl_h[14][2] || xdim15 != dims_bounds_kernel_eqR_yl_h[15][0] || ydim15 != dims_bounds_kernel_eqR_yl_h[15][1] || zdim15 != dims_bounds_kernel_eqR_yl_h[15][2] || xdim16 != dims_bounds_kernel_eqR_yl_h[16][0] || ydim16 != dims_bounds_kernel_eqR_yl_h[16][1] || zdim16 != dims_bounds_kernel_eqR_yl_h[16][2]) {
        dims_bounds_kernel_eqR_yl_h[0][0] = xdim0;
        dims_bounds_kernel_eqR_yl_h[0][1] = ydim0;
        dims_bounds_kernel_eqR_yl_h[0][2] = zdim0;
        dims_bounds_kernel_eqR_yl_h[1][0] = xdim1;
        dims_bounds_kernel_eqR_yl_h[1][1] = ydim1;
        dims_bounds_kernel_eqR_yl_h[1][2] = zdim1;
        dims_bounds_kernel_eqR_yl_h[2][0] = xdim2;
        dims_bounds_kernel_eqR_yl_h[2][1] = ydim2;
        dims_bounds_kernel_eqR_yl_h[2][2] = zdim2;
        dims_bounds_kernel_eqR_yl_h[3][0] = xdim3;
        dims_bounds_kernel_eqR_yl_h[3][1] = ydim3;
        dims_bounds_kernel_eqR_yl_h[3][2] = zdim3;
        dims_bounds_kernel_eqR_yl_h[4][0] = xdim4;
        dims_bounds_kernel_eqR_yl_h[4][1] = ydim4;
        dims_bounds_kernel_eqR_yl_h[4][2] = zdim4;
        dims_bounds_kernel_eqR_yl_h[5][0] = xdim5;
        dims_bounds_kernel_eqR_yl_h[5][1] = ydim5;
        dims_bounds_kernel_eqR_yl_h[5][2] = zdim5;
        dims_bounds_kernel_eqR_yl_h[6][0] = xdim6;
        dims_bounds_kernel_eqR_yl_h[6][1] = ydim6;
        dims_bounds_kernel_eqR_yl_h[6][2] = zdim6;
        dims_bounds_kernel_eqR_yl_h[7][0] = xdim7;
        dims_bounds_kernel_eqR_yl_h[7][1] = ydim7;
        dims_bounds_kernel_eqR_yl_h[7][2] = zdim7;
        dims_bounds_kernel_eqR_yl_h[8][0] = xdim8;
        dims_bounds_kernel_eqR_yl_h[8][1] = ydim8;
        dims_bounds_kernel_eqR_yl_h[8][2] = zdim8;
        dims_bounds_kernel_eqR_yl_h[9][0] = xdim9;
        dims_bounds_kernel_eqR_yl_h[9][1] = ydim9;
        dims_bounds_kernel_eqR_yl_h[9][2] = zdim9;
        dims_bounds_kernel_eqR_yl_h[10][0] = xdim10;
        dims_bounds_kernel_eqR_yl_h[10][1] = ydim10;
        dims_bounds_kernel_eqR_yl_h[10][2] = zdim10;
        dims_bounds_kernel_eqR_yl_h[11][0] = xdim11;
        dims_bounds_kernel_eqR_yl_h[11][1] = ydim11;
        dims_bounds_kernel_eqR_yl_h[11][2] = zdim11;
        dims_bounds_kernel_eqR_yl_h[12][0] = xdim12;
        dims_bounds_kernel_eqR_yl_h[12][1] = ydim12;
        dims_bounds_kernel_eqR_yl_h[12][2] = zdim12;
        dims_bounds_kernel_eqR_yl_h[13][0] = xdim13;
        dims_bounds_kernel_eqR_yl_h[13][1] = ydim13;
        dims_bounds_kernel_eqR_yl_h[13][2] = zdim13;
        dims_bounds_kernel_eqR_yl_h[14][0] = xdim14;
        dims_bounds_kernel_eqR_yl_h[14][1] = ydim14;
        dims_bounds_kernel_eqR_yl_h[14][2] = zdim14;
        dims_bounds_kernel_eqR_yl_h[15][0] = xdim15;
        dims_bounds_kernel_eqR_yl_h[15][1] = ydim15;
        dims_bounds_kernel_eqR_yl_h[15][2] = zdim15;
        dims_bounds_kernel_eqR_yl_h[16][0] = xdim16;
        dims_bounds_kernel_eqR_yl_h[16][1] = ydim16;
        dims_bounds_kernel_eqR_yl_h[16][2] = zdim16;

        hipSafeCall(block->instance->ostream(), hipMemcpyToSymbol( dims_bounds_kernel_eqR_yl, dims_bounds_kernel_eqR_yl_h, sizeof(dims_bounds_kernel_eqR_yl)));
    }

    int x_size = MAX(0,end_indx[0]-start_indx[0]+1);
    int y_size = MAX(0,end_indx[1]-start_indx[1]+1);
    int z_size = MAX(0,end_indx[2]-start_indx[2]+1);

    dim3 grid( (x_size-1)/block->instance->OPS_block_size_x + 1, (y_size-1)/block->instance->OPS_block_size_y + 1, (z_size-1)/block->instance->OPS_block_size_z + 1);

    dim3 tblock(block->instance->OPS_block_size_x,block->instance->OPS_block_size_y,block->instance->OPS_block_size_z);

    char *p_a[17];

//  =======================
//  Set up initial pointers
//  =======================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style
    p_a[0] = (char *)drhs_p;

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ urhs_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style
    p_a[1] = (char *)urhs_p;

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ vrhs_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style
    p_a[2] = (char *)vrhs_p;

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ wrhs_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style
    p_a[3] = (char *)wrhs_p;

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style
    p_a[4] = (char *)erhs_p;

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ bcl2yl_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style
    p_a[5] = (char *)bcl2yl_p;

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ bcl3yl_p = (double *)(args[6].data_d) + base6 - 1; // Subtracting 1 to convert to C-style
    p_a[6] = (char *)bcl3yl_p;

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ bcl4yl_p = (double *)(args[7].data_d) + base7 - 1; // Subtracting 1 to convert to C-style
    p_a[7] = (char *)bcl4yl_p;

    int base8 = getDatBaseFromOpsArg3D(&args[8], start_indx, 1);
    double * __restrict__ bcl5yl_p = (double *)(args[8].data_d) + base8 - 1; // Subtracting 1 to convert to C-style
    p_a[8] = (char *)bcl5yl_p;

    int base9 = getDatBaseFromOpsArg3D(&args[9], start_indx, 1);
    double * __restrict__ struyl_p = (double *)(args[9].data_d) + base9 - 1; // Subtracting 1 to convert to C-style
    p_a[9] = (char *)struyl_p;

    int base10 = getDatBaseFromOpsArg3D(&args[10], start_indx, 1);
    double * __restrict__ strvyl_p = (double *)(args[10].data_d) + base10 - 1; // Subtracting 1 to convert to C-style
    p_a[10] = (char *)strvyl_p;

    int base11 = getDatBaseFromOpsArg3D(&args[11], start_indx, 1);
    double * __restrict__ strwyl_p = (double *)(args[11].data_d) + base11 - 1; // Subtracting 1 to convert to C-style
    p_a[11] = (char *)strwyl_p;

    int base12 = getDatBaseFromOpsArg3D(&args[12], start_indx, 1);
    double * __restrict__ strdyl_p = (double *)(args[12].data_d) + base12 - 1; // Subtracting 1 to convert to C-style
    p_a[12] = (char *)strdyl_p;

    int base13 = getDatBaseFromOpsArg3D(&args[13], start_indx, 1);
    double * __restrict__ streyl_p = (double *)(args[13].data_d) + base13 - 1; // Subtracting 1 to convert to C-style
    p_a[13] = (char *)streyl_p;

    int base14 = getDatBaseFromOpsArg3D(&args[14], start_indx, 1);
    double * __restrict__ acouyl_p = (double *)(args[14].data_d) + base14 - 1; // Subtracting 1 to convert to C-style
    p_a[14] = (char *)acouyl_p;

    int base15 = getDatBaseFromOpsArg3D(&args[15], start_indx, 1);
    double * __restrict__ ova2yl_p = (double *)(args[15].data_d) + base15 - 1; // Subtracting 1 to convert to C-style
    p_a[15] = (char *)ova2yl_p;

    int base16 = getDatBaseFromOpsArg3D(&args[16], start_indx, 1);
    double * __restrict__ ovgmyl_p = (double *)(args[16].data_d) + base16 - 1; // Subtracting 1 to convert to C-style
    p_a[16] = (char *)ovgmyl_p;

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 17);
    ops_halo_exchanges(args, 17, range);
#endif

    if (block->instance->OPS_diags > 1) { 
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[402].mpi_time += __t2 - __t1;
    }

//  ==========================================================
//  ops_dat strides for offset calculation in wrapper function
//  ==========================================================
    int xstride_0, ystride_0, zstride_0;
    xstride_0 = args[0].stencil->stride[0];    ystride_0 = args[0].stencil->stride[1];
    zstride_0 = args[0].stencil->stride[2];
    int xstride_1, ystride_1, zstride_1;
    xstride_1 = args[1].stencil->stride[0];    ystride_1 = args[1].stencil->stride[1];
    zstride_1 = args[1].stencil->stride[2];
    int xstride_2, ystride_2, zstride_2;
    xstride_2 = args[2].stencil->stride[0];    ystride_2 = args[2].stencil->stride[1];
    zstride_2 = args[2].stencil->stride[2];
    int xstride_3, ystride_3, zstride_3;
    xstride_3 = args[3].stencil->stride[0];    ystride_3 = args[3].stencil->stride[1];
    zstride_3 = args[3].stencil->stride[2];
    int xstride_4, ystride_4, zstride_4;
    xstride_4 = args[4].stencil->stride[0];    ystride_4 = args[4].stencil->stride[1];
    zstride_4 = args[4].stencil->stride[2];
    int xstride_5, ystride_5, zstride_5;
    xstride_5 = args[5].stencil->stride[0];    ystride_5 = args[5].stencil->stride[1];
    zstride_5 = args[5].stencil->stride[2];
    int xstride_6, ystride_6, zstride_6;
    xstride_6 = args[6].stencil->stride[0];    ystride_6 = args[6].stencil->stride[1];
    zstride_6 = args[6].stencil->stride[2];
    int xstride_7, ystride_7, zstride_7;
    xstride_7 = args[7].stencil->stride[0];    ystride_7 = args[7].stencil->stride[1];
    zstride_7 = args[7].stencil->stride[2];
    int xstride_8, ystride_8, zstride_8;
    xstride_8 = args[8].stencil->stride[0];    ystride_8 = args[8].stencil->stride[1];
    zstride_8 = args[8].stencil->stride[2];
    int xstride_9, ystride_9, zstride_9;
    xstride_9 = args[9].stencil->stride[0];    ystride_9 = args[9].stencil->stride[1];
    zstride_9 = args[9].stencil->stride[2];
    int xstride_10, ystride_10, zstride_10;
    xstride_10 = args[10].stencil->stride[0];    ystride_10 = args[10].stencil->stride[1];
    zstride_10 = args[10].stencil->stride[2];
    int xstride_11, ystride_11, zstride_11;
    xstride_11 = args[11].stencil->stride[0];    ystride_11 = args[11].stencil->stride[1];
    zstride_11 = args[11].stencil->stride[2];
    int xstride_12, ystride_12, zstride_12;
    xstride_12 = args[12].stencil->stride[0];    ystride_12 = args[12].stencil->stride[1];
    zstride_12 = args[12].stencil->stride[2];
    int xstride_13, ystride_13, zstride_13;
    xstride_13 = args[13].stencil->stride[0];    ystride_13 = args[13].stencil->stride[1];
    zstride_13 = args[13].stencil->stride[2];
    int xstride_14, ystride_14, zstride_14;
    xstride_14 = args[14].stencil->stride[0];    ystride_14 = args[14].stencil->stride[1];
    zstride_14 = args[14].stencil->stride[2];
    int xstride_15, ystride_15, zstride_15;
    xstride_15 = args[15].stencil->stride[0];    ystride_15 = args[15].stencil->stride[1];
    zstride_15 = args[15].stencil->stride[2];
    int xstride_16, ystride_16, zstride_16;
    xstride_16 = args[16].stencil->stride[0];    ystride_16 = args[16].stencil->stride[1];
    zstride_16 = args[16].stencil->stride[2];

//  call kernel wrapper function, passing in pointers to data
    if (x_size > 0 && y_size > 0 && z_size > 0) {

        ops_bounds_kernel_eqR_yl<<<grid, tblock >>> (
                   (double *)p_a[0], xstride_0, ystride_0, zstride_0, 
                   (double *)p_a[1], xstride_1, ystride_1, zstride_1, 
                   (double *)p_a[2], xstride_2, ystride_2, zstride_2, 
                   (double *)p_a[3], xstride_3, ystride_3, zstride_3, 
                   (double *)p_a[4], xstride_4, ystride_4, zstride_4, 
                   (double *)p_a[5], xstride_5, ystride_5, zstride_5, 
                   (double *)p_a[6], xstride_6, ystride_6, zstride_6, 
                   (double *)p_a[7], xstride_7, ystride_7, zstride_7, 
                   (double *)p_a[8], xstride_8, ystride_8, zstride_8, 
                   (double *)p_a[9], xstride_9, ystride_9, zstride_9, 
                   (double *)p_a[10], xstride_10, ystride_10, zstride_10, 
                   (double *)p_a[11], xstride_11, ystride_11, zstride_11, 
                   (double *)p_a[12], xstride_12, ystride_12, zstride_12, 
                   (double *)p_a[13], xstride_13, ystride_13, zstride_13, 
                   (double *)p_a[14], xstride_14, ystride_14, zstride_14, 
                   (double *)p_a[15], xstride_15, ystride_15, zstride_15, 
                   (double *)p_a[16], xstride_16, ystride_16, zstride_16, 
                   x_size, y_size, z_size);

    }

    hipSafeCall(block->instance->ostream(), hipGetLastError());

    if(block->instance->OPS_diags > 1) {
        hipSafeCall(block->instance->ostream(), hipDeviceSynchronize());
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[402].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 17);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
    ops_set_halo_dirtybit3(&args[3], range);
    ops_set_halo_dirtybit3(&args[4], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[402].mpi_time += __t2 - __t1;
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg8);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg9);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg10);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg11);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg12);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg13);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg14);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg15);
        block->instance->OPS_kernels[402].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg16);
    }
}

#ifdef OPS_LAZY
extern "C"
void bounds_kernel_eqR_yl_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bounds_kernel_eqR_yl", args, 17, 402, dim, 1, range, block, bounds_kernel_eqR_yl_execute);
}
#endif

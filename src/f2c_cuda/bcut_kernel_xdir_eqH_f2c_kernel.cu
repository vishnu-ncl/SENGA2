// Auto-generated at 2026-04-28 18:43:56.434310 by ops-translator

__constant__ int dims_bcut_kernel_xdir_eqH[16][3];
static int dims_bcut_kernel_xdir_eqH_h[16][3] = {{0}};

//  =============
//  User function
//  =============
__device__ void bcut_kernel_xdir_eqH_gpu(ACC<double> &strux, ACC<double> &strvx, ACC<double> &strwx, ACC<double> &dudtx, ACC<double> &dvdtx, ACC<double> &dwdtx, const double *slope, const double *widthp, const double *ptly, const double *rxlprm_1, const double *rxlprc_7, const double *fornow, const double *argmnt, const double *deltagy, const double *ygdlen, const int *idx) {

    int jc;
    double hyptan;

    jc = idx[1];
    hyptan = (0.5 * f2c::tanh(slope[0] * (f2c::dble(jc - 1) * deltagy[0] - (ygdlen[0] / 2.0 - widthp[0] * ptly[0]))) + 0.5) * (-0.5 * f2c::tanh(slope[0] * (f2c::dble(jc - 1) * deltagy[0] - (ygdlen[0] / 2.0 + widthp[0] * ptly[0]))) + 0.5);
    strux(0, 0, 0) = rxlprm_1[0] + hyptan * rxlprc_7[0] * f2c::sin(argmnt[0]);
    strvx(0, 0, 0) = 0.0;
    strwx(0, 0, 0) = 0.0;
    dudtx(0, 0, 0) = hyptan * fornow[0] * rxlprc_7[0] * f2c::cos(argmnt[0]);
    dvdtx(0, 0, 0) = 0.0;
    dwdtx(0, 0, 0) = 0.0;

}

//  ============================
//  Cuda kernel wrapper function
//  ============================
__global__ void ops_bcut_kernel_xdir_eqH(double* __restrict arg0, int xstride_0, int ystride_0, int zstride_0, 
double* __restrict arg1, int xstride_1, int ystride_1, int zstride_1, 
double* __restrict arg2, int xstride_2, int ystride_2, int zstride_2, 
double* __restrict arg3, int xstride_3, int ystride_3, int zstride_3, 
double* __restrict arg4, int xstride_4, int ystride_4, int zstride_4, 
double* __restrict arg5, int xstride_5, int ystride_5, int zstride_5, 
const double arg6, 
const double arg7, 
const double arg8, 
const double arg9, 
const double arg10, 
const double arg11, 
const double arg12, 
const double arg13, 
const double arg14, 
int arg_idx0, int arg_idx1, int arg_idx2, 
int size0, int size1, int size2) {

    int idx_z = blockDim.z * blockIdx.z + threadIdx.z;
    int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
    int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

    int arg_idx[3];
    arg_idx[0] = arg_idx0+idx_x;
    arg_idx[1] = arg_idx1+idx_y;
    arg_idx[2] = arg_idx2+idx_z;

    arg0 += idx_x * xstride_0 + idx_y * ystride_0 * dims_bcut_kernel_xdir_eqH[0][0] + idx_z * zstride_0 * dims_bcut_kernel_xdir_eqH[0][0] * dims_bcut_kernel_xdir_eqH[0][1];
    arg1 += idx_x * xstride_1 + idx_y * ystride_1 * dims_bcut_kernel_xdir_eqH[1][0] + idx_z * zstride_1 * dims_bcut_kernel_xdir_eqH[1][0] * dims_bcut_kernel_xdir_eqH[1][1];
    arg2 += idx_x * xstride_2 + idx_y * ystride_2 * dims_bcut_kernel_xdir_eqH[2][0] + idx_z * zstride_2 * dims_bcut_kernel_xdir_eqH[2][0] * dims_bcut_kernel_xdir_eqH[2][1];
    arg3 += idx_x * xstride_3 + idx_y * ystride_3 * dims_bcut_kernel_xdir_eqH[3][0] + idx_z * zstride_3 * dims_bcut_kernel_xdir_eqH[3][0] * dims_bcut_kernel_xdir_eqH[3][1];
    arg4 += idx_x * xstride_4 + idx_y * ystride_4 * dims_bcut_kernel_xdir_eqH[4][0] + idx_z * zstride_4 * dims_bcut_kernel_xdir_eqH[4][0] * dims_bcut_kernel_xdir_eqH[4][1];
    arg5 += idx_x * xstride_5 + idx_y * ystride_5 * dims_bcut_kernel_xdir_eqH[5][0] + idx_z * zstride_5 * dims_bcut_kernel_xdir_eqH[5][0] * dims_bcut_kernel_xdir_eqH[5][1];

    if(idx_x < size0 && idx_y < size1 && idx_z < size2) {

        ACC<double> argp0(dims_bcut_kernel_xdir_eqH[0][0], dims_bcut_kernel_xdir_eqH[0][1], arg0);
        ACC<double> argp1(dims_bcut_kernel_xdir_eqH[1][0], dims_bcut_kernel_xdir_eqH[1][1], arg1);
        ACC<double> argp2(dims_bcut_kernel_xdir_eqH[2][0], dims_bcut_kernel_xdir_eqH[2][1], arg2);
        ACC<double> argp3(dims_bcut_kernel_xdir_eqH[3][0], dims_bcut_kernel_xdir_eqH[3][1], arg3);
        ACC<double> argp4(dims_bcut_kernel_xdir_eqH[4][0], dims_bcut_kernel_xdir_eqH[4][1], arg4);
        ACC<double> argp5(dims_bcut_kernel_xdir_eqH[5][0], dims_bcut_kernel_xdir_eqH[5][1], arg5);

        bcut_kernel_xdir_eqH_gpu(argp0, argp1, argp2, argp3, argp4, argp5, &arg6, &arg7, &arg8, &arg9, &arg10, &arg11, &arg12, &arg13, &arg14, arg_idx);

    }// End of cuda index in_range check

}// End of cuda kernel wrapper function

//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bcut_kernel_xdir_eqH_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bcut_kernel_xdir_eqH_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[16];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 16, range, 325)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 325, "bcut_kernel_xdir_eqH");
        block->instance->OPS_kernels[325].count++;
        ops_timers_core(&__c1, &__t1);
    }

//  =================================================
//  compute locally allocated range for the sub-block
//  =================================================
    int start_indx[3];
    int end_indx[3];
    int arg_idx[3];

//  Range here is in C-style while start_indx and end_indx is Fortran style
#if defined(OPS_MPI) && !defined(OPS_LAZY)
    if ( getRange(block, start_indx, end_indx, range) < 0 ) return;
#else
    for (int n = 0; n < 3; n++) {
        start_indx[n] = range[2*n] + 1;
        end_indx[n]   = range[2*n+1];
    }
#endif

#ifdef OPS_MPI
    getIdx(block, start_indx, arg_idx);
#else
    arg_idx[0] = start_indx[0];
    arg_idx[1] = start_indx[1];
    arg_idx[2] = start_indx[2];
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

    if (xdim0 != dims_bcut_kernel_xdir_eqH_h[0][0] || ydim0 != dims_bcut_kernel_xdir_eqH_h[0][1] || zdim0 != dims_bcut_kernel_xdir_eqH_h[0][2] || xdim1 != dims_bcut_kernel_xdir_eqH_h[1][0] || ydim1 != dims_bcut_kernel_xdir_eqH_h[1][1] || zdim1 != dims_bcut_kernel_xdir_eqH_h[1][2] || xdim2 != dims_bcut_kernel_xdir_eqH_h[2][0] || ydim2 != dims_bcut_kernel_xdir_eqH_h[2][1] || zdim2 != dims_bcut_kernel_xdir_eqH_h[2][2] || xdim3 != dims_bcut_kernel_xdir_eqH_h[3][0] || ydim3 != dims_bcut_kernel_xdir_eqH_h[3][1] || zdim3 != dims_bcut_kernel_xdir_eqH_h[3][2] || xdim4 != dims_bcut_kernel_xdir_eqH_h[4][0] || ydim4 != dims_bcut_kernel_xdir_eqH_h[4][1] || zdim4 != dims_bcut_kernel_xdir_eqH_h[4][2] || xdim5 != dims_bcut_kernel_xdir_eqH_h[5][0] || ydim5 != dims_bcut_kernel_xdir_eqH_h[5][1] || zdim5 != dims_bcut_kernel_xdir_eqH_h[5][2]) {
        dims_bcut_kernel_xdir_eqH_h[0][0] = xdim0;
        dims_bcut_kernel_xdir_eqH_h[0][1] = ydim0;
        dims_bcut_kernel_xdir_eqH_h[0][2] = zdim0;
        dims_bcut_kernel_xdir_eqH_h[1][0] = xdim1;
        dims_bcut_kernel_xdir_eqH_h[1][1] = ydim1;
        dims_bcut_kernel_xdir_eqH_h[1][2] = zdim1;
        dims_bcut_kernel_xdir_eqH_h[2][0] = xdim2;
        dims_bcut_kernel_xdir_eqH_h[2][1] = ydim2;
        dims_bcut_kernel_xdir_eqH_h[2][2] = zdim2;
        dims_bcut_kernel_xdir_eqH_h[3][0] = xdim3;
        dims_bcut_kernel_xdir_eqH_h[3][1] = ydim3;
        dims_bcut_kernel_xdir_eqH_h[3][2] = zdim3;
        dims_bcut_kernel_xdir_eqH_h[4][0] = xdim4;
        dims_bcut_kernel_xdir_eqH_h[4][1] = ydim4;
        dims_bcut_kernel_xdir_eqH_h[4][2] = zdim4;
        dims_bcut_kernel_xdir_eqH_h[5][0] = xdim5;
        dims_bcut_kernel_xdir_eqH_h[5][1] = ydim5;
        dims_bcut_kernel_xdir_eqH_h[5][2] = zdim5;

        cutilSafeCall(block->instance->ostream(), cudaMemcpyToSymbol( dims_bcut_kernel_xdir_eqH, dims_bcut_kernel_xdir_eqH_h, sizeof(dims_bcut_kernel_xdir_eqH)));
    }

    int x_size = MAX(0,end_indx[0]-start_indx[0]+1);
    int y_size = MAX(0,end_indx[1]-start_indx[1]+1);
    int z_size = MAX(0,end_indx[2]-start_indx[2]+1);

    dim3 grid( (x_size-1)/block->instance->OPS_block_size_x + 1, (y_size-1)/block->instance->OPS_block_size_y + 1, (z_size-1)/block->instance->OPS_block_size_z + 1);

    dim3 tblock(block->instance->OPS_block_size_x,block->instance->OPS_block_size_y,block->instance->OPS_block_size_z);

    char *p_a[16];

//  =======================
//  Set up initial pointers
//  =======================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ strux_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style
    p_a[0] = (char *)strux_p;

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ strvx_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style
    p_a[1] = (char *)strvx_p;

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ strwx_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style
    p_a[2] = (char *)strwx_p;

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ dudtx_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style
    p_a[3] = (char *)dudtx_p;

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ dvdtx_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style
    p_a[4] = (char *)dvdtx_p;

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ dwdtx_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style
    p_a[5] = (char *)dwdtx_p;

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 16);
    ops_halo_exchanges(args, 16, range);
#endif

    if (block->instance->OPS_diags > 1) { 
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[325].mpi_time += __t2 - __t1;
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

//  call kernel wrapper function, passing in pointers to data
    if (x_size > 0 && y_size > 0 && z_size > 0) {

        ops_bcut_kernel_xdir_eqH<<<grid, tblock >>> (
                   (double *)p_a[0], xstride_0, ystride_0, zstride_0, 
                   (double *)p_a[1], xstride_1, ystride_1, zstride_1, 
                   (double *)p_a[2], xstride_2, ystride_2, zstride_2, 
                   (double *)p_a[3], xstride_3, ystride_3, zstride_3, 
                   (double *)p_a[4], xstride_4, ystride_4, zstride_4, 
                   (double *)p_a[5], xstride_5, ystride_5, zstride_5, 
                   *(double *)arg6.data, 
                   *(double *)arg7.data, 
                   *(double *)arg8.data, 
                   *(double *)arg9.data, 
                   *(double *)arg10.data, 
                   *(double *)arg11.data, 
                   *(double *)arg12.data, 
                   *(double *)arg13.data, 
                   *(double *)arg14.data, 
                   arg_idx[0], arg_idx[1], arg_idx[2], 
                   x_size, y_size, z_size);

    }

    cutilSafeCall(block->instance->ostream(), cudaGetLastError());

    if(block->instance->OPS_diags > 1) {
        cutilSafeCall(block->instance->ostream(), cudaDeviceSynchronize());
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[325].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 16);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
    ops_set_halo_dirtybit3(&args[3], range);
    ops_set_halo_dirtybit3(&args[4], range);
    ops_set_halo_dirtybit3(&args[5], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[325].mpi_time += __t2 - __t1;
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[325].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
    }
}

#ifdef OPS_LAZY
extern "C"
void bcut_kernel_xdir_eqH_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bcut_kernel_xdir_eqH", args, 16, 325, dim, 1, range, block, bcut_kernel_xdir_eqH_execute);
}
#endif

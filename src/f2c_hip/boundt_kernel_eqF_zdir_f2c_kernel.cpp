// Auto-generated at 2026-04-28 18:44:27.177270 by ops-translator

__constant__ int dims_boundt_kernel_eqF_zdir[14][3];
static int dims_boundt_kernel_eqF_zdir_h[14][3] = {{0}};

//  =============
//  User function
//  =============
__device__ void boundt_kernel_eqF_zdir_gpu(ACC<double> &erhs, ACC<double> &yrhs, const ACC<int> &itndex, const ACC<double> &drhs, const ACC<double> &strtz, const ACC<double> &stryz, const double *amasch, const double *rgspec, const int *ncpoly, const int *ncpom1, const int *ncenth, const int *ispec, const int *icoef1, const int *icoef2) {

    double fornow;
    int itint;
    int icp;

    itint = 1 + f2c::mod(itndex(0, 0, 0), icoef1[0]) / icoef2[0];
    fornow = amasch[(ncpoly[(itint-1)+(ispec[0]-1)*ntinmx]-1)+(itint-1)*ncofmx+(ispec[0]-1)*ncofmx*ntinmx];
    for (icp = ncpom1[(itint-1)+(ispec[0]-1)*ntinmx]; icp >= 1; icp -= 1) {
        fornow = fornow * strtz(0, 0, 0) + amasch[(icp-1)+(itint-1)*ncofmx+(ispec[0]-1)*ncofmx*ntinmx];
    }
    fornow = amasch[(ncenth[(itint-1)+(ispec[0]-1)*ntinmx]-1)+(itint-1)*ncofmx+(ispec[0]-1)*ncofmx*ntinmx] + fornow * strtz(0, 0, 0);
    yrhs(0, 0, 0) = drhs(0, 0, 0) * stryz(0, 0, 0);
    erhs(0, 0, 0) = erhs(0, 0, 0) + (fornow - rgspec[ispec[0]-1] * strtz(0, 0, 0)) * yrhs(0, 0, 0);

}

//  ============================
//  Cuda kernel wrapper function
//  ============================
__global__ void ops_boundt_kernel_eqF_zdir(double* __restrict arg0, int xstride_0, int ystride_0, int zstride_0, 
double* __restrict arg1, int xstride_1, int ystride_1, int zstride_1, 
int* __restrict arg2, int xstride_2, int ystride_2, int zstride_2, 
double* __restrict arg3, int xstride_3, int ystride_3, int zstride_3, 
double* __restrict arg4, int xstride_4, int ystride_4, int zstride_4, 
double* __restrict arg5, int xstride_5, int ystride_5, int zstride_5, 
const double* __restrict arg6, 
const double* __restrict arg7, 
const int* __restrict arg8, 
const int* __restrict arg9, 
const int* __restrict arg10, 
const int arg11, 
const int arg12, 
const int arg13, 
int size0, int size1, int size2) {

    int idx_z = blockDim.z * blockIdx.z + threadIdx.z;
    int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
    int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

    arg0 += idx_x * xstride_0 + idx_y * ystride_0 * dims_boundt_kernel_eqF_zdir[0][0] + idx_z * zstride_0 * dims_boundt_kernel_eqF_zdir[0][0] * dims_boundt_kernel_eqF_zdir[0][1];
    arg1 += idx_x * xstride_1 + idx_y * ystride_1 * dims_boundt_kernel_eqF_zdir[1][0] + idx_z * zstride_1 * dims_boundt_kernel_eqF_zdir[1][0] * dims_boundt_kernel_eqF_zdir[1][1];
    arg2 += idx_x * xstride_2 + idx_y * ystride_2 * dims_boundt_kernel_eqF_zdir[2][0] + idx_z * zstride_2 * dims_boundt_kernel_eqF_zdir[2][0] * dims_boundt_kernel_eqF_zdir[2][1];
    arg3 += idx_x * xstride_3 + idx_y * ystride_3 * dims_boundt_kernel_eqF_zdir[3][0] + idx_z * zstride_3 * dims_boundt_kernel_eqF_zdir[3][0] * dims_boundt_kernel_eqF_zdir[3][1];
    arg4 += idx_x * xstride_4 + idx_y * ystride_4 * dims_boundt_kernel_eqF_zdir[4][0] + idx_z * zstride_4 * dims_boundt_kernel_eqF_zdir[4][0] * dims_boundt_kernel_eqF_zdir[4][1];
    arg5 += idx_x * xstride_5 + idx_y * ystride_5 * dims_boundt_kernel_eqF_zdir[5][0] + idx_z * zstride_5 * dims_boundt_kernel_eqF_zdir[5][0] * dims_boundt_kernel_eqF_zdir[5][1];

    if(idx_x < size0 && idx_y < size1 && idx_z < size2) {

        ACC<double> argp0(dims_boundt_kernel_eqF_zdir[0][0], dims_boundt_kernel_eqF_zdir[0][1], arg0);
        ACC<double> argp1(dims_boundt_kernel_eqF_zdir[1][0], dims_boundt_kernel_eqF_zdir[1][1], arg1);
        const ACC<int> argp2(dims_boundt_kernel_eqF_zdir[2][0], dims_boundt_kernel_eqF_zdir[2][1], arg2);
        const ACC<double> argp3(dims_boundt_kernel_eqF_zdir[3][0], dims_boundt_kernel_eqF_zdir[3][1], arg3);
        const ACC<double> argp4(dims_boundt_kernel_eqF_zdir[4][0], dims_boundt_kernel_eqF_zdir[4][1], arg4);
        const ACC<double> argp5(dims_boundt_kernel_eqF_zdir[5][0], dims_boundt_kernel_eqF_zdir[5][1], arg5);

        boundt_kernel_eqF_zdir_gpu(argp0, argp1, argp2, argp3, argp4, argp5, arg6, arg7, arg8, arg9, arg10, &arg11, &arg12, &arg13);

    }// End of cuda index in_range check

}// End of cuda kernel wrapper function

//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void boundt_kernel_eqF_zdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void boundt_kernel_eqF_zdir_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 14, range, 531)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 531, "boundt_kernel_eqF_zdir");
        block->instance->OPS_kernels[531].count++;
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

    if (xdim0 != dims_boundt_kernel_eqF_zdir_h[0][0] || ydim0 != dims_boundt_kernel_eqF_zdir_h[0][1] || zdim0 != dims_boundt_kernel_eqF_zdir_h[0][2] || xdim1 != dims_boundt_kernel_eqF_zdir_h[1][0] || ydim1 != dims_boundt_kernel_eqF_zdir_h[1][1] || zdim1 != dims_boundt_kernel_eqF_zdir_h[1][2] || xdim2 != dims_boundt_kernel_eqF_zdir_h[2][0] || ydim2 != dims_boundt_kernel_eqF_zdir_h[2][1] || zdim2 != dims_boundt_kernel_eqF_zdir_h[2][2] || xdim3 != dims_boundt_kernel_eqF_zdir_h[3][0] || ydim3 != dims_boundt_kernel_eqF_zdir_h[3][1] || zdim3 != dims_boundt_kernel_eqF_zdir_h[3][2] || xdim4 != dims_boundt_kernel_eqF_zdir_h[4][0] || ydim4 != dims_boundt_kernel_eqF_zdir_h[4][1] || zdim4 != dims_boundt_kernel_eqF_zdir_h[4][2] || xdim5 != dims_boundt_kernel_eqF_zdir_h[5][0] || ydim5 != dims_boundt_kernel_eqF_zdir_h[5][1] || zdim5 != dims_boundt_kernel_eqF_zdir_h[5][2]) {
        dims_boundt_kernel_eqF_zdir_h[0][0] = xdim0;
        dims_boundt_kernel_eqF_zdir_h[0][1] = ydim0;
        dims_boundt_kernel_eqF_zdir_h[0][2] = zdim0;
        dims_boundt_kernel_eqF_zdir_h[1][0] = xdim1;
        dims_boundt_kernel_eqF_zdir_h[1][1] = ydim1;
        dims_boundt_kernel_eqF_zdir_h[1][2] = zdim1;
        dims_boundt_kernel_eqF_zdir_h[2][0] = xdim2;
        dims_boundt_kernel_eqF_zdir_h[2][1] = ydim2;
        dims_boundt_kernel_eqF_zdir_h[2][2] = zdim2;
        dims_boundt_kernel_eqF_zdir_h[3][0] = xdim3;
        dims_boundt_kernel_eqF_zdir_h[3][1] = ydim3;
        dims_boundt_kernel_eqF_zdir_h[3][2] = zdim3;
        dims_boundt_kernel_eqF_zdir_h[4][0] = xdim4;
        dims_boundt_kernel_eqF_zdir_h[4][1] = ydim4;
        dims_boundt_kernel_eqF_zdir_h[4][2] = zdim4;
        dims_boundt_kernel_eqF_zdir_h[5][0] = xdim5;
        dims_boundt_kernel_eqF_zdir_h[5][1] = ydim5;
        dims_boundt_kernel_eqF_zdir_h[5][2] = zdim5;

        hipSafeCall(block->instance->ostream(), hipMemcpyToSymbol( dims_boundt_kernel_eqF_zdir, dims_boundt_kernel_eqF_zdir_h, sizeof(dims_boundt_kernel_eqF_zdir)));
    }

    double *arg6h = (double*)arg6.data;
    double *arg7h = (double*)arg7.data;
    int *arg8h = (int*)arg8.data;
    int *arg9h = (int*)arg9.data;
    int *arg10h = (int*)arg10.data;

    int x_size = MAX(0,end_indx[0]-start_indx[0]+1);
    int y_size = MAX(0,end_indx[1]-start_indx[1]+1);
    int z_size = MAX(0,end_indx[2]-start_indx[2]+1);

    dim3 grid( (x_size-1)/block->instance->OPS_block_size_x + 1, (y_size-1)/block->instance->OPS_block_size_y + 1, (z_size-1)/block->instance->OPS_block_size_z + 1);

    dim3 tblock(block->instance->OPS_block_size_x,block->instance->OPS_block_size_y,block->instance->OPS_block_size_z);

    int consts_bytes = 0;

    consts_bytes += ROUND_UP(arg6.dim*sizeof(double));
    consts_bytes += ROUND_UP(arg7.dim*sizeof(double));
    consts_bytes += ROUND_UP(arg8.dim*sizeof(int));
    consts_bytes += ROUND_UP(arg9.dim*sizeof(int));
    consts_bytes += ROUND_UP(arg10.dim*sizeof(int));

    reallocConstArrays(block->instance, consts_bytes);
    consts_bytes = 0;

    arg6.data = block->instance->OPS_consts_h + consts_bytes;
    arg6.data_d = block->instance->OPS_consts_d + consts_bytes;
    for (int d = 0; d < arg6.dim; d++)    ((double *)arg6.data)[d] = arg6h[d];
    consts_bytes += ROUND_UP(arg6.dim*sizeof(double));
    arg7.data = block->instance->OPS_consts_h + consts_bytes;
    arg7.data_d = block->instance->OPS_consts_d + consts_bytes;
    for (int d = 0; d < arg7.dim; d++)    ((double *)arg7.data)[d] = arg7h[d];
    consts_bytes += ROUND_UP(arg7.dim*sizeof(double));
    arg8.data = block->instance->OPS_consts_h + consts_bytes;
    arg8.data_d = block->instance->OPS_consts_d + consts_bytes;
    for (int d = 0; d < arg8.dim; d++)    ((int *)arg8.data)[d] = arg8h[d];
    consts_bytes += ROUND_UP(arg8.dim*sizeof(int));
    arg9.data = block->instance->OPS_consts_h + consts_bytes;
    arg9.data_d = block->instance->OPS_consts_d + consts_bytes;
    for (int d = 0; d < arg9.dim; d++)    ((int *)arg9.data)[d] = arg9h[d];
    consts_bytes += ROUND_UP(arg9.dim*sizeof(int));
    arg10.data = block->instance->OPS_consts_h + consts_bytes;
    arg10.data_d = block->instance->OPS_consts_d + consts_bytes;
    for (int d = 0; d < arg10.dim; d++)    ((int *)arg10.data)[d] = arg10h[d];
    consts_bytes += ROUND_UP(arg10.dim*sizeof(int));

    mvConstArraysToDevice(block->instance, consts_bytes);

    char *p_a[14];

//  =======================
//  Set up initial pointers
//  =======================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style
    p_a[0] = (char *)erhs_p;

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ yrhs_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style
    p_a[1] = (char *)yrhs_p;

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    int * __restrict__ itndex_p = (int *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style
    p_a[2] = (char *)itndex_p;

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style
    p_a[3] = (char *)drhs_p;

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ strtz_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style
    p_a[4] = (char *)strtz_p;

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ stryz_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style
    p_a[5] = (char *)stryz_p;

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 14);
    ops_halo_exchanges(args, 14, range);
#endif

    if (block->instance->OPS_diags > 1) { 
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[531].mpi_time += __t2 - __t1;
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

        ops_boundt_kernel_eqF_zdir<<<grid, tblock >>> (
                   (double *)p_a[0], xstride_0, ystride_0, zstride_0, 
                   (double *)p_a[1], xstride_1, ystride_1, zstride_1, 
                   (int *)p_a[2], xstride_2, ystride_2, zstride_2, 
                   (double *)p_a[3], xstride_3, ystride_3, zstride_3, 
                   (double *)p_a[4], xstride_4, ystride_4, zstride_4, 
                   (double *)p_a[5], xstride_5, ystride_5, zstride_5, 
                   (double *)arg6.data_d, 
                   (double *)arg7.data_d, 
                   (int *)arg8.data_d, 
                   (int *)arg9.data_d, 
                   (int *)arg10.data_d, 
                   *(int *)arg11.data, 
                   *(int *)arg12.data, 
                   *(int *)arg13.data, 
                   x_size, y_size, z_size);

    }

    hipSafeCall(block->instance->ostream(), hipGetLastError());

    if(block->instance->OPS_diags > 1) {
        hipSafeCall(block->instance->ostream(), hipDeviceSynchronize());
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[531].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 14);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[531].mpi_time += __t2 - __t1;
        block->instance->OPS_kernels[531].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[531].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[531].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[531].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[531].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[531].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
    }
}

#ifdef OPS_LAZY
extern "C"
void boundt_kernel_eqF_zdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("boundt_kernel_eqF_zdir", args, 14, 531, dim, 1, range, block, boundt_kernel_eqF_zdir_execute);
}
#endif

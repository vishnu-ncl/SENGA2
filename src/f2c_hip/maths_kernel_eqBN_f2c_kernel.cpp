// Auto-generated at 2026-04-28 18:44:27.524718 by ops-translator

__constant__ int dims_maths_kernel_eqBN[15][3];
static int dims_maths_kernel_eqBN_h[15][3] = {{0}};

//  =============
//  User function
//  =============
__device__ void maths_kernel_eqBN_gpu(ACC<double> &store2, const ACC<double> &trun, const ACC<int> &itndex, const double *amolgb, const int *ncpoly, const int *ncpom1, const int *ncenth, const int *ncenpy, const double *diffmu, const int *isspec, const int *ispec, const int *istep, const int *ipower, const int *icoef1, const int *icoef2) {

    int itint;
    int icp;
    double gibbsp;

    itint = 1 + f2c::mod(itndex(0, 0, 0), icoef1[0]) / icoef2[0];
    gibbsp = amolgb[(ncpoly[(itint-1)+(ispec[0]-1)*ntinmx]-1)+(itint-1)*ncofmx+(ispec[0]-1)*ncofmx*ntinmx];
    for (icp = ncpom1[(itint-1)+(ispec[0]-1)*ntinmx]; icp >= 1; icp -= 1) {
        gibbsp = amolgb[(icp-1)+(itint-1)*ncofmx+(ispec[0]-1)*ncofmx*ntinmx] + gibbsp * trun(0, 0, 0);
    }
    gibbsp = amolgb[(ncenth[(itint-1)+(ispec[0]-1)*ntinmx]-1)+(itint-1)*ncofmx+(ispec[0]-1)*ncofmx*ntinmx] / trun(0, 0, 0) - amolgb[(ncenpy[(itint-1)+(ispec[0]-1)*ntinmx]-1)+(itint-1)*ncofmx+(ispec[0]-1)*ncofmx*ntinmx] * f2c::log(trun(0, 0, 0)) - gibbsp;
    store2(0, 0, 0) = store2(0, 0, 0) + diffmu[(isspec[0]-1)+(istep[0]-1)*nssmax] * gibbsp;

}

//  ============================
//  Cuda kernel wrapper function
//  ============================
__global__ void ops_maths_kernel_eqBN(double* __restrict arg0, int xstride_0, int ystride_0, int zstride_0, 
double* __restrict arg1, int xstride_1, int ystride_1, int zstride_1, 
int* __restrict arg2, int xstride_2, int ystride_2, int zstride_2, 
const double* __restrict arg3, 
const int* __restrict arg4, 
const int* __restrict arg5, 
const int* __restrict arg6, 
const int* __restrict arg7, 
const double* __restrict arg8, 
const int arg9, 
const int arg10, 
const int arg11, 
const int arg12, 
const int arg13, 
const int arg14, 
int size0, int size1, int size2) {

    int idx_z = blockDim.z * blockIdx.z + threadIdx.z;
    int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
    int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

    arg0 += idx_x * xstride_0 + idx_y * ystride_0 * dims_maths_kernel_eqBN[0][0] + idx_z * zstride_0 * dims_maths_kernel_eqBN[0][0] * dims_maths_kernel_eqBN[0][1];
    arg1 += idx_x * xstride_1 + idx_y * ystride_1 * dims_maths_kernel_eqBN[1][0] + idx_z * zstride_1 * dims_maths_kernel_eqBN[1][0] * dims_maths_kernel_eqBN[1][1];
    arg2 += idx_x * xstride_2 + idx_y * ystride_2 * dims_maths_kernel_eqBN[2][0] + idx_z * zstride_2 * dims_maths_kernel_eqBN[2][0] * dims_maths_kernel_eqBN[2][1];

    if(idx_x < size0 && idx_y < size1 && idx_z < size2) {

        ACC<double> argp0(dims_maths_kernel_eqBN[0][0], dims_maths_kernel_eqBN[0][1], arg0);
        const ACC<double> argp1(dims_maths_kernel_eqBN[1][0], dims_maths_kernel_eqBN[1][1], arg1);
        const ACC<int> argp2(dims_maths_kernel_eqBN[2][0], dims_maths_kernel_eqBN[2][1], arg2);

        maths_kernel_eqBN_gpu(argp0, argp1, argp2, arg3, arg4, arg5, arg6, arg7, arg8, &arg9, &arg10, &arg11, &arg12, &arg13, &arg14);

    }// End of cuda index in_range check

}// End of cuda kernel wrapper function

//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void maths_kernel_eqBN_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void maths_kernel_eqBN_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[15];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 15, range, 544)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 544, "maths_kernel_eqBN");
        block->instance->OPS_kernels[544].count++;
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

    if (xdim0 != dims_maths_kernel_eqBN_h[0][0] || ydim0 != dims_maths_kernel_eqBN_h[0][1] || zdim0 != dims_maths_kernel_eqBN_h[0][2] || xdim1 != dims_maths_kernel_eqBN_h[1][0] || ydim1 != dims_maths_kernel_eqBN_h[1][1] || zdim1 != dims_maths_kernel_eqBN_h[1][2] || xdim2 != dims_maths_kernel_eqBN_h[2][0] || ydim2 != dims_maths_kernel_eqBN_h[2][1] || zdim2 != dims_maths_kernel_eqBN_h[2][2]) {
        dims_maths_kernel_eqBN_h[0][0] = xdim0;
        dims_maths_kernel_eqBN_h[0][1] = ydim0;
        dims_maths_kernel_eqBN_h[0][2] = zdim0;
        dims_maths_kernel_eqBN_h[1][0] = xdim1;
        dims_maths_kernel_eqBN_h[1][1] = ydim1;
        dims_maths_kernel_eqBN_h[1][2] = zdim1;
        dims_maths_kernel_eqBN_h[2][0] = xdim2;
        dims_maths_kernel_eqBN_h[2][1] = ydim2;
        dims_maths_kernel_eqBN_h[2][2] = zdim2;

        hipSafeCall(block->instance->ostream(), hipMemcpyToSymbol( dims_maths_kernel_eqBN, dims_maths_kernel_eqBN_h, sizeof(dims_maths_kernel_eqBN)));
    }

    double *arg3h = (double*)arg3.data;
    int *arg4h = (int*)arg4.data;
    int *arg5h = (int*)arg5.data;
    int *arg6h = (int*)arg6.data;
    int *arg7h = (int*)arg7.data;
    double *arg8h = (double*)arg8.data;

    int x_size = MAX(0,end_indx[0]-start_indx[0]+1);
    int y_size = MAX(0,end_indx[1]-start_indx[1]+1);
    int z_size = MAX(0,end_indx[2]-start_indx[2]+1);

    dim3 grid( (x_size-1)/block->instance->OPS_block_size_x + 1, (y_size-1)/block->instance->OPS_block_size_y + 1, (z_size-1)/block->instance->OPS_block_size_z + 1);

    dim3 tblock(block->instance->OPS_block_size_x,block->instance->OPS_block_size_y,block->instance->OPS_block_size_z);

    int consts_bytes = 0;

    consts_bytes += ROUND_UP(arg3.dim*sizeof(double));
    consts_bytes += ROUND_UP(arg4.dim*sizeof(int));
    consts_bytes += ROUND_UP(arg5.dim*sizeof(int));
    consts_bytes += ROUND_UP(arg6.dim*sizeof(int));
    consts_bytes += ROUND_UP(arg7.dim*sizeof(int));
    consts_bytes += ROUND_UP(arg8.dim*sizeof(double));

    reallocConstArrays(block->instance, consts_bytes);
    consts_bytes = 0;

    arg3.data = block->instance->OPS_consts_h + consts_bytes;
    arg3.data_d = block->instance->OPS_consts_d + consts_bytes;
    for (int d = 0; d < arg3.dim; d++)    ((double *)arg3.data)[d] = arg3h[d];
    consts_bytes += ROUND_UP(arg3.dim*sizeof(double));
    arg4.data = block->instance->OPS_consts_h + consts_bytes;
    arg4.data_d = block->instance->OPS_consts_d + consts_bytes;
    for (int d = 0; d < arg4.dim; d++)    ((int *)arg4.data)[d] = arg4h[d];
    consts_bytes += ROUND_UP(arg4.dim*sizeof(int));
    arg5.data = block->instance->OPS_consts_h + consts_bytes;
    arg5.data_d = block->instance->OPS_consts_d + consts_bytes;
    for (int d = 0; d < arg5.dim; d++)    ((int *)arg5.data)[d] = arg5h[d];
    consts_bytes += ROUND_UP(arg5.dim*sizeof(int));
    arg6.data = block->instance->OPS_consts_h + consts_bytes;
    arg6.data_d = block->instance->OPS_consts_d + consts_bytes;
    for (int d = 0; d < arg6.dim; d++)    ((int *)arg6.data)[d] = arg6h[d];
    consts_bytes += ROUND_UP(arg6.dim*sizeof(int));
    arg7.data = block->instance->OPS_consts_h + consts_bytes;
    arg7.data_d = block->instance->OPS_consts_d + consts_bytes;
    for (int d = 0; d < arg7.dim; d++)    ((int *)arg7.data)[d] = arg7h[d];
    consts_bytes += ROUND_UP(arg7.dim*sizeof(int));
    arg8.data = block->instance->OPS_consts_h + consts_bytes;
    arg8.data_d = block->instance->OPS_consts_d + consts_bytes;
    for (int d = 0; d < arg8.dim; d++)    ((double *)arg8.data)[d] = arg8h[d];
    consts_bytes += ROUND_UP(arg8.dim*sizeof(double));

    mvConstArraysToDevice(block->instance, consts_bytes);

    char *p_a[15];

//  =======================
//  Set up initial pointers
//  =======================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ store2_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style
    p_a[0] = (char *)store2_p;

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ trun_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style
    p_a[1] = (char *)trun_p;

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    int * __restrict__ itndex_p = (int *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style
    p_a[2] = (char *)itndex_p;

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 15);
    ops_halo_exchanges(args, 15, range);
#endif

    if (block->instance->OPS_diags > 1) { 
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[544].mpi_time += __t2 - __t1;
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

//  call kernel wrapper function, passing in pointers to data
    if (x_size > 0 && y_size > 0 && z_size > 0) {

        ops_maths_kernel_eqBN<<<grid, tblock >>> (
                   (double *)p_a[0], xstride_0, ystride_0, zstride_0, 
                   (double *)p_a[1], xstride_1, ystride_1, zstride_1, 
                   (int *)p_a[2], xstride_2, ystride_2, zstride_2, 
                   (double *)arg3.data_d, 
                   (int *)arg4.data_d, 
                   (int *)arg5.data_d, 
                   (int *)arg6.data_d, 
                   (int *)arg7.data_d, 
                   (double *)arg8.data_d, 
                   *(int *)arg9.data, 
                   *(int *)arg10.data, 
                   *(int *)arg11.data, 
                   *(int *)arg12.data, 
                   *(int *)arg13.data, 
                   *(int *)arg14.data, 
                   x_size, y_size, z_size);

    }

    hipSafeCall(block->instance->ostream(), hipGetLastError());

    if(block->instance->OPS_diags > 1) {
        hipSafeCall(block->instance->ostream(), hipDeviceSynchronize());
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[544].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 15);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[544].mpi_time += __t2 - __t1;
        block->instance->OPS_kernels[544].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[544].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[544].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
    }
}

#ifdef OPS_LAZY
extern "C"
void maths_kernel_eqBN_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("maths_kernel_eqBN", args, 15, 544, dim, 1, range, block, maths_kernel_eqBN_execute);
}
#endif

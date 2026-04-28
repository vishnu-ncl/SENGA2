// Auto-generated at 2026-04-28 18:43:57.774888 by ops-translator

__constant__ int dims_bounds_kernel_eqAG_xdir[3][3];
static int dims_bounds_kernel_eqAG_xdir_h[3][3] = {{0}};

//  =============
//  User function
//  =============
__device__ void bounds_kernel_eqAG_xdir_gpu(ACC<double> &sorpx, const ACC<double> &strhx, const ACC<double> &ratex) {

    sorpx(0, 0, 0) = sorpx(0, 0, 0) + strhx(0, 0, 0) * ratex(0, 0, 0);

}

//  ============================
//  Cuda kernel wrapper function
//  ============================
__global__ void ops_bounds_kernel_eqAG_xdir(double* __restrict arg0, int xstride_0, int ystride_0, int zstride_0, 
double* __restrict arg1, int xstride_1, int ystride_1, int zstride_1, 
double* __restrict arg2, int xstride_2, int ystride_2, int zstride_2, 
int size0, int size1, int size2) {

    int idx_z = blockDim.z * blockIdx.z + threadIdx.z;
    int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
    int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

    arg0 += idx_x * xstride_0 + idx_y * ystride_0 * dims_bounds_kernel_eqAG_xdir[0][0] + idx_z * zstride_0 * dims_bounds_kernel_eqAG_xdir[0][0] * dims_bounds_kernel_eqAG_xdir[0][1];
    arg1 += idx_x * xstride_1 + idx_y * ystride_1 * dims_bounds_kernel_eqAG_xdir[1][0] + idx_z * zstride_1 * dims_bounds_kernel_eqAG_xdir[1][0] * dims_bounds_kernel_eqAG_xdir[1][1];
    arg2 += idx_x * xstride_2 + idx_y * ystride_2 * dims_bounds_kernel_eqAG_xdir[2][0] + idx_z * zstride_2 * dims_bounds_kernel_eqAG_xdir[2][0] * dims_bounds_kernel_eqAG_xdir[2][1];

    if(idx_x < size0 && idx_y < size1 && idx_z < size2) {

        ACC<double> argp0(dims_bounds_kernel_eqAG_xdir[0][0], dims_bounds_kernel_eqAG_xdir[0][1], arg0);
        const ACC<double> argp1(dims_bounds_kernel_eqAG_xdir[1][0], dims_bounds_kernel_eqAG_xdir[1][1], arg1);
        const ACC<double> argp2(dims_bounds_kernel_eqAG_xdir[2][0], dims_bounds_kernel_eqAG_xdir[2][1], arg2);

        bounds_kernel_eqAG_xdir_gpu(argp0, argp1, argp2);

    }// End of cuda index in_range check

}// End of cuda kernel wrapper function

//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void bounds_kernel_eqAG_xdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void bounds_kernel_eqAG_xdir_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[3];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
#endif

//  ======
//  Timing
//  ======
    double __t1, __t2, __c1, __c2;

    ops_arg arg0 = args[0];
    ops_arg arg1 = args[1];
    ops_arg arg2 = args[2];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 3, range, 361)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 361, "bounds_kernel_eqAG_xdir");
        block->instance->OPS_kernels[361].count++;
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

    if (xdim0 != dims_bounds_kernel_eqAG_xdir_h[0][0] || ydim0 != dims_bounds_kernel_eqAG_xdir_h[0][1] || zdim0 != dims_bounds_kernel_eqAG_xdir_h[0][2] || xdim1 != dims_bounds_kernel_eqAG_xdir_h[1][0] || ydim1 != dims_bounds_kernel_eqAG_xdir_h[1][1] || zdim1 != dims_bounds_kernel_eqAG_xdir_h[1][2] || xdim2 != dims_bounds_kernel_eqAG_xdir_h[2][0] || ydim2 != dims_bounds_kernel_eqAG_xdir_h[2][1] || zdim2 != dims_bounds_kernel_eqAG_xdir_h[2][2]) {
        dims_bounds_kernel_eqAG_xdir_h[0][0] = xdim0;
        dims_bounds_kernel_eqAG_xdir_h[0][1] = ydim0;
        dims_bounds_kernel_eqAG_xdir_h[0][2] = zdim0;
        dims_bounds_kernel_eqAG_xdir_h[1][0] = xdim1;
        dims_bounds_kernel_eqAG_xdir_h[1][1] = ydim1;
        dims_bounds_kernel_eqAG_xdir_h[1][2] = zdim1;
        dims_bounds_kernel_eqAG_xdir_h[2][0] = xdim2;
        dims_bounds_kernel_eqAG_xdir_h[2][1] = ydim2;
        dims_bounds_kernel_eqAG_xdir_h[2][2] = zdim2;

        cutilSafeCall(block->instance->ostream(), cudaMemcpyToSymbol( dims_bounds_kernel_eqAG_xdir, dims_bounds_kernel_eqAG_xdir_h, sizeof(dims_bounds_kernel_eqAG_xdir)));
    }

    int x_size = MAX(0,end_indx[0]-start_indx[0]+1);
    int y_size = MAX(0,end_indx[1]-start_indx[1]+1);
    int z_size = MAX(0,end_indx[2]-start_indx[2]+1);

    dim3 grid( (x_size-1)/block->instance->OPS_block_size_x + 1, (y_size-1)/block->instance->OPS_block_size_y + 1, (z_size-1)/block->instance->OPS_block_size_z + 1);

    dim3 tblock(block->instance->OPS_block_size_x,block->instance->OPS_block_size_y,block->instance->OPS_block_size_z);

    char *p_a[3];

//  =======================
//  Set up initial pointers
//  =======================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ sorpx_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style
    p_a[0] = (char *)sorpx_p;

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ strhx_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style
    p_a[1] = (char *)strhx_p;

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ ratex_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style
    p_a[2] = (char *)ratex_p;

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 3);
    ops_halo_exchanges(args, 3, range);
#endif

    if (block->instance->OPS_diags > 1) { 
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[361].mpi_time += __t2 - __t1;
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

        ops_bounds_kernel_eqAG_xdir<<<grid, tblock >>> (
                   (double *)p_a[0], xstride_0, ystride_0, zstride_0, 
                   (double *)p_a[1], xstride_1, ystride_1, zstride_1, 
                   (double *)p_a[2], xstride_2, ystride_2, zstride_2, 
                   x_size, y_size, z_size);

    }

    cutilSafeCall(block->instance->ostream(), cudaGetLastError());

    if(block->instance->OPS_diags > 1) {
        cutilSafeCall(block->instance->ostream(), cudaDeviceSynchronize());
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[361].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 3);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[361].mpi_time += __t2 - __t1;
        block->instance->OPS_kernels[361].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[361].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[361].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
    }
}

#ifdef OPS_LAZY
extern "C"
void bounds_kernel_eqAG_xdir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("bounds_kernel_eqAG_xdir", args, 3, 361, dim, 1, range, block, bounds_kernel_eqAG_xdir_execute);
}
#endif

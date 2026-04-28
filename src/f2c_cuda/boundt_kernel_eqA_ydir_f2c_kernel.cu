// Auto-generated at 2026-04-28 18:44:04.352443 by ops-translator

__constant__ int dims_boundt_kernel_eqA_ydir[7][3];
static int dims_boundt_kernel_eqA_ydir_h[7][3] = {{0}};

//  =============
//  User function
//  =============
__device__ void boundt_kernel_eqA_ydir_gpu(ACC<double> &urhs, ACC<double> &vrhs, ACC<double> &wrhs, const ACC<double> &drhs, const ACC<double> &struy, const ACC<double> &strvy, const ACC<double> &strwy) {

    urhs(0, 0, 0) = drhs(0, 0, 0) * struy(0, 0, 0);
    vrhs(0, 0, 0) = drhs(0, 0, 0) * strvy(0, 0, 0);
    wrhs(0, 0, 0) = drhs(0, 0, 0) * strwy(0, 0, 0);

}

//  ============================
//  Cuda kernel wrapper function
//  ============================
__global__ void ops_boundt_kernel_eqA_ydir(double* __restrict arg0, int xstride_0, int ystride_0, int zstride_0, 
double* __restrict arg1, int xstride_1, int ystride_1, int zstride_1, 
double* __restrict arg2, int xstride_2, int ystride_2, int zstride_2, 
double* __restrict arg3, int xstride_3, int ystride_3, int zstride_3, 
double* __restrict arg4, int xstride_4, int ystride_4, int zstride_4, 
double* __restrict arg5, int xstride_5, int ystride_5, int zstride_5, 
double* __restrict arg6, int xstride_6, int ystride_6, int zstride_6, 
int size0, int size1, int size2) {

    int idx_z = blockDim.z * blockIdx.z + threadIdx.z;
    int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
    int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

    arg0 += idx_x * xstride_0 + idx_y * ystride_0 * dims_boundt_kernel_eqA_ydir[0][0] + idx_z * zstride_0 * dims_boundt_kernel_eqA_ydir[0][0] * dims_boundt_kernel_eqA_ydir[0][1];
    arg1 += idx_x * xstride_1 + idx_y * ystride_1 * dims_boundt_kernel_eqA_ydir[1][0] + idx_z * zstride_1 * dims_boundt_kernel_eqA_ydir[1][0] * dims_boundt_kernel_eqA_ydir[1][1];
    arg2 += idx_x * xstride_2 + idx_y * ystride_2 * dims_boundt_kernel_eqA_ydir[2][0] + idx_z * zstride_2 * dims_boundt_kernel_eqA_ydir[2][0] * dims_boundt_kernel_eqA_ydir[2][1];
    arg3 += idx_x * xstride_3 + idx_y * ystride_3 * dims_boundt_kernel_eqA_ydir[3][0] + idx_z * zstride_3 * dims_boundt_kernel_eqA_ydir[3][0] * dims_boundt_kernel_eqA_ydir[3][1];
    arg4 += idx_x * xstride_4 + idx_y * ystride_4 * dims_boundt_kernel_eqA_ydir[4][0] + idx_z * zstride_4 * dims_boundt_kernel_eqA_ydir[4][0] * dims_boundt_kernel_eqA_ydir[4][1];
    arg5 += idx_x * xstride_5 + idx_y * ystride_5 * dims_boundt_kernel_eqA_ydir[5][0] + idx_z * zstride_5 * dims_boundt_kernel_eqA_ydir[5][0] * dims_boundt_kernel_eqA_ydir[5][1];
    arg6 += idx_x * xstride_6 + idx_y * ystride_6 * dims_boundt_kernel_eqA_ydir[6][0] + idx_z * zstride_6 * dims_boundt_kernel_eqA_ydir[6][0] * dims_boundt_kernel_eqA_ydir[6][1];

    if(idx_x < size0 && idx_y < size1 && idx_z < size2) {

        ACC<double> argp0(dims_boundt_kernel_eqA_ydir[0][0], dims_boundt_kernel_eqA_ydir[0][1], arg0);
        ACC<double> argp1(dims_boundt_kernel_eqA_ydir[1][0], dims_boundt_kernel_eqA_ydir[1][1], arg1);
        ACC<double> argp2(dims_boundt_kernel_eqA_ydir[2][0], dims_boundt_kernel_eqA_ydir[2][1], arg2);
        const ACC<double> argp3(dims_boundt_kernel_eqA_ydir[3][0], dims_boundt_kernel_eqA_ydir[3][1], arg3);
        const ACC<double> argp4(dims_boundt_kernel_eqA_ydir[4][0], dims_boundt_kernel_eqA_ydir[4][1], arg4);
        const ACC<double> argp5(dims_boundt_kernel_eqA_ydir[5][0], dims_boundt_kernel_eqA_ydir[5][1], arg5);
        const ACC<double> argp6(dims_boundt_kernel_eqA_ydir[6][0], dims_boundt_kernel_eqA_ydir[6][1], arg6);

        boundt_kernel_eqA_ydir_gpu(argp0, argp1, argp2, argp3, argp4, argp5, argp6);

    }// End of cuda index in_range check

}// End of cuda kernel wrapper function

//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void boundt_kernel_eqA_ydir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void boundt_kernel_eqA_ydir_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[7];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
    args[4] = desc->args[4];
    args[5] = desc->args[5];
    args[6] = desc->args[6];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 7, range, 524)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 524, "boundt_kernel_eqA_ydir");
        block->instance->OPS_kernels[524].count++;
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

    if (xdim0 != dims_boundt_kernel_eqA_ydir_h[0][0] || ydim0 != dims_boundt_kernel_eqA_ydir_h[0][1] || zdim0 != dims_boundt_kernel_eqA_ydir_h[0][2] || xdim1 != dims_boundt_kernel_eqA_ydir_h[1][0] || ydim1 != dims_boundt_kernel_eqA_ydir_h[1][1] || zdim1 != dims_boundt_kernel_eqA_ydir_h[1][2] || xdim2 != dims_boundt_kernel_eqA_ydir_h[2][0] || ydim2 != dims_boundt_kernel_eqA_ydir_h[2][1] || zdim2 != dims_boundt_kernel_eqA_ydir_h[2][2] || xdim3 != dims_boundt_kernel_eqA_ydir_h[3][0] || ydim3 != dims_boundt_kernel_eqA_ydir_h[3][1] || zdim3 != dims_boundt_kernel_eqA_ydir_h[3][2] || xdim4 != dims_boundt_kernel_eqA_ydir_h[4][0] || ydim4 != dims_boundt_kernel_eqA_ydir_h[4][1] || zdim4 != dims_boundt_kernel_eqA_ydir_h[4][2] || xdim5 != dims_boundt_kernel_eqA_ydir_h[5][0] || ydim5 != dims_boundt_kernel_eqA_ydir_h[5][1] || zdim5 != dims_boundt_kernel_eqA_ydir_h[5][2] || xdim6 != dims_boundt_kernel_eqA_ydir_h[6][0] || ydim6 != dims_boundt_kernel_eqA_ydir_h[6][1] || zdim6 != dims_boundt_kernel_eqA_ydir_h[6][2]) {
        dims_boundt_kernel_eqA_ydir_h[0][0] = xdim0;
        dims_boundt_kernel_eqA_ydir_h[0][1] = ydim0;
        dims_boundt_kernel_eqA_ydir_h[0][2] = zdim0;
        dims_boundt_kernel_eqA_ydir_h[1][0] = xdim1;
        dims_boundt_kernel_eqA_ydir_h[1][1] = ydim1;
        dims_boundt_kernel_eqA_ydir_h[1][2] = zdim1;
        dims_boundt_kernel_eqA_ydir_h[2][0] = xdim2;
        dims_boundt_kernel_eqA_ydir_h[2][1] = ydim2;
        dims_boundt_kernel_eqA_ydir_h[2][2] = zdim2;
        dims_boundt_kernel_eqA_ydir_h[3][0] = xdim3;
        dims_boundt_kernel_eqA_ydir_h[3][1] = ydim3;
        dims_boundt_kernel_eqA_ydir_h[3][2] = zdim3;
        dims_boundt_kernel_eqA_ydir_h[4][0] = xdim4;
        dims_boundt_kernel_eqA_ydir_h[4][1] = ydim4;
        dims_boundt_kernel_eqA_ydir_h[4][2] = zdim4;
        dims_boundt_kernel_eqA_ydir_h[5][0] = xdim5;
        dims_boundt_kernel_eqA_ydir_h[5][1] = ydim5;
        dims_boundt_kernel_eqA_ydir_h[5][2] = zdim5;
        dims_boundt_kernel_eqA_ydir_h[6][0] = xdim6;
        dims_boundt_kernel_eqA_ydir_h[6][1] = ydim6;
        dims_boundt_kernel_eqA_ydir_h[6][2] = zdim6;

        cutilSafeCall(block->instance->ostream(), cudaMemcpyToSymbol( dims_boundt_kernel_eqA_ydir, dims_boundt_kernel_eqA_ydir_h, sizeof(dims_boundt_kernel_eqA_ydir)));
    }

    int x_size = MAX(0,end_indx[0]-start_indx[0]+1);
    int y_size = MAX(0,end_indx[1]-start_indx[1]+1);
    int z_size = MAX(0,end_indx[2]-start_indx[2]+1);

    dim3 grid( (x_size-1)/block->instance->OPS_block_size_x + 1, (y_size-1)/block->instance->OPS_block_size_y + 1, (z_size-1)/block->instance->OPS_block_size_z + 1);

    dim3 tblock(block->instance->OPS_block_size_x,block->instance->OPS_block_size_y,block->instance->OPS_block_size_z);

    char *p_a[7];

//  =======================
//  Set up initial pointers
//  =======================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ urhs_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style
    p_a[0] = (char *)urhs_p;

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ vrhs_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style
    p_a[1] = (char *)vrhs_p;

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ wrhs_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style
    p_a[2] = (char *)wrhs_p;

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style
    p_a[3] = (char *)drhs_p;

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ struy_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style
    p_a[4] = (char *)struy_p;

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ strvy_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style
    p_a[5] = (char *)strvy_p;

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ strwy_p = (double *)(args[6].data_d) + base6 - 1; // Subtracting 1 to convert to C-style
    p_a[6] = (char *)strwy_p;

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 7);
    ops_halo_exchanges(args, 7, range);
#endif

    if (block->instance->OPS_diags > 1) { 
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[524].mpi_time += __t2 - __t1;
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

//  call kernel wrapper function, passing in pointers to data
    if (x_size > 0 && y_size > 0 && z_size > 0) {

        ops_boundt_kernel_eqA_ydir<<<grid, tblock >>> (
                   (double *)p_a[0], xstride_0, ystride_0, zstride_0, 
                   (double *)p_a[1], xstride_1, ystride_1, zstride_1, 
                   (double *)p_a[2], xstride_2, ystride_2, zstride_2, 
                   (double *)p_a[3], xstride_3, ystride_3, zstride_3, 
                   (double *)p_a[4], xstride_4, ystride_4, zstride_4, 
                   (double *)p_a[5], xstride_5, ystride_5, zstride_5, 
                   (double *)p_a[6], xstride_6, ystride_6, zstride_6, 
                   x_size, y_size, z_size);

    }

    cutilSafeCall(block->instance->ostream(), cudaGetLastError());

    if(block->instance->OPS_diags > 1) {
        cutilSafeCall(block->instance->ostream(), cudaDeviceSynchronize());
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[524].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 7);
    ops_set_halo_dirtybit3(&args[0], range);
    ops_set_halo_dirtybit3(&args[1], range);
    ops_set_halo_dirtybit3(&args[2], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[524].mpi_time += __t2 - __t1;
        block->instance->OPS_kernels[524].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[524].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[524].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[524].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[524].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[524].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[524].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
    }
}

#ifdef OPS_LAZY
extern "C"
void boundt_kernel_eqA_ydir_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("boundt_kernel_eqA_ydir", args, 7, 524, dim, 1, range, block, boundt_kernel_eqA_ydir_execute);
}
#endif

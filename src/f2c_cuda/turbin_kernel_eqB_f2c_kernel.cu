// Auto-generated at 2026-04-28 18:44:05.131254 by ops-translator

__constant__ int dims_turbin_kernel_eqB[2][3];
static int dims_turbin_kernel_eqB_h[2][3] = {{0}};

//  =============
//  User function
//  =============
__device__ void turbin_kernel_eqB_gpu(const ACC<double> &in_arr, double *reduction_var) {

    *reduction_var = *reduction_var + in_arr(0, 0, 0);

}

//  ============================
//  Cuda kernel wrapper function
//  ============================
__global__ void ops_turbin_kernel_eqB(double* __restrict arg0, int xstride_0, int ystride_0, int zstride_0, 
double* __restrict arg1, 
int size0, int size1, int size2) {
    double arg1_l[1];
    for (int d = 0; d < 1; d++) arg1_l[d] = ZERO_double;

    int idx_z = blockDim.z * blockIdx.z + threadIdx.z;
    int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
    int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

    arg0 += idx_x * xstride_0 + idx_y * ystride_0 * dims_turbin_kernel_eqB[0][0] + idx_z * zstride_0 * dims_turbin_kernel_eqB[0][0] * dims_turbin_kernel_eqB[0][1];

    if(idx_x < size0 && idx_y < size1 && idx_z < size2) {

        const ACC<double> argp0(dims_turbin_kernel_eqB[0][0], dims_turbin_kernel_eqB[0][1], arg0);

        turbin_kernel_eqB_gpu(argp0, arg1_l);

    }// End of cuda index in_range check

//  ==============================
//  Reduction across thread blocks
//  ==============================
    for(int d = 0; d < 1; d++)
        ops_reduction_cuda<OPS_INC>(&arg1[d+(blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y)*1],arg1_l[d]);

}// End of cuda kernel wrapper function

//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void turbin_kernel_eqB_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void turbin_kernel_eqB_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 2, range, 556)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 556, "turbin_kernel_eqB");
        block->instance->OPS_kernels[556].count++;
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

    if (xdim0 != dims_turbin_kernel_eqB_h[0][0] || ydim0 != dims_turbin_kernel_eqB_h[0][1] || zdim0 != dims_turbin_kernel_eqB_h[0][2]) {
        dims_turbin_kernel_eqB_h[0][0] = xdim0;
        dims_turbin_kernel_eqB_h[0][1] = ydim0;
        dims_turbin_kernel_eqB_h[0][2] = zdim0;

        cutilSafeCall(block->instance->ostream(), cudaMemcpyToSymbol( dims_turbin_kernel_eqB, dims_turbin_kernel_eqB_h, sizeof(dims_turbin_kernel_eqB)));
    }

#ifdef OPS_MPI
    double *arg1h = (double*)(((ops_reduction)args[1].data)->data + ((ops_reduction)args[1].data)->size*block->index);
#else
    double *arg1h = (double*)(((ops_reduction)args[1].data)->data);
#endif

    int x_size = MAX(0,end_indx[0]-start_indx[0]+1);
    int y_size = MAX(0,end_indx[1]-start_indx[1]+1);
    int z_size = MAX(0,end_indx[2]-start_indx[2]+1);

    dim3 grid( (x_size-1)/block->instance->OPS_block_size_x + 1, (y_size-1)/block->instance->OPS_block_size_y + 1, (z_size-1)/block->instance->OPS_block_size_z + 1);

    dim3 tblock(block->instance->OPS_block_size_x,block->instance->OPS_block_size_y,block->instance->OPS_block_size_z);

    int nblocks = ((x_size-1)/block->instance->OPS_block_size_x + 1)*((y_size-1)/block->instance->OPS_block_size_y + 1)*((z_size-1)/block->instance->OPS_block_size_z + 1);

    int maxblocks = nblocks;
    int reduct_bytes = 0;
    size_t reduct_size = 0;

    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    reduct_size = MAX(reduct_size,1*sizeof(double));

    reallocReductArrays(block->instance, reduct_bytes);
    reduct_bytes = 0;

    arg1.data = block->instance->OPS_reduct_h + reduct_bytes;
    arg1.data_d = block->instance->OPS_reduct_d + reduct_bytes;
    for (int b = 0; b < maxblocks; b++) {
        for (int d = 0; d < 1; d++)   ((double *)arg1.data)[d+b*1] = ZERO_double;
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));

    mvReductArraysToDevice(block->instance, reduct_bytes);
    char *p_a[2];

//  =======================
//  Set up initial pointers
//  =======================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ in_arr_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style
    p_a[0] = (char *)in_arr_p;

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 2);
    ops_halo_exchanges(args, 2, range);
#endif

    if (block->instance->OPS_diags > 1) { 
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[556].mpi_time += __t2 - __t1;
    }

    size_t nshared = 0;
    int nthread = block->instance->OPS_block_size_x*block->instance->OPS_block_size_y*block->instance->OPS_block_size_z;

    nshared = MAX(nshared,sizeof(double)*1);

    nshared = MAX(nshared*nthread,reduct_size*nthread);

//  ==========================================================
//  ops_dat strides for offset calculation in wrapper function
//  ==========================================================
    int xstride_0, ystride_0, zstride_0;
    xstride_0 = args[0].stencil->stride[0];    ystride_0 = args[0].stencil->stride[1];
    zstride_0 = args[0].stencil->stride[2];

//  call kernel wrapper function, passing in pointers to data
    if (x_size > 0 && y_size > 0 && z_size > 0) {

        ops_turbin_kernel_eqB<<<grid, tblock, nshared >>> (
                   (double *)p_a[0], xstride_0, ystride_0, zstride_0, 
                   (double *)arg1.data_d, 
                   x_size, y_size, z_size);

    }

    cutilSafeCall(block->instance->ostream(), cudaGetLastError());

    mvReductArraysToHost(block->instance, reduct_bytes);

    for (int b = 0; b < maxblocks; b++)
        for (int d = 0; d < 1; d++)
           arg1h[d] = arg1h[d] + ((double *)arg1.data)[d+b*1];
    arg1.data = (char *)arg1h;

    if(block->instance->OPS_diags > 1) {
        cutilSafeCall(block->instance->ostream(), cudaDeviceSynchronize());
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[556].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 2);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[556].mpi_time += __t2 - __t1;
        block->instance->OPS_kernels[556].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
    }
}

#ifdef OPS_LAZY
extern "C"
void turbin_kernel_eqB_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("turbin_kernel_eqB", args, 2, 556, dim, 1, range, block, turbin_kernel_eqB_execute);
}
#endif

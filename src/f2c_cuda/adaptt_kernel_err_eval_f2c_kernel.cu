// Auto-generated at 2026-04-28 18:43:56.227183 by ops-translator

__constant__ int dims_adaptt_kernel_err_eval[4][3];
static int dims_adaptt_kernel_err_eval_h[4][3] = {{0}};

//  =============
//  User function
//  =============
__device__ void adaptt_kernel_err_eval_gpu(const ACC<double> &err_arr, const ACC<double> &run_arr, const double *ernrm, double *ertot) {

    double fornow;

    fornow = f2c::abs(err_arr(0, 0, 0)) / (f2c::abs(run_arr(0, 0, 0)) + ernrm[0]);
    *ertot = f2c::max(*ertot, fornow);

}

//  ============================
//  Cuda kernel wrapper function
//  ============================
__global__ void ops_adaptt_kernel_err_eval(double* __restrict arg0, int xstride_0, int ystride_0, int zstride_0, 
double* __restrict arg1, int xstride_1, int ystride_1, int zstride_1, 
const double arg2, 
double* __restrict arg3, 
int size0, int size1, int size2) {
    double arg3_l[1];
    for (int d = 0; d < 1; d++) arg3_l[d] = -INFINITY_double;

    int idx_z = blockDim.z * blockIdx.z + threadIdx.z;
    int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
    int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

    arg0 += idx_x * xstride_0 + idx_y * ystride_0 * dims_adaptt_kernel_err_eval[0][0] + idx_z * zstride_0 * dims_adaptt_kernel_err_eval[0][0] * dims_adaptt_kernel_err_eval[0][1];
    arg1 += idx_x * xstride_1 + idx_y * ystride_1 * dims_adaptt_kernel_err_eval[1][0] + idx_z * zstride_1 * dims_adaptt_kernel_err_eval[1][0] * dims_adaptt_kernel_err_eval[1][1];

    if(idx_x < size0 && idx_y < size1 && idx_z < size2) {

        const ACC<double> argp0(dims_adaptt_kernel_err_eval[0][0], dims_adaptt_kernel_err_eval[0][1], arg0);
        const ACC<double> argp1(dims_adaptt_kernel_err_eval[1][0], dims_adaptt_kernel_err_eval[1][1], arg1);

        adaptt_kernel_err_eval_gpu(argp0, argp1, &arg2, arg3_l);

    }// End of cuda index in_range check

//  ==============================
//  Reduction across thread blocks
//  ==============================
    for(int d = 0; d < 1; d++)
        ops_reduction_cuda<OPS_MAX>(&arg3[d+(blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y)*1],arg3_l[d]);

}// End of cuda kernel wrapper function

//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void adaptt_kernel_err_eval_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void adaptt_kernel_err_eval_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[4];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
#endif

//  ======
//  Timing
//  ======
    double __t1, __t2, __c1, __c2;

    ops_arg arg0 = args[0];
    ops_arg arg1 = args[1];
    ops_arg arg2 = args[2];
    ops_arg arg3 = args[3];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 4, range, 315)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 315, "adaptt_kernel_err_eval");
        block->instance->OPS_kernels[315].count++;
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

    if (xdim0 != dims_adaptt_kernel_err_eval_h[0][0] || ydim0 != dims_adaptt_kernel_err_eval_h[0][1] || zdim0 != dims_adaptt_kernel_err_eval_h[0][2] || xdim1 != dims_adaptt_kernel_err_eval_h[1][0] || ydim1 != dims_adaptt_kernel_err_eval_h[1][1] || zdim1 != dims_adaptt_kernel_err_eval_h[1][2]) {
        dims_adaptt_kernel_err_eval_h[0][0] = xdim0;
        dims_adaptt_kernel_err_eval_h[0][1] = ydim0;
        dims_adaptt_kernel_err_eval_h[0][2] = zdim0;
        dims_adaptt_kernel_err_eval_h[1][0] = xdim1;
        dims_adaptt_kernel_err_eval_h[1][1] = ydim1;
        dims_adaptt_kernel_err_eval_h[1][2] = zdim1;

        cutilSafeCall(block->instance->ostream(), cudaMemcpyToSymbol( dims_adaptt_kernel_err_eval, dims_adaptt_kernel_err_eval_h, sizeof(dims_adaptt_kernel_err_eval)));
    }

#ifdef OPS_MPI
    double *arg3h = (double*)(((ops_reduction)args[3].data)->data + ((ops_reduction)args[3].data)->size*block->index);
#else
    double *arg3h = (double*)(((ops_reduction)args[3].data)->data);
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

    arg3.data = block->instance->OPS_reduct_h + reduct_bytes;
    arg3.data_d = block->instance->OPS_reduct_d + reduct_bytes;
    for (int b = 0; b < maxblocks; b++) {
        for (int d = 0; d < 1; d++)   ((double *)arg3.data)[d+b*1] = -INFINITY_double;
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));

    mvReductArraysToDevice(block->instance, reduct_bytes);
    char *p_a[4];

//  =======================
//  Set up initial pointers
//  =======================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ err_arr_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style
    p_a[0] = (char *)err_arr_p;

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ run_arr_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style
    p_a[1] = (char *)run_arr_p;

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 4);
    ops_halo_exchanges(args, 4, range);
#endif

    if (block->instance->OPS_diags > 1) { 
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[315].mpi_time += __t2 - __t1;
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
    int xstride_1, ystride_1, zstride_1;
    xstride_1 = args[1].stencil->stride[0];    ystride_1 = args[1].stencil->stride[1];
    zstride_1 = args[1].stencil->stride[2];

//  call kernel wrapper function, passing in pointers to data
    if (x_size > 0 && y_size > 0 && z_size > 0) {

        ops_adaptt_kernel_err_eval<<<grid, tblock, nshared >>> (
                   (double *)p_a[0], xstride_0, ystride_0, zstride_0, 
                   (double *)p_a[1], xstride_1, ystride_1, zstride_1, 
                   *(double *)arg2.data, 
                   (double *)arg3.data_d, 
                   x_size, y_size, z_size);

    }

    cutilSafeCall(block->instance->ostream(), cudaGetLastError());

    mvReductArraysToHost(block->instance, reduct_bytes);

    for (int b = 0; b < maxblocks; b++)
        for (int d = 0; d < 1; d++)
            arg3h[d] = MAX(arg3h[d],((double *)arg3.data)[d+b*1]);
    arg3.data = (char *)arg3h;

    if(block->instance->OPS_diags > 1) {
        cutilSafeCall(block->instance->ostream(), cudaDeviceSynchronize());
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[315].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 4);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[315].mpi_time += __t2 - __t1;
        block->instance->OPS_kernels[315].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[315].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void adaptt_kernel_err_eval_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("adaptt_kernel_err_eval", args, 4, 315, dim, 1, range, block, adaptt_kernel_err_eval_execute);
}
#endif

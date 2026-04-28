// Auto-generated at 2026-04-28 18:44:06.246313 by ops-translator

__constant__ int dims_d2fdx2_kernel_main[7][3];
static int dims_d2fdx2_kernel_main_h[7][3] = {{0}};

//  =============
//  User function
//  =============
__device__ void d2fdx2_kernel_main_gpu(const ACC<double> &functn, ACC<double> &fderiv, const int *nxglbl, const int *nendxl, const int *nendxr, const int *nbound, const int *idx) {

    double fdifap;
    double fdifbp;
    double fdifcp;
    double fdifdp;
    double fdifep;
    double fdifam;
    double fdifbm;
    double fdifcm;
    double fdifdm;
    double fdifem;
    int ic;
    int istart;
    int ifinis;

    istart = 1;
    ifinis = nxglbl[0];
    if (nendxl[0] == nbound[0]) {
        istart = 6;
    }
    if (nendxr[0] == nbound[0]) {
        ifinis = nxglbl[0] - 5;
    }
    ic = idx[0];
    if (ic >= istart && ic <= ifinis) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(-3, 0, 0);
        fdifdp = functn(4, 0, 0) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(-4, 0, 0);
        fdifep = functn(5, 0, 0) - functn(0, 0, 0);
        fdifem = functn(0, 0, 0) - functn(-5, 0, 0);
        fderiv(0, 0, 0) = acofsx * (fdifap - fdifam) + bcofsx * (fdifbp - fdifbm) + ccofsx * (fdifcp - fdifcm) + dcofsx * (fdifdp - fdifdm) + ecofsx * (fdifep - fdifem);
    } else if (ic == 1) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifdp = functn(4, 0, 0) - functn(0, 0, 0);
        fdifep = functn(5, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs1x * fdifap + bcfs1x * fdifbp + ccfs1x * fdifcp + dcfs1x * fdifdp + ecfs1x * fdifep;
    } else if (ic == 2) {
        fdifap = functn(-1, 0, 0) - functn(0, 0, 0);
        fdifbp = functn(1, 0, 0) - functn(0, 0, 0);
        fdifcp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifdp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifep = functn(4, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs2x * fdifap + bcfs2x * fdifbp + ccfs2x * fdifcp + dcfs2x * fdifdp + ecfs2x * fdifep;
    } else if (ic == 3) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fderiv(0, 0, 0) = acfs3x * (fdifap - fdifam) + bcfs3x * (fdifbp - fdifbm);
    } else if (ic == 4) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(-3, 0, 0);
        fderiv(0, 0, 0) = acfs4x * (fdifap - fdifam) + bcfs4x * (fdifbp - fdifbm) + ccfs4x * (fdifcp - fdifcm);
    } else if (ic == 5) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(-3, 0, 0);
        fdifdp = functn(4, 0, 0) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(-4, 0, 0);
        fderiv(0, 0, 0) = acfs5x * (fdifap - fdifam) + bcfs5x * (fdifbp - fdifbm) + ccfs5x * (fdifcp - fdifcm) + dcfs5x * (fdifdp - fdifdm);
    } else if (ic == nxglbl[0] - 4) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(-3, 0, 0);
        fdifdp = functn(4, 0, 0) - functn(0, 0, 0);
        fdifdm = functn(0, 0, 0) - functn(-4, 0, 0);
        fderiv(0, 0, 0) = acfs5x * (fdifap - fdifam) + bcfs5x * (fdifbp - fdifbm) + ccfs5x * (fdifcp - fdifcm) + dcfs5x * (fdifdp - fdifdm);
    } else if (ic == nxglbl[0] - 3) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fdifcp = functn(3, 0, 0) - functn(0, 0, 0);
        fdifcm = functn(0, 0, 0) - functn(-3, 0, 0);
        fderiv(0, 0, 0) = acfs4x * (fdifap - fdifam) + bcfs4x * (fdifbp - fdifbm) + ccfs4x * (fdifcp - fdifcm);
    } else if (ic == nxglbl[0] - 2) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifam = functn(0, 0, 0) - functn(-1, 0, 0);
        fdifbp = functn(2, 0, 0) - functn(0, 0, 0);
        fdifbm = functn(0, 0, 0) - functn(-2, 0, 0);
        fderiv(0, 0, 0) = acfs3x * (fdifap - fdifam) + bcfs3x * (fdifbp - fdifbm);
    } else if (ic == nxglbl[0] - 1) {
        fdifap = functn(1, 0, 0) - functn(0, 0, 0);
        fdifbp = functn(-1, 0, 0) - functn(0, 0, 0);
        fdifcp = functn(-2, 0, 0) - functn(0, 0, 0);
        fdifdp = functn(-3, 0, 0) - functn(0, 0, 0);
        fdifep = functn(-4, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs2x * fdifap + bcfs2x * fdifbp + ccfs2x * fdifcp + dcfs2x * fdifdp + ecfs2x * fdifep;
    } else if (ic == nxglbl[0]) {
        fdifap = functn(-1, 0, 0) - functn(0, 0, 0);
        fdifbp = functn(-2, 0, 0) - functn(0, 0, 0);
        fdifcp = functn(-3, 0, 0) - functn(0, 0, 0);
        fdifdp = functn(-4, 0, 0) - functn(0, 0, 0);
        fdifep = functn(-5, 0, 0) - functn(0, 0, 0);
        fderiv(0, 0, 0) = acfs1x * fdifap + bcfs1x * fdifbp + ccfs1x * fdifcp + dcfs1x * fdifdp + ecfs1x * fdifep;
    }
    fderiv(0, 0, 0) = fderiv(0, 0, 0) * ovdlx2;

}

//  ============================
//  Cuda kernel wrapper function
//  ============================
__global__ void ops_d2fdx2_kernel_main(double* __restrict arg0, int xstride_0, int ystride_0, int zstride_0, 
double* __restrict arg1, int xstride_1, int ystride_1, int zstride_1, 
const int arg2, 
const int arg3, 
const int arg4, 
const int arg5, 
int arg_idx0, int arg_idx1, int arg_idx2, 
int size0, int size1, int size2) {

    int idx_z = blockDim.z * blockIdx.z + threadIdx.z;
    int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
    int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

    int arg_idx[3];
    arg_idx[0] = arg_idx0+idx_x;
    arg_idx[1] = arg_idx1+idx_y;
    arg_idx[2] = arg_idx2+idx_z;

    arg0 += idx_x * xstride_0 + idx_y * ystride_0 * dims_d2fdx2_kernel_main[0][0] + idx_z * zstride_0 * dims_d2fdx2_kernel_main[0][0] * dims_d2fdx2_kernel_main[0][1];
    arg1 += idx_x * xstride_1 + idx_y * ystride_1 * dims_d2fdx2_kernel_main[1][0] + idx_z * zstride_1 * dims_d2fdx2_kernel_main[1][0] * dims_d2fdx2_kernel_main[1][1];

    if(idx_x < size0 && idx_y < size1 && idx_z < size2) {

        const ACC<double> argp0(dims_d2fdx2_kernel_main[0][0], dims_d2fdx2_kernel_main[0][1], arg0);
        ACC<double> argp1(dims_d2fdx2_kernel_main[1][0], dims_d2fdx2_kernel_main[1][1], arg1);

        d2fdx2_kernel_main_gpu(argp0, argp1, &arg2, &arg3, &arg4, &arg5, arg_idx);

    }// End of cuda index in_range check

}// End of cuda kernel wrapper function

//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void d2fdx2_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void d2fdx2_kernel_main_execute(ops_kernel_descriptor *desc)
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
    if (!ops_checkpointing_before(args, 7, range, 7)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 7, "d2fdx2_kernel_main");
        block->instance->OPS_kernels[7].count++;
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

    if (xdim0 != dims_d2fdx2_kernel_main_h[0][0] || ydim0 != dims_d2fdx2_kernel_main_h[0][1] || zdim0 != dims_d2fdx2_kernel_main_h[0][2] || xdim1 != dims_d2fdx2_kernel_main_h[1][0] || ydim1 != dims_d2fdx2_kernel_main_h[1][1] || zdim1 != dims_d2fdx2_kernel_main_h[1][2]) {
        dims_d2fdx2_kernel_main_h[0][0] = xdim0;
        dims_d2fdx2_kernel_main_h[0][1] = ydim0;
        dims_d2fdx2_kernel_main_h[0][2] = zdim0;
        dims_d2fdx2_kernel_main_h[1][0] = xdim1;
        dims_d2fdx2_kernel_main_h[1][1] = ydim1;
        dims_d2fdx2_kernel_main_h[1][2] = zdim1;

        hipSafeCall(block->instance->ostream(), hipMemcpyToSymbol( dims_d2fdx2_kernel_main, dims_d2fdx2_kernel_main_h, sizeof(dims_d2fdx2_kernel_main)));
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
    double * __restrict__ functn_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style
    p_a[0] = (char *)functn_p;

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ fderiv_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style
    p_a[1] = (char *)fderiv_p;

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 7);
    ops_halo_exchanges(args, 7, range);
#endif

    if (block->instance->OPS_diags > 1) { 
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[7].mpi_time += __t2 - __t1;
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

//  call kernel wrapper function, passing in pointers to data
    if (x_size > 0 && y_size > 0 && z_size > 0) {

        ops_d2fdx2_kernel_main<<<grid, tblock >>> (
                   (double *)p_a[0], xstride_0, ystride_0, zstride_0, 
                   (double *)p_a[1], xstride_1, ystride_1, zstride_1, 
                   *(int *)arg2.data, 
                   *(int *)arg3.data, 
                   *(int *)arg4.data, 
                   *(int *)arg5.data, 
                   arg_idx[0], arg_idx[1], arg_idx[2], 
                   x_size, y_size, z_size);

    }

    hipSafeCall(block->instance->ostream(), hipGetLastError());

    if(block->instance->OPS_diags > 1) {
        hipSafeCall(block->instance->ostream(), hipDeviceSynchronize());
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[7].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 7);
    ops_set_halo_dirtybit3(&args[1], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[7].mpi_time += __t2 - __t1;
        block->instance->OPS_kernels[7].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[7].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
    }
}

#ifdef OPS_LAZY
extern "C"
void d2fdx2_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("d2fdx2_kernel_main", args, 7, 7, dim, 1, range, block, d2fdx2_kernel_main_execute);
}
#endif

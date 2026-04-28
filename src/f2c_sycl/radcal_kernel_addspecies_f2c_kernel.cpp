// Auto-generated at 2026-04-28 18:44:58.588894 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void radcal_kernel_addspecies_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void radcal_kernel_addspecies_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[6];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
    args[4] = desc->args[4];
    args[5] = desc->args[5];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 6, range, 535)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 535, "radcal_kernel_addspecies");
        block->instance->OPS_kernels[535].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "radcal_kernel_addspecies");
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
    int xdim0_radcal_kernel_addspecies = args[0].dat->size[0];
    int ydim0_radcal_kernel_addspecies = args[0].dat->size[1];
    int xdim1_radcal_kernel_addspecies = args[1].dat->size[0];
    int ydim1_radcal_kernel_addspecies = args[1].dat->size[1];
    int xdim2_radcal_kernel_addspecies = args[2].dat->size[0];
    int ydim2_radcal_kernel_addspecies = args[2].dat->size[1];
    int xdim3_radcal_kernel_addspecies = args[3].dat->size[0];
    int ydim3_radcal_kernel_addspecies = args[3].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ store1_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ yrhs_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ trun_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ store2_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    double *arg4h = (double *)args[4].data;

    int ispec_val = *(int *)args[5].data;

//  Subtracting 1 here as start_indx and end_indx is in Fortran style - converting it to c-style range
    int start_0 = start_indx[0]-1;
    int end_0 = end_indx[0];
    int start_1 = start_indx[1]-1;
    int end_1 = end_indx[1];
    int start_2 = start_indx[2]-1;
    int end_2 = end_indx[2];

    int consts_bytes = 0;

    consts_bytes += ROUND_UP(arg4.dim*sizeof(double));

    reallocConstArrays(block->instance, consts_bytes);
    consts_bytes = 0;

    arg4.data = block->instance->OPS_consts_h + consts_bytes;
    double* arg4_data_d = (double*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg4.dim; d++)     ((double *)arg4.data)[d] = arg4h[d];
    consts_bytes += ROUND_UP(arg4.dim*sizeof(double));

    mvConstArraysToDevice(block->instance, consts_bytes);

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 6);
    ops_halo_exchanges(args, 6, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[535].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            cgh.parallel_for<class radcal_kernel_addspecies_sycl>(
                cl::sycl::nd_range<3>(
                    cl::sycl::range<3>(
                        ((end_2-start_2-1)/block->instance->OPS_block_size_z+1)*block->instance->OPS_block_size_z,
                        ((end_1-start_1-1)/block->instance->OPS_block_size_y+1)*block->instance->OPS_block_size_y,
                        ((end_0-start_0-1)/block->instance->OPS_block_size_x+1)*block->instance->OPS_block_size_x),
                    cl::sycl::range<3>(
                        block->instance->OPS_block_size_z,
                        block->instance->OPS_block_size_y,
                        block->instance->OPS_block_size_x)
            )
            , [=](cl::sycl::nd_item<3> item
            ) [[intel::kernel_args_restrict]] {

                int n_z = item.get_global_id()[0];
                int n_y = item.get_global_id()[1];
                int n_x = item.get_global_id()[2];

                 ACC<double> store1(xdim0_radcal_kernel_addspecies, ydim0_radcal_kernel_addspecies, store1_p + (n_x * 1) + (n_y * xdim0_radcal_kernel_addspecies * 1) + (n_z * xdim0_radcal_kernel_addspecies * ydim0_radcal_kernel_addspecies * 1));
                const  ACC<double> yrhs(xdim1_radcal_kernel_addspecies, ydim1_radcal_kernel_addspecies, yrhs_p + (n_x * 1) + (n_y * xdim1_radcal_kernel_addspecies * 1) + (n_z * xdim1_radcal_kernel_addspecies * ydim1_radcal_kernel_addspecies * 1));
                const  ACC<double> trun(xdim2_radcal_kernel_addspecies, ydim2_radcal_kernel_addspecies, trun_p + (n_x * 1) + (n_y * xdim2_radcal_kernel_addspecies * 1) + (n_z * xdim2_radcal_kernel_addspecies * ydim2_radcal_kernel_addspecies * 1));
                const  ACC<double> store2(xdim3_radcal_kernel_addspecies, ydim3_radcal_kernel_addspecies, store2_p + (n_x * 1) + (n_y * xdim3_radcal_kernel_addspecies * 1) + (n_z * xdim3_radcal_kernel_addspecies * ydim3_radcal_kernel_addspecies * 1));

                const double *rgspec = arg4_data_d;
                const int *ispec = &ispec_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double fornow;

    fornow = yrhs(0, 0, 0) * rgspec[ispec[0]-1] * trun(0, 0, 0);
    store1(0, 0, 0) = store1(0, 0, 0) + store2(0, 0, 0) * fornow;

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[535].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 6);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[535].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[535].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[535].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[535].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[535].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
    }
}

#ifdef OPS_LAZY
extern "C"
void radcal_kernel_addspecies_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("radcal_kernel_addspecies", args, 6, 535, dim, 1, range, block, radcal_kernel_addspecies_execute);
}
#endif

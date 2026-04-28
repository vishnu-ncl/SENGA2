// Auto-generated at 2026-04-28 18:44:59.416564 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void temper_kernel_eqD_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void temper_kernel_eqD_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[9];
    args[0] = desc->args[0];
    args[1] = desc->args[1];
    args[2] = desc->args[2];
    args[3] = desc->args[3];
    args[4] = desc->args[4];
    args[5] = desc->args[5];
    args[6] = desc->args[6];
    args[7] = desc->args[7];
    args[8] = desc->args[8];
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

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 9, range, 562)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 562, "temper_kernel_eqD");
        block->instance->OPS_kernels[562].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "temper_kernel_eqD");
#endif

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
    int xdim0_temper_kernel_eqD = args[0].dat->size[0];
    int ydim0_temper_kernel_eqD = args[0].dat->size[1];
    int xdim1_temper_kernel_eqD = args[1].dat->size[0];
    int ydim1_temper_kernel_eqD = args[1].dat->size[1];
    int zdim1_temper_kernel_eqD = args[1].dat->size[2];
    int xdim2_temper_kernel_eqD = args[2].dat->size[0];
    int ydim2_temper_kernel_eqD = args[2].dat->size[1];
    int zdim2_temper_kernel_eqD = args[2].dat->size[2];
    int xdim3_temper_kernel_eqD = args[3].dat->size[0];
    int ydim3_temper_kernel_eqD = args[3].dat->size[1];
    int xdim4_temper_kernel_eqD = args[4].dat->size[0];
    int ydim4_temper_kernel_eqD = args[4].dat->size[1];
    int xdim5_temper_kernel_eqD = args[5].dat->size[0];
    int ydim5_temper_kernel_eqD = args[5].dat->size[1];
    int xdim6_temper_kernel_eqD = args[6].dat->size[0];
    int ydim6_temper_kernel_eqD = args[6].dat->size[1];
    int xdim7_temper_kernel_eqD = args[7].dat->size[0];
    int ydim7_temper_kernel_eqD = args[7].dat->size[1];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ trun_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int multi_d1 = getDatDimFromOpsArg(&args[1]);
    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, multi_d1);
    double * __restrict__ tcoeff_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int multi_d2 = getDatDimFromOpsArg(&args[2]);
    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, multi_d2);
    double * __restrict__ tderiv_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ urhs_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ vrhs_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style

    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, 1);
    double * __restrict__ wrhs_p = (double *)(args[6].data_d) + base6 - 1; // Subtracting 1 to convert to C-style

    int base7 = getDatBaseFromOpsArg3D(&args[7], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[7].data_d) + base7 - 1; // Subtracting 1 to convert to C-style

//  Subtracting 1 here as start_indx and end_indx is in Fortran style - converting it to c-style range
    int start_0 = start_indx[0]-1;
    int end_0 = end_indx[0];
    int arg_idx_0 = arg_idx[0];
    int start_1 = start_indx[1]-1;
    int end_1 = end_indx[1];
    int arg_idx_1 = arg_idx[1];
    int start_2 = start_indx[2]-1;
    int end_2 = end_indx[2];
    int arg_idx_2 = arg_idx[2];

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 9);
    ops_halo_exchanges(args, 9, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[562].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto nctmax_sycl = (*nctmax_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto nctmm1_sycl = (*nctmm1_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class temper_kernel_eqD_sycl>(
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

                int idx[] = {arg_idx_0+n_x, arg_idx_1+n_y, arg_idx_2+n_z};

                 ACC<double> trun(xdim0_temper_kernel_eqD, ydim0_temper_kernel_eqD, trun_p + (n_x * 1) + (n_y * xdim0_temper_kernel_eqD * 1) + (n_z * xdim0_temper_kernel_eqD * ydim0_temper_kernel_eqD * 1));
#ifdef OPS_SOA
                const  ACC<double> tcoeff(6,xdim1_temper_kernel_eqD, ydim1_temper_kernel_eqD, zdim1_temper_kernel_eqD, tcoeff_p + (n_x * 1) + (n_y * xdim1_temper_kernel_eqD * 1) + (n_z * xdim1_temper_kernel_eqD * ydim1_temper_kernel_eqD * 1));
#else
                const  ACC<double> tcoeff(6,xdim1_temper_kernel_eqD, ydim1_temper_kernel_eqD, zdim1_temper_kernel_eqD, tcoeff_p + 6 * ((n_x * 1) + (n_y * xdim1_temper_kernel_eqD * 1) + (n_z * xdim1_temper_kernel_eqD * ydim1_temper_kernel_eqD * 1)));
#endif
#ifdef OPS_SOA
                const  ACC<double> tderiv(5,xdim2_temper_kernel_eqD, ydim2_temper_kernel_eqD, zdim2_temper_kernel_eqD, tderiv_p + (n_x * 1) + (n_y * xdim2_temper_kernel_eqD * 1) + (n_z * xdim2_temper_kernel_eqD * ydim2_temper_kernel_eqD * 1));
#else
                const  ACC<double> tderiv(5,xdim2_temper_kernel_eqD, ydim2_temper_kernel_eqD, zdim2_temper_kernel_eqD, tderiv_p + 5 * ((n_x * 1) + (n_y * xdim2_temper_kernel_eqD * 1) + (n_z * xdim2_temper_kernel_eqD * ydim2_temper_kernel_eqD * 1)));
#endif
                const  ACC<double> drhs(xdim3_temper_kernel_eqD, ydim3_temper_kernel_eqD, drhs_p + (n_x * 1) + (n_y * xdim3_temper_kernel_eqD * 1) + (n_z * xdim3_temper_kernel_eqD * ydim3_temper_kernel_eqD * 1));
                const  ACC<double> urhs(xdim4_temper_kernel_eqD, ydim4_temper_kernel_eqD, urhs_p + (n_x * 1) + (n_y * xdim4_temper_kernel_eqD * 1) + (n_z * xdim4_temper_kernel_eqD * ydim4_temper_kernel_eqD * 1));
                const  ACC<double> vrhs(xdim5_temper_kernel_eqD, ydim5_temper_kernel_eqD, vrhs_p + (n_x * 1) + (n_y * xdim5_temper_kernel_eqD * 1) + (n_z * xdim5_temper_kernel_eqD * ydim5_temper_kernel_eqD * 1));
                const  ACC<double> wrhs(xdim6_temper_kernel_eqD, ydim6_temper_kernel_eqD, wrhs_p + (n_x * 1) + (n_y * xdim6_temper_kernel_eqD * 1) + (n_z * xdim6_temper_kernel_eqD * ydim6_temper_kernel_eqD * 1));
                const  ACC<double> erhs(xdim7_temper_kernel_eqD, ydim7_temper_kernel_eqD, erhs_p + (n_x * 1) + (n_y * xdim7_temper_kernel_eqD * 1) + (n_z * xdim7_temper_kernel_eqD * ydim7_temper_kernel_eqD * 1));

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double tempor;
    double tfpoly;
    double tdpoly;
    double deltmp;
    int ititrs;
    int icp;
    int ic;
    int jc;
    int kc;

double toltmp = 1.0e-10f;
int ntitrs = 100;

    ic = idx[0];
    jc = idx[1];
    kc = idx[2];
    tempor = trun(0, 0, 0);
    ititrs = 1;
    tfpoly = tcoeff(nctmax_sycl[0] + 1-1, 0, 0, 0);
    tdpoly = tderiv(nctmax_sycl[0]-1, 0, 0, 0);
    for (icp = nctmm1_sycl[0]; icp >= 1; icp -= 1) {
        tfpoly = tcoeff(icp + 1-1, 0, 0, 0) + tfpoly * tempor;
        tdpoly = tderiv(icp-1, 0, 0, 0) + tdpoly * tempor;
    }
    tfpoly = tcoeff(1-1, 0, 0, 0) + tfpoly * tempor;
    deltmp = -tfpoly / tdpoly;
    while (f2c::abs(deltmp) > toltmp) {
        if (ititrs < ntitrs) {
            tempor = tempor + deltmp;
            ititrs = ititrs + 1;
            tfpoly = tcoeff(nctmax_sycl[0] + 1-1, 0, 0, 0);
            tdpoly = tderiv(nctmax_sycl[0]-1, 0, 0, 0);
            for (icp = nctmm1_sycl[0]; icp >= 1; icp -= 1) {
                tfpoly = tcoeff(icp + 1-1, 0, 0, 0) + tfpoly * tempor;
                tdpoly = tderiv(icp-1, 0, 0, 0) + tdpoly * tempor;
            }
            tfpoly = tcoeff(1-1, 0, 0, 0) + tfpoly * tempor;
            deltmp = -tfpoly / tdpoly;
        } else {
            //std::cout << "'Fatal: TEMPER: T iteration failed to converge'" << std::endl;
            //std::cout << "'at point:'" << " " << ic << " " << jc << " " << kc << std::endl;
            //std::cout << "'with values:'" << " " << tempor << " " << deltmp << std::endl;
            //std::cout << drhs(0, 0, 0) << std::endl;
            //std::cout << urhs(0, 0, 0) << std::endl;
            //std::cout << vrhs(0, 0, 0) << std::endl;
            //std::cout << wrhs(0, 0, 0) << std::endl;
            //std::cout << erhs(0, 0, 0) << std::endl;
            f2c::trap();
        }
    }
    trun(0, 0, 0) = tempor;

                }

            });
        });
    }

    if (block->instance->OPS_diags > 1) {
        block->instance->sycl_instance->queue->wait();
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[562].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 9);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[562].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[562].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[562].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[562].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[562].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[562].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[562].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[562].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
        block->instance->OPS_kernels[562].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg7);
    }
}

#ifdef OPS_LAZY
extern "C"
void temper_kernel_eqD_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("temper_kernel_eqD", args, 9, 562, dim, 1, range, block, temper_kernel_eqD_execute);
}
#endif

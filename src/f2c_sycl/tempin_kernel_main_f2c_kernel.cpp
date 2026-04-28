// Auto-generated at 2026-04-28 18:44:59.774019 by ops-translator


//  ==================
//  Host stub function
//  ==================
#ifndef OPS_LAZY
extern "C"
void tempin_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
)
{
#else
void tempin_kernel_main_execute(ops_kernel_descriptor *desc)
{
    ops_block block = desc->block;
    int dim = desc->dim;
    int *range = desc->range;
    ops_arg args[17];
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
    args[16] = desc->args[16];
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
    ops_arg arg16 = args[16];

#if defined(CHECKPOINTING) && !defined(OPS_LAZY)
    if (!ops_checkpointing_before(args, 17, range, 565)) return;
#endif

    if (block->instance->OPS_diags > 1)
    {
        ops_timing_realloc(block->instance, 565, "tempin_kernel_main");
        block->instance->OPS_kernels[565].count++;
        ops_timers_core(&__c1, &__t1);
    }

#ifdef OPS_DEBUG
    ops_register_args(block->instance, args, "tempin_kernel_main");
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
    int xdim0_tempin_kernel_main = args[0].dat->size[0];
    int ydim0_tempin_kernel_main = args[0].dat->size[1];
    int xdim1_tempin_kernel_main = args[1].dat->size[0];
    int ydim1_tempin_kernel_main = args[1].dat->size[1];
    int xdim2_tempin_kernel_main = args[2].dat->size[0];
    int ydim2_tempin_kernel_main = args[2].dat->size[1];
    int xdim3_tempin_kernel_main = args[3].dat->size[0];
    int ydim3_tempin_kernel_main = args[3].dat->size[1];
    int xdim4_tempin_kernel_main = args[4].dat->size[0];
    int ydim4_tempin_kernel_main = args[4].dat->size[1];
    int xdim5_tempin_kernel_main = args[5].dat->size[0];
    int ydim5_tempin_kernel_main = args[5].dat->size[1];
    int xdim6_tempin_kernel_main = args[6].dat->size[0];
    int ydim6_tempin_kernel_main = args[6].dat->size[1];
    int zdim6_tempin_kernel_main = args[6].dat->size[2];

//  =======================================================
//  Set up initial pointers and exchange halos if necessary
//  =======================================================
    int base0 = getDatBaseFromOpsArg3D(&args[0], start_indx, 1);
    double * __restrict__ trun_p = (double *)(args[0].data_d) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[1].data_d) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ urhs_p = (double *)(args[2].data_d) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ vrhs_p = (double *)(args[3].data_d) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ wrhs_p = (double *)(args[4].data_d) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[5].data_d) + base5 - 1; // Subtracting 1 to convert to C-style

    int multi_d6 = getDatDimFromOpsArg(&args[6]);
    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, multi_d6);
    double * __restrict__ yrhs_p = (double *)(args[6].data_d) + base6 - 1; // Subtracting 1 to convert to C-style

    double *arg7h = (double *)args[7].data;

    double *arg8h = (double *)args[8].data;

    int *arg9h = (int *)args[9].data;

    int *arg10h = (int *)args[10].data;

    double *arg11h = (double *)args[11].data;

    int *arg12h = (int *)args[12].data;

    double trin_val = *(double *)args[13].data;

    int nspec_val = *(int *)args[14].data;

    int iproc_val = *(int *)args[15].data;

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

    int consts_bytes = 0;

    consts_bytes += ROUND_UP(arg7.dim*sizeof(double));

    consts_bytes += ROUND_UP(arg8.dim*sizeof(double));

    consts_bytes += ROUND_UP(arg9.dim*sizeof(int));

    consts_bytes += ROUND_UP(arg10.dim*sizeof(int));

    consts_bytes += ROUND_UP(arg11.dim*sizeof(double));

    consts_bytes += ROUND_UP(arg12.dim*sizeof(int));

    reallocConstArrays(block->instance, consts_bytes);
    consts_bytes = 0;

    arg7.data = block->instance->OPS_consts_h + consts_bytes;
    double* arg7_data_d = (double*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg7.dim; d++)     ((double *)arg7.data)[d] = arg7h[d];
    consts_bytes += ROUND_UP(arg7.dim*sizeof(double));
    arg8.data = block->instance->OPS_consts_h + consts_bytes;
    double* arg8_data_d = (double*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg8.dim; d++)     ((double *)arg8.data)[d] = arg8h[d];
    consts_bytes += ROUND_UP(arg8.dim*sizeof(double));
    arg9.data = block->instance->OPS_consts_h + consts_bytes;
    int* arg9_data_d = (int*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg9.dim; d++)     ((int *)arg9.data)[d] = arg9h[d];
    consts_bytes += ROUND_UP(arg9.dim*sizeof(int));
    arg10.data = block->instance->OPS_consts_h + consts_bytes;
    int* arg10_data_d = (int*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg10.dim; d++)     ((int *)arg10.data)[d] = arg10h[d];
    consts_bytes += ROUND_UP(arg10.dim*sizeof(int));
    arg11.data = block->instance->OPS_consts_h + consts_bytes;
    double* arg11_data_d = (double*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg11.dim; d++)     ((double *)arg11.data)[d] = arg11h[d];
    consts_bytes += ROUND_UP(arg11.dim*sizeof(double));
    arg12.data = block->instance->OPS_consts_h + consts_bytes;
    int* arg12_data_d = (int*)(block->instance->OPS_consts_d + consts_bytes);
    for (int d = 0; d < arg12.dim; d++)     ((int *)arg12.data)[d] = arg12h[d];
    consts_bytes += ROUND_UP(arg12.dim*sizeof(int));

    mvConstArraysToDevice(block->instance, consts_bytes);

//  =============
//  Halo exchange
//  =============
#ifndef OPS_LAZY
    ops_H_D_exchanges_device(args, 17);
    ops_halo_exchanges(args, 17, range);
#endif

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[565].mpi_time += __t2 - __t1;
    }

    if ((end_0-start_0)>0 && (end_1-start_1)>0 && (end_2-start_2)>0) {
        block->instance->sycl_instance->queue->submit([&](cl::sycl::handler &cgh) {

            auto ncofmx_sycl = (*ncofmx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto ntinmx_sycl = (*ntinmx_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto nctmax_sycl = (*nctmax_p).template get_access<cl::sycl::access::mode::read>(cgh);
            auto nctmm1_sycl = (*nctmm1_p).template get_access<cl::sycl::access::mode::read>(cgh);

            cgh.parallel_for<class tempin_kernel_main_sycl>(
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

                 ACC<double> trun(xdim0_tempin_kernel_main, ydim0_tempin_kernel_main, trun_p + (n_x * 1) + (n_y * xdim0_tempin_kernel_main * 1) + (n_z * xdim0_tempin_kernel_main * ydim0_tempin_kernel_main * 1));
                const  ACC<double> drhs(xdim1_tempin_kernel_main, ydim1_tempin_kernel_main, drhs_p + (n_x * 1) + (n_y * xdim1_tempin_kernel_main * 1) + (n_z * xdim1_tempin_kernel_main * ydim1_tempin_kernel_main * 1));
                const  ACC<double> urhs(xdim2_tempin_kernel_main, ydim2_tempin_kernel_main, urhs_p + (n_x * 1) + (n_y * xdim2_tempin_kernel_main * 1) + (n_z * xdim2_tempin_kernel_main * ydim2_tempin_kernel_main * 1));
                const  ACC<double> vrhs(xdim3_tempin_kernel_main, ydim3_tempin_kernel_main, vrhs_p + (n_x * 1) + (n_y * xdim3_tempin_kernel_main * 1) + (n_z * xdim3_tempin_kernel_main * ydim3_tempin_kernel_main * 1));
                const  ACC<double> wrhs(xdim4_tempin_kernel_main, ydim4_tempin_kernel_main, wrhs_p + (n_x * 1) + (n_y * xdim4_tempin_kernel_main * 1) + (n_z * xdim4_tempin_kernel_main * ydim4_tempin_kernel_main * 1));
                const  ACC<double> erhs(xdim5_tempin_kernel_main, ydim5_tempin_kernel_main, erhs_p + (n_x * 1) + (n_y * xdim5_tempin_kernel_main * 1) + (n_z * xdim5_tempin_kernel_main * ydim5_tempin_kernel_main * 1));
#ifdef OPS_SOA
                const  ACC<double> yrhs(22,xdim6_tempin_kernel_main, ydim6_tempin_kernel_main, zdim6_tempin_kernel_main, yrhs_p + (n_x * 1) + (n_y * xdim6_tempin_kernel_main * 1) + (n_z * xdim6_tempin_kernel_main * ydim6_tempin_kernel_main * 1));
#else
                const  ACC<double> yrhs(22,xdim6_tempin_kernel_main, ydim6_tempin_kernel_main, zdim6_tempin_kernel_main, yrhs_p + 22 * ((n_x * 1) + (n_y * xdim6_tempin_kernel_main * 1) + (n_z * xdim6_tempin_kernel_main * ydim6_tempin_kernel_main * 1)));
#endif

                const double *amascp = arg7_data_d;
                const double *amasct = arg8_data_d;
                const int *ncpoly = arg9_data_d;
                const int *ncenth = arg10_data_d;
                const double *tinthi = arg11_data_d;
                const int *ntint = arg12_data_d;
                const double *trin = &trin_val;
                const int *nspec = &nspec_val;
                const int *iproc = &iproc_val;

// =========
// User code
// =========
                if ((n_x < end_0-start_0) && (n_y < end_1-start_1) && (n_z < end_2-start_2)) {

    double tcoeff[(1 + 5 - (0))];
    double ukuk;
    double tempor;
    double tupper;
    double tlower;
    double tresid;
    double tuk2me;
    double cpfory;
    int ic;
    int jc;
    int kc;
    int ispec;
    int itint;
    int icp;
    int iindex;
    int ipower;
    bool fnconv;

double toltmp = 0.00010;
double tininc = 50.0;
double tlimlo = 200.0;
double tlimhi = 3000.0;

    ic = idx[0];
    jc = idx[1];
    kc = idx[2];
    ukuk = (urhs(0, 0, 0) * urhs(0, 0, 0) + vrhs(0, 0, 0) * vrhs(0, 0, 0) + wrhs(0, 0, 0) * wrhs(0, 0, 0)) / drhs(0, 0, 0);
    tuk2me = 0.5 * ukuk - erhs(0, 0, 0);
    tlower = tlimlo;
    tupper = tlimhi;
    fnconv = true;
    tempor = trin[0];
    tcoeff[0] = tuk2me;
    for (icp = 1; icp <= nctmax_sycl[0]; ++icp) {
        tcoeff[icp] = 0.0;
    }
    for (ispec = 1; ispec <= nspec[0]; ++ispec) {
        itint = 1;
        while (tempor > tinthi[(itint-1)+(ispec-1)*ntinmx_sycl[0]] && itint < ntint[ispec-1]) {
            itint = itint + 1;
        }
        tcoeff[0] = tcoeff[0] + yrhs(ispec-1, 0, 0, 0) * amascp[(ncenth[(itint-1)+(ispec-1)*ntinmx_sycl[0]]-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
        tcoeff[1] = tcoeff[1] + yrhs(ispec-1, 0, 0, 0) * amasct[(1-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
        for (icp = 2; icp <= ncpoly[(itint-1)+(ispec-1)*ntinmx_sycl[0]]; ++icp) {
            tcoeff[icp] = tcoeff[icp] + yrhs(ispec-1, 0, 0, 0) * amasct[(icp-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
        }
    }
    tresid = tcoeff[nctmax_sycl[0]];
    for (icp = nctmm1_sycl[0]; icp >= 1; icp -= 1) {
        tresid = tcoeff[icp] + tresid * tempor;
    }
    tresid = tcoeff[0] + tresid * tempor;
    if (f2c::abs(tresid) < toltmp) {
        fnconv = false;
    } else if (tresid < 0.0) {
        while (true) {
            tlower = tempor;
            tempor = tempor + tininc;
            tcoeff[0] = tuk2me;
            for (icp = 1; icp <= nctmax_sycl[0]; ++icp) {
                tcoeff[icp] = 0.0;
            }
            for (ispec = 1; ispec <= nspec[0]; ++ispec) {
                itint = 1;
                while (tempor > tinthi[(itint-1)+(ispec-1)*ntinmx_sycl[0]] && itint < ntint[ispec-1]) {
                    itint = itint + 1;
                }
                tcoeff[0] = tcoeff[0] + yrhs(ispec-1, 0, 0, 0) * amascp[(ncenth[(itint-1)+(ispec-1)*ntinmx_sycl[0]]-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
                tcoeff[1] = tcoeff[1] + yrhs(ispec-1, 0, 0, 0) * amasct[(1-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
                for (icp = 2; icp <= ncpoly[(itint-1)+(ispec-1)*ntinmx_sycl[0]]; ++icp) {
                    tcoeff[icp] = tcoeff[icp] + yrhs(ispec-1, 0, 0, 0) * amasct[(icp-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
                }
            }
            tresid = tcoeff[nctmax_sycl[0]];
            for (icp = nctmm1_sycl[0]; icp >= 1; icp -= 1) {
                tresid = tcoeff[icp] + tresid * tempor;
            }
            tresid = tcoeff[0] + tresid * tempor;
            if (f2c::abs(tresid) < toltmp) {
                fnconv = false;
                break;
            } else if (tresid < 0.0) {
                if (tempor >= tlimhi) {
                    //std::cout << "'Fatal: TEMPIN: T upper bracket failed to converge'" << std::endl;
                    //std::cout << "'processor:'" << " " << iproc[0] << std::endl;
                    //std::cout << "'at point:'" << " " << ic << " " << jc << " " << kc << std::endl;
                    //std::cout << "'with values:'" << " " << tempor << " " << tresid << std::endl;
                    //std::cout << drhs(0, 0, 0) << std::endl;
                    //std::cout << urhs(0, 0, 0) << std::endl;
                    //std::cout << vrhs(0, 0, 0) << std::endl;
                    //std::cout << wrhs(0, 0, 0) << std::endl;
                    //std::cout << erhs(0, 0, 0) << std::endl;
                    for (ispec = 1; ispec <= nspec[0]; ++ispec) {
                        //std::cout << yrhs(ispec-1, 0, 0, 0) << std::endl;
                    }
                    f2c::trap();
                }
            } else if (tresid > 0.0) {
                tupper = tempor;
                break;
            }
        }
    } else if (tresid > 0.0) {
        while (true) {
            tupper = tempor;
            tempor = tempor - tininc;
            tcoeff[0] = tuk2me;
            for (icp = 1; icp <= nctmax_sycl[0]; ++icp) {
                tcoeff[icp] = 0.0;
            }
            for (ispec = 1; ispec <= nspec[0]; ++ispec) {
                itint = 1;
                while (tempor > tinthi[(itint-1)+(ispec-1)*ntinmx_sycl[0]] && itint < ntint[ispec-1]) {
                    itint = itint + 1;
                }
                tcoeff[0] = tcoeff[0] + yrhs(ispec-1, 0, 0, 0) * amascp[(ncenth[(itint-1)+(ispec-1)*ntinmx_sycl[0]]-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
                tcoeff[1] = tcoeff[1] + yrhs(ispec-1, 0, 0, 0) * amasct[(1-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
                for (icp = 2; icp <= ncpoly[(itint-1)+(ispec-1)*ntinmx_sycl[0]]; ++icp) {
                    tcoeff[icp] = tcoeff[icp] + yrhs(ispec-1, 0, 0, 0) * amasct[(icp-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
                }
            }
            tresid = tcoeff[nctmax_sycl[0]];
            for (icp = nctmm1_sycl[0]; icp >= 1; icp -= 1) {
                tresid = tcoeff[icp] + tresid * tempor;
            }
            tresid = tcoeff[0] + tresid * tempor;
            if (f2c::abs(tresid) < toltmp) {
                fnconv = false;
                break;
            } else if (tresid > 0.0) {
                if (tempor <= tlimlo) {
                    //std::cout << "'Fatal: TEMPIN: T lower bracket failed to converge'" << std::endl;
                    //std::cout << "'processor:'" << " " << iproc[0] << std::endl;
                    //std::cout << "'at point:'" << " " << ic << " " << jc << " " << kc << std::endl;
                    //std::cout << "'with values:'" << " " << tempor << " " << tresid << std::endl;
                    //std::cout << drhs(0, 0, 0) << std::endl;
                    //std::cout << urhs(0, 0, 0) << std::endl;
                    //std::cout << vrhs(0, 0, 0) << std::endl;
                    //std::cout << wrhs(0, 0, 0) << std::endl;
                    //std::cout << erhs(0, 0, 0) << std::endl;
                    for (ispec = 1; ispec <= nspec[0]; ++ispec) {
                        //std::cout << yrhs(ispec-1, 0, 0, 0) << std::endl;
                    }
                    f2c::trap();
                }
            } else if (tresid < 0.0) {
                tlower = tempor;
                break;
            }
        }
    }
    if (fnconv) {
        while (true) {
            tempor = 0.5 * (tlower + tupper);
            tcoeff[0] = tuk2me;
            for (icp = 1; icp <= nctmax_sycl[0]; ++icp) {
                tcoeff[icp] = 0.0;
            }
            for (ispec = 1; ispec <= nspec[0]; ++ispec) {
                itint = 1;
                while (tempor > tinthi[(itint-1)+(ispec-1)*ntinmx_sycl[0]] && itint < ntint[ispec-1]) {
                    itint = itint + 1;
                }
                tcoeff[0] = tcoeff[0] + yrhs(ispec-1, 0, 0, 0) * amascp[(ncenth[(itint-1)+(ispec-1)*ntinmx_sycl[0]]-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
                tcoeff[1] = tcoeff[1] + yrhs(ispec-1, 0, 0, 0) * amasct[(1-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
                for (icp = 2; icp <= ncpoly[(itint-1)+(ispec-1)*ntinmx_sycl[0]]; ++icp) {
                    tcoeff[icp] = tcoeff[icp] + yrhs(ispec-1, 0, 0, 0) * amasct[(icp-1)+(itint-1)*ncofmx_sycl[0]+(ispec-1)*ncofmx_sycl[0]*ntinmx_sycl[0]];
                }
            }
            tresid = tcoeff[nctmax_sycl[0]];
            for (icp = nctmm1_sycl[0]; icp >= 1; icp -= 1) {
                tresid = tcoeff[icp] + tresid * tempor;
            }
            tresid = tcoeff[0] + tresid * tempor;
            if (f2c::abs(tresid) < toltmp) {
                trun(0, 0, 0) = tempor;
                break;
            } else if (tresid < 0.0) {
                tlower = tempor;
            } else if (tresid > 0.0) {
                tupper = tempor;
            }
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
        block->instance->OPS_kernels[565].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_device(args, 17);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

    if (block->instance->OPS_diags > 1) {
//      ====================
//      Update kernel record
//      ====================
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[565].mpi_time += __t2 -__t1;
        block->instance->OPS_kernels[565].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg0);
        block->instance->OPS_kernels[565].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg1);
        block->instance->OPS_kernels[565].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg2);
        block->instance->OPS_kernels[565].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg3);
        block->instance->OPS_kernels[565].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg4);
        block->instance->OPS_kernels[565].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg5);
        block->instance->OPS_kernels[565].transfer += ops_compute_transfer(dim, start_indx, end_indx, &arg6);
    }
}

#ifdef OPS_LAZY
extern "C"
void tempin_kernel_main_host_c(
    const char *name,
    ops_arg *args,
    int nargs,
    int dim,
    int *range,
    ops_block block
    )
{

    create_kerneldesc_and_enque("tempin_kernel_main", args, 17, 565, dim, 1, range, block, tempin_kernel_main_execute);
}
#endif

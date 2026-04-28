// Auto-generated at 2026-04-28 18:43:38.321941 by ops-translator


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
    double * __restrict__ trun_p = (double *)(args[0].data) + base0 - 1; // Subtracting 1 to convert to C-style

    int base1 = getDatBaseFromOpsArg3D(&args[1], start_indx, 1);
    double * __restrict__ drhs_p = (double *)(args[1].data) + base1 - 1; // Subtracting 1 to convert to C-style

    int base2 = getDatBaseFromOpsArg3D(&args[2], start_indx, 1);
    double * __restrict__ urhs_p = (double *)(args[2].data) + base2 - 1; // Subtracting 1 to convert to C-style

    int base3 = getDatBaseFromOpsArg3D(&args[3], start_indx, 1);
    double * __restrict__ vrhs_p = (double *)(args[3].data) + base3 - 1; // Subtracting 1 to convert to C-style

    int base4 = getDatBaseFromOpsArg3D(&args[4], start_indx, 1);
    double * __restrict__ wrhs_p = (double *)(args[4].data) + base4 - 1; // Subtracting 1 to convert to C-style

    int base5 = getDatBaseFromOpsArg3D(&args[5], start_indx, 1);
    double * __restrict__ erhs_p = (double *)(args[5].data) + base5 - 1; // Subtracting 1 to convert to C-style

    int multi_d6 = getDatDimFromOpsArg(&args[6]);
    int base6 = getDatBaseFromOpsArg3D(&args[6], start_indx, multi_d6);
    double * __restrict__ yrhs_p = (double *)(args[6].data) + base6 - 1; // Subtracting 1 to convert to C-style

    double * __restrict__  amascp = (double *)args[7].data;
    double * __restrict__  amasct = (double *)args[8].data;
    int * __restrict__  ncpoly = (int *)args[9].data;
    int * __restrict__  ncenth = (int *)args[10].data;
    double * __restrict__  tinthi = (double *)args[11].data;
    int * __restrict__  ntint = (int *)args[12].data;
    double * __restrict__  trin = (double *)args[13].data;
    int * __restrict__  nspec = (int *)args[14].data;
    int * __restrict__  iproc = (int *)args[15].data;

//  ==============
//  Halo exchanges
//  ==============
#ifndef OPS_LAZY
    ops_H_D_exchanges_host(args, 17);
    ops_halo_exchanges(args, 17, range);
    ops_H_D_exchanges_host(args, 17);
#endif //OPS_LAZY

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c2, &__t2);
        block->instance->OPS_kernels[565].mpi_time += __t2 - __t1;
    }

    for (int n_z = 0; n_z < end_indx[2]-start_indx[2] +1; n_z++)
    {
        for (int n_y = 0; n_y < end_indx[1]-start_indx[1] +1; n_y++)
        {
            for(int n_x = 0; n_x < end_indx[0]-start_indx[0] +1; n_x++)
            {
                int idx[] = {arg_idx[0] + n_x, arg_idx[1] + n_y, arg_idx[2] + n_z};

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
    for (icp = 1; icp <= nctmax; ++icp) {
        tcoeff[icp] = 0.0;
    }
    for (ispec = 1; ispec <= nspec[0]; ++ispec) {
        itint = 1;
        while (tempor > tinthi[(itint-1)+(ispec-1)*ntinmx] && itint < ntint[ispec-1]) {
            itint = itint + 1;
        }
        tcoeff[0] = tcoeff[0] + yrhs(ispec-1, 0, 0, 0) * amascp[(ncenth[(itint-1)+(ispec-1)*ntinmx]-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
        tcoeff[1] = tcoeff[1] + yrhs(ispec-1, 0, 0, 0) * amasct[(1-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
        for (icp = 2; icp <= ncpoly[(itint-1)+(ispec-1)*ntinmx]; ++icp) {
            tcoeff[icp] = tcoeff[icp] + yrhs(ispec-1, 0, 0, 0) * amasct[(icp-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
        }
    }
    tresid = tcoeff[nctmax];
    for (icp = nctmm1; icp >= 1; icp -= 1) {
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
            for (icp = 1; icp <= nctmax; ++icp) {
                tcoeff[icp] = 0.0;
            }
            for (ispec = 1; ispec <= nspec[0]; ++ispec) {
                itint = 1;
                while (tempor > tinthi[(itint-1)+(ispec-1)*ntinmx] && itint < ntint[ispec-1]) {
                    itint = itint + 1;
                }
                tcoeff[0] = tcoeff[0] + yrhs(ispec-1, 0, 0, 0) * amascp[(ncenth[(itint-1)+(ispec-1)*ntinmx]-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
                tcoeff[1] = tcoeff[1] + yrhs(ispec-1, 0, 0, 0) * amasct[(1-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
                for (icp = 2; icp <= ncpoly[(itint-1)+(ispec-1)*ntinmx]; ++icp) {
                    tcoeff[icp] = tcoeff[icp] + yrhs(ispec-1, 0, 0, 0) * amasct[(icp-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
                }
            }
            tresid = tcoeff[nctmax];
            for (icp = nctmm1; icp >= 1; icp -= 1) {
                tresid = tcoeff[icp] + tresid * tempor;
            }
            tresid = tcoeff[0] + tresid * tempor;
            if (f2c::abs(tresid) < toltmp) {
                fnconv = false;
                break;
            } else if (tresid < 0.0) {
                if (tempor >= tlimhi) {
                    std::cout << "'Fatal: TEMPIN: T upper bracket failed to converge'" << std::endl;
                    std::cout << "'processor:'" << " " << iproc[0] << std::endl;
                    std::cout << "'at point:'" << " " << ic << " " << jc << " " << kc << std::endl;
                    std::cout << "'with values:'" << " " << tempor << " " << tresid << std::endl;
                    std::cout << drhs(0, 0, 0) << std::endl;
                    std::cout << urhs(0, 0, 0) << std::endl;
                    std::cout << vrhs(0, 0, 0) << std::endl;
                    std::cout << wrhs(0, 0, 0) << std::endl;
                    std::cout << erhs(0, 0, 0) << std::endl;
                    for (ispec = 1; ispec <= nspec[0]; ++ispec) {
                        std::cout << yrhs(ispec-1, 0, 0, 0) << std::endl;
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
            for (icp = 1; icp <= nctmax; ++icp) {
                tcoeff[icp] = 0.0;
            }
            for (ispec = 1; ispec <= nspec[0]; ++ispec) {
                itint = 1;
                while (tempor > tinthi[(itint-1)+(ispec-1)*ntinmx] && itint < ntint[ispec-1]) {
                    itint = itint + 1;
                }
                tcoeff[0] = tcoeff[0] + yrhs(ispec-1, 0, 0, 0) * amascp[(ncenth[(itint-1)+(ispec-1)*ntinmx]-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
                tcoeff[1] = tcoeff[1] + yrhs(ispec-1, 0, 0, 0) * amasct[(1-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
                for (icp = 2; icp <= ncpoly[(itint-1)+(ispec-1)*ntinmx]; ++icp) {
                    tcoeff[icp] = tcoeff[icp] + yrhs(ispec-1, 0, 0, 0) * amasct[(icp-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
                }
            }
            tresid = tcoeff[nctmax];
            for (icp = nctmm1; icp >= 1; icp -= 1) {
                tresid = tcoeff[icp] + tresid * tempor;
            }
            tresid = tcoeff[0] + tresid * tempor;
            if (f2c::abs(tresid) < toltmp) {
                fnconv = false;
                break;
            } else if (tresid > 0.0) {
                if (tempor <= tlimlo) {
                    std::cout << "'Fatal: TEMPIN: T lower bracket failed to converge'" << std::endl;
                    std::cout << "'processor:'" << " " << iproc[0] << std::endl;
                    std::cout << "'at point:'" << " " << ic << " " << jc << " " << kc << std::endl;
                    std::cout << "'with values:'" << " " << tempor << " " << tresid << std::endl;
                    std::cout << drhs(0, 0, 0) << std::endl;
                    std::cout << urhs(0, 0, 0) << std::endl;
                    std::cout << vrhs(0, 0, 0) << std::endl;
                    std::cout << wrhs(0, 0, 0) << std::endl;
                    std::cout << erhs(0, 0, 0) << std::endl;
                    for (ispec = 1; ispec <= nspec[0]; ++ispec) {
                        std::cout << yrhs(ispec-1, 0, 0, 0) << std::endl;
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
            for (icp = 1; icp <= nctmax; ++icp) {
                tcoeff[icp] = 0.0;
            }
            for (ispec = 1; ispec <= nspec[0]; ++ispec) {
                itint = 1;
                while (tempor > tinthi[(itint-1)+(ispec-1)*ntinmx] && itint < ntint[ispec-1]) {
                    itint = itint + 1;
                }
                tcoeff[0] = tcoeff[0] + yrhs(ispec-1, 0, 0, 0) * amascp[(ncenth[(itint-1)+(ispec-1)*ntinmx]-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
                tcoeff[1] = tcoeff[1] + yrhs(ispec-1, 0, 0, 0) * amasct[(1-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
                for (icp = 2; icp <= ncpoly[(itint-1)+(ispec-1)*ntinmx]; ++icp) {
                    tcoeff[icp] = tcoeff[icp] + yrhs(ispec-1, 0, 0, 0) * amasct[(icp-1)+(itint-1)*ncofmx+(ispec-1)*ncofmx*ntinmx];
                }
            }
            tresid = tcoeff[nctmax];
            for (icp = nctmm1; icp >= 1; icp -= 1) {
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
      }

    }

    if (block->instance->OPS_diags > 1) {
        ops_timers_core(&__c1, &__t1);
        block->instance->OPS_kernels[565].time += __t1 - __t2;
    }

#ifndef OPS_LAZY
    ops_set_dirtybit_host(args, 17);
    ops_set_halo_dirtybit3(&args[0], range);
#endif

//  ====================
//  Update kernel record
//  ====================
    if (block->instance->OPS_diags > 1) {
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

    create_kerneldesc_and_enque("tempin_kernel_main", args, 17, 565, dim, 0, range, block, tempin_kernel_main_execute);
}
#endif

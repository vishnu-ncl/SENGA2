#!/bin/bash

rm senga2

mpif90 -g -fcheck=all -cpp -O0 -Wall -mcmodel=medium -fp-model strict -fp-model source -prec-div -prec-sqrt -module /ext-home/asl/OPS/ops/fortran/mod/intel -I/ext-home/asl/OPS/ops/fortran/mod/intel  -L/ext-home/asl/OPS/ops/fortran/lib/intel -lstdc++ data_type.F90 constants.F90 com_espect.F90 com_senga2.F90 pardom.F90 parfer.F90 tempin.F90 flinit.F90 output.F90 dfmstr.F90 exhalo.F90 inhalo.F90 exhals.F90 inhals.F90 ardums.F90 ardump.F90 buftxl.F90 chkarr.F90 chkary.F90 cxhalo.F90 cxhals.F90 cyhalo.F90 cyhals.F90 czhalo.F90 czhals.F90 erfunc.F90 espect.F90 espksq.F90 espovk.F90 finish.F90 diffin.F90 hdf5io.F90 contin.F90 bcinit.F90 turbft.F90 fftsym.F90 fftf3d.F90 integf.F90 fftixl.F90 MPI/*_seq_kernel.F90 senga2_ops.F90 com_ops_senga2_ops.F90 ops_data_init_ops.F90 dfinit_ops.F90 dfbydx_ops.F90 d2fdx2_ops.F90 dfbydy_ops.F90 dfbydz_ops.F90 d2fdy2_ops.F90 d2fdz2_ops.F90 zeroxl_ops.F90 zeroxr_ops.F90 zeroyl_ops.F90 zeroyr_ops.F90 zerozl_ops.F90 zerozr_ops.F90 d2fdxy_ops.F90 d2fdxz_ops.F90 d2fdyz_ops.F90 rhscal_ops.F90 rhsvel_ops.F90 lincom_ops.F90 fincom_ops.F90 adaptt_ops.F90 bcdtxl_ops.F90 bcdtxr_ops.F90 bcdtyl_ops.F90 bcdtyr_ops.F90 bcdtzl_ops.F90 bcdtzr_ops.F90 bcttxl_ops.F90 bcttxr_ops.F90 bcttyl_ops.F90 bcttyr_ops.F90 bcttzl_ops.F90 bcttzr_ops.F90 bcutxl_ops.F90 bcutxr_ops.F90 bcutyl_ops.F90 bcutyr_ops.F90 bcutzl_ops.F90 bcutzr_ops.F90 bcytxl_ops.F90 bcytxr_ops.F90 bcytyl_ops.F90 bcytyr_ops.F90 bcytzl_ops.F90 bcytzr_ops.F90 bounds_ops.F90 bountt_ops.F90 boundt_ops.F90 radcal_ops.F90 radcin_ops.F90 bctixl_ops.F90 flamin_ops.F90 chrate_ops.F90 indata_ops.F90 dtinit_ops.F90 chemin_ops.F90 espini_ops.F90 turbin_ops.F90 temper_ops.F90 \
    -o senga2 fftlib/fftlib.a rndlib/rndlib.a mpi_stub/mpi_stub.a -lm -lops_for_seq

# parlib/parlib.a mpi_stub/mpi_stub.a
#rm senga2_mpi

#mpif90 -g -fcheck=all -cpp -O0 -Wall -mcmodel=medium -fp-model strict -fp-model source -prec-div -prec-sqrt -module /ext-home/asl/OPS/ops/fortran/mod/intel -I/ext-home/asl/OPS/ops/fortran/mod/intel  -L/ext-home/asl/OPS/ops/fortran/lib/intel -lstdc++ -DOPS_MPI data_type.F90 constants.F90 com_espect.F90 com_senga2.F90 pardom.F90 parfer.F90 tempin.F90 flinit.F90 output.F90 dfmstr.F90 exhalo.F90 inhalo.F90 exhals.F90 inhals.F90 ardums.F90 ardump.F90 buftxl.F90 chkarr.F90 chkary.F90 cxhalo.F90 cxhals.F90 cyhalo.F90 cyhals.F90 czhalo.F90 czhals.F90 erfunc.F90 espect.F90 espksq.F90 espovk.F90 finish.F90 diffin.F90 hdf5io.F90 contin.F90 bcinit.F90 turbft.F90 fftsym.F90 fftf3d.F90 integf.F90 fftixl.F90 MPI/*_seq_kernel.F90 senga2_ops.F90 com_ops_senga2_ops.F90 ops_data_init_ops.F90 dfinit_ops.F90 dfbydx_ops.F90 d2fdx2_ops.F90 dfbydy_ops.F90 dfbydz_ops.F90 d2fdy2_ops.F90 d2fdz2_ops.F90 zeroxl_ops.F90 zeroxr_ops.F90 zeroyl_ops.F90 zeroyr_ops.F90 zerozl_ops.F90 zerozr_ops.F90 d2fdxy_ops.F90 d2fdxz_ops.F90 d2fdyz_ops.F90 rhscal_ops.F90 rhsvel_ops.F90 lincom_ops.F90 fincom_ops.F90 adaptt_ops.F90 bcdtxl_ops.F90 bcdtxr_ops.F90 bcdtyl_ops.F90 bcdtyr_ops.F90 bcdtzl_ops.F90 bcdtzr_ops.F90 bcttxl_ops.F90 bcttxr_ops.F90 bcttyl_ops.F90 bcttyr_ops.F90 bcttzl_ops.F90 bcttzr_ops.F90 bcutxl_ops.F90 bcutxr_ops.F90 bcutyl_ops.F90 bcutyr_ops.F90 bcutzl_ops.F90 bcutzr_ops.F90 bcytxl_ops.F90 bcytxr_ops.F90 bcytyl_ops.F90 bcytyr_ops.F90 bcytzl_ops.F90 bcytzr_ops.F90 bounds_ops.F90 bountt_ops.F90 boundt_ops.F90 radcal_ops.F90 radcin_ops.F90 bctixl_ops.F90 flamin_ops.F90 chrate_ops.F90 indata_ops.F90 dtinit_ops.F90 chemin_ops.F90 espini_ops.F90 turbin_ops.F90 temper_ops.F90 \
#        -o senga2_mpi fftlib/fftlib.a rndlib/rndlib.a parlib/parlib.a -lm -lops_for_mpi








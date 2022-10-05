
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.common
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.mpi
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.cuda
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hdf5

OTHER_FILES = data_type.F90 constants.F90 com_espect.F90 com_senga2.F90 \
	pardom.F90 indata.F90 \
    dtinit.F90 parfer.F90 tempin.F90 flinit.F90 output.F90 \
    dfmstr.F90 \
    exhalo.F90 inhalo.F90 exhals.F90 inhals.F90 \
    adaptt.F90 ardums.F90 bcdtxl.F90 bcdtxr.F90 bcdtyl.F90 bcdtyr.F90 \
    bcdtzl.F90 bcdtzr.F90 bctixl.F90 bcttxl.F90 bcttxr.F90 bcttyl.F90 \
    bcttyr.F90 bcttzl.F90 bcttzr.F90 ardump.F90 \
    bcutxl.F90 bcutxr.F90 bcutyl.F90 bcutyr.F90 \
    bcutzl.F90 bcutzr.F90 bcytxl.F90 bcytxr.F90 bcytyl.F90 bcytyr.F90 \
    bcytzr.F90 buftxl.F90 chkarr.F90 chkary.F90 cxhalo.F90 cxhals.F90 \
    cyhalo.F90 cyhals.F90 czhalo.F90 czhals.F90 erfunc.F90 espect.F90 \
    espksq.F90 espovk.F90 finish.F90 bcytzl.F90 \
    temper.F90 chrate.F90 chemin.F90 diffin.F90 radcal.F90 radcin.F90 \
    hdf5io.F90 contin.F90 \
    bcinit.F90 bounds.F90 bountt.F90 boundt.F90 \
    turbin.F90 turbft.F90 fftsym.F90 fftf3d.F90 \
    espini.F90 integf.F90 \
    flamin.F90 \
    fftixl.F90

OPS_FILES = senga2.F90 com_ops_senga2.F90 ops_data_init.F90 dfinit.F90 dfbydx.F90 d2fdx2.F90 \
    dfbydy.F90 dfbydz.F90 d2fdy2.F90 d2fdz2.F90 zeroxl.F90 zeroxr.F90 zeroyl.F90 zeroyr.F90 \
    zerozl.F90 zerozr.F90 d2fdxy.F90 d2fdxz.F90 d2fdyz.F90 rhscal.F90 rhsvel.F90 \
    lincom.F90 fincom.F90

OPS_GENERATED = senga2_ops.F90 com_ops_senga2_ops.F90 ops_data_init_ops.F90 dfinit_ops.F90 dfbydx_ops.F90 d2fdx2_ops.F90 \
    dfbydy_ops.F90 dfbydz_ops.F90 d2fdy2_ops.F90 d2fdz2_ops.F90 zeroxl_ops.F90 zeroxr_ops.F90 zeroyl_ops.F90 zeroyr_ops.F90 \
    zerozl_ops.F90 zerozr_ops.F90 d2fdxy_ops.F90 d2fdxz_ops.F90 d2fdyz_ops.F90 rhscal_ops.F90 rhsvel_ops.F90 \
    lincom_ops.F90 fincom_ops.F90

APP_LIBS = fftlib/fftlib.a rndlib/rndlib.a parlib/parlib.a -lm

APP = senga2
MAIN_SRC=senga2

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.f_app


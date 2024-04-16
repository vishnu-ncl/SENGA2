include $(OPS_INSTALL_PATH)/../makefiles/Makefile.common
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.mpi
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.cuda
USE_HDF5=1
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hdf5

OTHER_FILES=constants.F90

OPS_FILES=verify.F90 

OPS_GENERATED_SUFFIX=_ops

APP = verify
MAIN_SRC=verify

TARGETS=seq

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.f_app


include $(OPS_INSTALL_PATH)/../makefiles/Makefile.common
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.mpi
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.cuda
USE_HDF5=1
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hdf5

OTHER_FILES=constants.F90

OPS_FILES=flamin.F90 

OPS_GENERATED_SUFFIX=_ops

APP = flamin
MAIN_SRC=flamin

TARGETS=seq

APP_FLAGS = -fPIC

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.f_app


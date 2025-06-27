# Makefile for building flamin.x using GNU Fortran and HDF5

FC = ftn
# Adjust as needed
HDF5_INSTALL_PATH=/opt/cray/pe/hdf5-parallel/1.12.2.11/crayclang/17.0/include
HDF5_INC = -I$(HDF5_INSTALL_PATH)/include
HDF5_LIBS = -L$(HDF5_INSTALL_PATH)/lib -lhdf5_hl -lhdf5 -lz
FFLAGS = -O2 -g -Wall $(HDF5_INC) -fallow-argument-mismatch 

# Source files
SRC = com_senga2.F90 write_ops_hdf5.F90 flamin.F90

# Output executable
TARGET = write_flame

# Build rule
$(TARGET): $(SRC)
	$(FC) $(FFLAGS) -o $(TARGET) $(SRC) $(HDF5_LIBS)

# Clean rule
clean:
	rm -f $(TARGET) *.o *.mod


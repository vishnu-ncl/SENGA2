#
# make file for SENGA2
#
# RSC 08-JAN-2005 LAM/MPI implementation
# RSC 12-JUN-2008 SUSE/MPICH implementation
#

# Fortran compiler
FC = mpif90 -g -fcheck=all -cpp
#FC = h5pfc -g -fcheck=all -DHDF5 -cpp

# Quick non-optimized compile
#FFLAGS = -g -fno-silent -ftrapping-math -ftrapv -Wimplicit -Wall

# Non-optimized compile with profiling
#FFLAGS = -pg -fno-silent -ftrapping-math -ftrapv -Wimplicit -Wall

# Optimized compile
FFLAGS = -O3 -Wall -mcmodel=medium

# Optimized compile with profiling
#FFLAGS = -pg -O3 -fno-silent -Wimplicit -Wall

# Source files
SRC =   rhscal.f lincom.f fincom.f \
	temper.f chrate.f chemin.f diffin.f radcal.f radcin.f \
	hdf5io.f contin.f finish.f \
	dfbydx.f d2fdx2.f \
	adaptt.f \
	turbin.f turbft.f fftsym.f fftf3d.f \
	espect.f espini.f espovk.f espksq.f integf.f \
	flamin.f erfunc.f \
	bcinit.f bounds.f boundt.f bountt.f \
	bcutxl.f bcttxl.f bcdtxl.f bcytxl.f bctixl.f buftxl.f fftixl.f \
	bcutxr.f bcttxr.f bcdtxr.f bcytxr.f \
	bcutyl.f bcttyl.f bcdtyl.f bcytyl.f \
	bcutyr.f bcttyr.f bcdtyr.f bcytyr.f \
	bcutzl.f bcttzl.f bcdtzl.f bcytzl.f \
	bcutzr.f bcttzr.f bcdtzr.f bcytzr.f \
	cxhalo.f cyhalo.f czhalo.f cxhals.f cyhals.f czhals.f \
	ardump.f ardums.f chkarr.f chkary.f

SRC_F90 = data_type.F90 com_espect.F90 com_senga2.F90 senga2.F90 \
	rhsvel.F90 dfinit.F90 dfbydy.F90 dfbydz.F90 pardom.F90 indata.F90 \
    dtinit.F90 parfer.F90 tempin.F90 flinit.F90 output.F90 \
    d2fdy2.F90 d2fdz2.F90 \
    d2fdxy.F90 d2fdxz.F90 d2fdyz.F90 dfmstr.F90 \
    exhalo.F90 inhalo.F90 exhals.F90 inhals.F90 \
    zeroxl.F90 zeroxr.F90 zeroyl.F90 zeroyr.F90 zerozl.F90 zerozr.F90 

# Libraries: need separate compilation in their own subdirectories
# FFTs; random number generators
# MPI libraries
LIBS = fftlib/fftlib.a rndlib/rndlib.a parlib/parlib.a -lm

# Object files are the same as the source files but with .o rather than .f
OBJS=$(SRC:.f=.o)
OBJS_F90=$(SRC_F90:.F90=.o)

# Default compilation of .f files to .o
.f.o: $(FC) $(INC) -c $<

%.o: %.F90
	$(FC) -c $<

# Full make of SENGA2
all: senga2

senga2: $(OBJS_F90) $(OBJS) 
	$(FC) $(OBJS) $(OBJS_F90) -o senga2 $(LIBS)

# Add -pg for profiling
#senga2: $(OBJS) 
#	$(FC) $(OBJS) -pg -o senga2 $(LIBS)

# Command to clean up after
clean: 
	rm -f *.o senga2 core *.mod

#
# make file for SENGA2
#
# RSC 08-JAN-2005 LAM/MPI implementation
# RSC 12-JUN-2008 SUSE/MPICH implementation
#

# Fortran compiler
#FC = ftn #-g -fcheck=all -cpp
#FC = h5pfc -g -fcheck=all -DHDF5 -cpp
FC = ftn

# GNU
#FFLAGS = -O3 -g -fopenmp -ftree-vectorize -ffree-form -ffree-line-length-none -fallow-argument-mismatch -mcmodel=medium

# GNU-ARCHER2
FFLAGS = -O3 -march=znver2 -mtune=znver2 -g -ftree-vectorize -fopenmp -ffree-form -ffree-line-length-none -fallow-argument-mismatch -mcmodel=medium

#FFLAGS = -O3 -march=znver2 -mtune=znver2 -g -ftree-vectorize -fopenmp -ffree-form -ffree-line-length-none -fallow-argument-mismatch -mcmodel=medium -fopt-info-vec -fopt-info-vec-optimized -fopt-info-vec-missed -funroll-loops -finline-functions -fstrict-aliasing

# INTEL
#FFLAGS = -O3 -g -mcmodel=large -fp-model precise -fp-model source -no-fma -fimf-arch-consistency=true -fimf-precision=high

# CRAY
#FFLAGS = -O3 -g -fopenmp -fPIC -M 878 -f free -N 1023 -h vector3 
#-h flex_mp=strict -h nofma -hfp1=noapprox -K trap=divz,fp,inv,unf

# Source files hdf5io.F90
SRC = 

SRC_F90 = com_espect.F90 com_senga2.F90 senga2.F90 \
	rhsvel.F90 dfinit.F90 dfbydy.F90 dfbydz.F90 pardom.F90 indata.F90 \
    dtinit.F90 parfer.F90 tempin.F90 flinit.F90 output.F90 \
    d2fdy2.F90 d2fdz2.F90 \
    d2fdxy.F90 d2fdxz.F90 d2fdyz.F90 dfmstr.F90 \
    exhalo.F90 inhalo.F90 exhals.F90 inhals.F90 \
    zeroxl.F90 zeroxr.F90 zeroyl.F90 zeroyr.F90 zerozl.F90 zerozr.F90 \
    adaptt.F90 ardums.F90 bcdtxl.F90 bcdtxr.F90 bcdtyl.F90 bcdtyr.F90 \
    bcdtzl.F90 bcdtzr.F90 bctixl.F90 bcttxl.F90 bcttxr.F90 bcttyl.F90 \
    bcttyr.F90 bcttzl.F90 bcttzr.F90 ardump.F90 \
    bcutxl.F90 bcutxr.F90 bcutyl.F90 bcutyr.F90 \
    bcutzl.F90 bcutzr.F90 bcytxl.F90 bcytxr.F90 bcytyl.F90 bcytyr.F90 \
    bcytzr.F90 buftxl.F90 chkarr.F90 chkary.F90 cxhalo.F90 cxhals.F90 \
    cyhalo.F90 cyhals.F90 czhalo.F90 czhals.F90 erfunc.F90 espect.F90 \
    espksq.F90 espovk.F90 finish.F90 bcytzl.F90 dfbydx.F90 d2fdx2.F90 \
    rhscal.F90 lincom.F90 fincom.F90 \
    temper.F90 chrate.F90 chemin.F90 diffin.F90 radcal.F90 radcin.F90 \
    contin.F90 \
    bcinit.F90 bounds.F90 bountt.F90 boundt.F90 \
    turbin.F90 turbft.F90 fftsym.F90 fftf3d.F90 \
    espini.F90 integf.F90 \
    flamin.F90 \
    fftixl.F90

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
	$(FC) $(FFLAGS) -c $<

# Full make of SENGA2
all: senga2

senga2: $(OBJS_F90) $(OBJS) 
	$(FC) $(FFLAGS) $(OBJS) $(OBJS_F90) -o senga2_gnu $(LIBS)

# Add -pg for profiling
#senga2: $(OBJS) 
#	$(FC) $(OBJS) -pg -o senga2 $(LIBS)

# Command to clean up after
clean: 
	rm -f *.o core *.mod

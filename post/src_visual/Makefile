FC = h5pfc
FFLAGS = -O3 -g -Wall -fcheck=all -mcmodel=medium -fbacktrace
SRC = decl_var.f95 chmech.f95 read_cont.f95 read_h5.f95 XDMF_write.f95 HDF5_write.f95 visual.f95
OBJ = ${SRC:.f95=.o}

a.out: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o: %.f95
	$(FC) $(FFLAGS) -o $@ -c $<

clean:
	rm *.mod *.o a.out

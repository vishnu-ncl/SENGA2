#==========================================================================#
#       COMBINE FOR SENGA2 makefile 
#	Created by J LAI
#       j.lai@ncl.ac.uk
#	Last time modified MAY 2015
#==========================================================================#
#   This makefiles creates several programs
#   use command: make ARG=n
#   where n is the program to create
#==========================================================================#
#FC      = /usr/bin/f95
FC      = gfortran

DEPEND	= makedepend
SRCS	= decomp.f\

OBJECTS= decomp.o


LIBES	= -lm

FFLAGS	  = -u -03

FFLAGS  =-mcmodel=large 

accum:	 		$(OBJECTS)
		$(FC) $(FFLAGS) -o decomp $(OBJECTS) $(LIBES)

clean:
		$(RM) $(OBJECTS)

distclean:
		$(RM) $(OBJECTS) core core.* combine

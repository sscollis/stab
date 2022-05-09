#=============================================================================
#  Makefile for stab (IBM XLF on Max OS X)
#
#  Author:  Scott Collis
#
#  Revised: 10-3-96 
#=============================================================================
NPROC    = 1
NAME     = stab 
DEBUG    = -O
FFLAGS   = -qrealsize=8 -qsuppress=cmpmsg -qsuffix=f=f:cpp=f -qfixed=120 \
	   -c $(DEBUG)
F90FLAGS = -qrealsize=8 -qsuppress=cmpmsg -qsuffix=f=f90:cpp=f90 -c $(DEBUG)
OFLAGS   = 
LIB      = -L/Users/sscoll/dist/atlas/lib -llapack -latlas -lg2c \
	   -Wl,-framework,accelerate
F90COMP  = xlf90
FCOMP    = xlf

.SUFFIXES: .f90 

MODS = stuff.o stencils.o

OBJECTS = input.o getmean.o stab.o getmat.o grad.o temporal.o chebyd.o	\
sgengrid.o spatial.o gengrid.o getmean2.o ogengrid.o stokes.o bump.o	\
sgengrid_s.o noreflect.o mtemporal.o mgengrid.o mgetmean.o mspatial.o	\
calch.o mbump.o circh.o

#OBJS1 = fd_temporal.o fd_spatial_m.o
OBJS1 = fd_temporal.o fd_spatial.o

OBJS2 = spline.o piksr2.o rtsafe.o

$(NAME): $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2)
	$(F90COMP) $(OFLAGS) $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2) \
	-o $(NAME) $(LIB)

$(OBJS1): $(MODS)

$(OBJECTS): stuff.o

.f90.o:
	$(F90COMP) $(F90FLAGS) $*.f90

.f.o:
	$(FCOMP) $(FFLAGS) $*.f

clean:
	/bin/rm *.o *.mod

getevec: getevec.o
	$(F90COMP) $(OFLAGS) getevec.o -o getevec

getalpha: getalpha.o
	$(F90COMP) $(OFLAGS) getalpha.o -o getalpha

getab: getab.o
	$(F90COMP) $(OFLAGS) getab.o -o getab

getax: getax.o
	$(F90COMP) $(OFLAGS) getax.o -o getax

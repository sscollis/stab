#=============================================================================
#  Makefile for stab (GCC on Max OS X)
#
#  Author:  Scott Collis
#
#  Revised: 1-17-2020
#=============================================================================
NAME     = stab
DEBUG    =
OPT      = -O2
FFLAGS   = -cpp -ffixed-line-length-120 -fdefault-real-8 -fdefault-double-8 \
           -std=legacy $(DEFINES) $(OPT) $(DEBUG)
F90FLAGS = -cpp -fdefault-real-8 -fdefault-double-8 $(DEFINES) $(OPT) $(DEBUG)
OFLAGS   =
LIB      = -L$(HOME)/local/OpenBLAS/lib -lopenblas
FC       = gfortran
F77      = gfortran

.SUFFIXES: .f90

MODS = stuff.o stencils.o

OBJECTS = input.o getmean.o stab.o getmat.o grad.o temporal.o chebyd.o	\
sgengrid.o spatial.o gengrid.o getmean2.o ogengrid.o stokes.o bump.o	\
sgengrid_s.o noreflect.o mtemporal.o mgengrid.o mgetmean.o mspatial.o	\
calch.o mbump.o circh.o

#OBJS1 = fd_temporal.o fd_spatial_m.o
OBJS1 = fd_temporal.o fd_spatial.o

OBJS2 = spline.o
#
# Optionally use Numerical-Recipes (commercial licensed)
#
ifdef USE_NR
  ifeq ($(LIBNR_DIR),)
    LIBNR_DIR = $(HOME)/git/NR-utilities
  endif
  LIB += -L$(LIBNR_DIR) -lnr
else
  $(info STAB currently requires Numerical-Recipes in FORTRAN routines)
  $(error See README.md for details, build with USE_NR=1)
endif
#
#ifdef USE_NR
#  OBJS2 += nr_rtsafe.o nr_rtflsp.o nr_piksr2.o
#endif
#
ALL:  $(NAME) getevec getalpha getab getax

$(NAME): $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2)
	$(FC) $(OFLAGS) $(LIB) $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2) -o $(NAME)

all: $(NAME) getevec getalpha getab getax

$(OBJS1): $(MODS)

$(OBJECTS): stuff.o

.f90.o:
	$(FC) $(F90FLAGS) -c $*.f90

.f.o:
	$(F77) $(FFLAGS) -c $*.f

clean:
	$(RM) -f *.o *.mod $(NAME) getax getab getevec getalpha

getevec: getevec.o
	$(FC) $(OFLAGS) getevec.o -o getevec

getalpha: getalpha.o
	$(FC) $(OFLAGS) getalpha.o -o getalpha

getab: getab.o
	$(FC) $(OFLAGS) getab.o -o getab

getax: getax.o
	$(FC) $(OFLAGS) getax.o -o getax

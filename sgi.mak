#=============================================================================
#  Makefile for stab on the SGI
#
#  Author:  Scott Collis
#
#  Revised: 10-3-96 
#=============================================================================
NPROC  = 1
NAME   = stab 
DEBUG  = -cpp -OPT:Olimit=0
FFLAGS = -r8 -n32 -O1 -c $(DEBUG)
OFLAGS = -r8 -n32 -O1 $(DEBUG)
LIB    = -lcomplib.sgimath 
COMP   = f90

.SUFFIXES: .f90 

MODS = stuff.o stencils.o

OBJECTS =  input.o getmean.o stab.o getmat.o grad.o \
temporal.o chebyd.o sgengrid.o spatial.o  gengrid.o \
getmean2.o ogengrid.o stokes.o bump.o sgengrid_s.o  \
noreflect.o mtemporal.o mgengrid.o mgetmean.o       \
mspatial.o calch.o mbump.o circh.o exit.o

OBJS1 = fd_temporal.o fd_spatial_m.o

OBJS2 = spline.o piksr2.o rtsafe.o

$(NAME): $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2)
	 $(COMP) $(OFLAGS) $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2) -o $(NAME) $(LIB)

$(MODS):
	$(COMP) $(FFLAGS) $*.f90

$(OBJS1): $(MODS)
	$(COMP) $(FFLAGS) $*.f90

$(OBJECTS): stuff.o
	$(COMP) $(FFLAGS) $*.f90

.f90.o:
	$(COMP) $(FFLAGS) $*.f90 

.f.o:
	$(COMP) $(FFLAGS) $*.f

clean:
	/bin/rm *.o

getevec: getevec.o
	f90 $(OFLAGS) getevec.o -o getevec

getalpha: getalpha.o
	f90 $(OFLAGS) getalpha.o -o getalpha

getab: getab.o
	f90 $(OFLAGS) getab.o -o getab

getax: getax.o
	f90 $(OFLAGS) getax.o -o getax

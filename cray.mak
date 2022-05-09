#==============================================================================
#  Makefile for stab on the Cray C90
#
#  Author:  Scott Collis
#
#  Revised: 10-3-96 
#==============================================================================
NPROC  = 1
NAME   = stab 
DEBUG  =
FFLAGS = -eZ -DCRAY -c $(DEBUG)
OFLAGS = $(DEBUG) -o $(NAME)
LIB    = -L/usr/local/lib -limsl -llapack
COMP   = f90

.SUFFIXES: .f90 

MODS = stuff.o stencils.o

OBJECTS =  input.o getmean.o stab.o getmat.o grad.o \
temporal.o chebyd.o sgengrid.o spatial.o  gengrid.o \
getmean2.o ogengrid.o stokes.o bump.o sgengrid_s.o  \
noreflect.o mtemporal.o mgengrid.o mgetmean.o       \
mspatial.o calch.o mbump.o circh.o

OBJS1 = fd_temporal.o fd_spatial.o

OBJS2 = spline.o piksr2.o rtsafe.o

$(NAME): $(MODS) $(OBJECTS) $(OBJS1) $(OBJS2)
	 $(COMP) $(OFLAGS) $(OBJECTS) $(OBJS1) $(OBJS2) $(LIB)

($MODS):
	$(COMP) $(FFLAGS) $*.f90

$(OBJS1): $(MODS)
	$(COMP) $(FFLAGS) -p stuff.o -p stencils.o $*.f90

$(OBJECTS): stuff.o
	$(COMP) $(FFLAGS) -p stuff.o $*.f90

.f90.o:
	 $(COMP) $(FFLAGS) $*.f90

.f.o:
	 cf77 -c $(DEBUG) $*.f

clean:
	/bin/rm *.o

getevec: getevec.o
	f90 getevec.o -o getevec
	
getalpha: getalpha.o
	f90 getalpha.o -o getalpha

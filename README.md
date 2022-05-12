# STAB:  Compressible stability solver

Updated to build with gfortran on Darwin

Currently uses a couple of Numerical-Recepies in FORTRAN routines that are
commercially licensed and not distributed here.

To build use:

    ln -s gcc.mak Makefile
    make USE_NR=1

It would be easy to remove the dependency on NR and I encourage someone to do
so.

S. Scott Collis

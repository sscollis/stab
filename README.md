# STAB:  Compressible stability solver

This uses a direct approach using either Chebyshev collocation or high-order
finite difference.   This is complemented by the `shoot` solver which uses
and ODE solver (shooting).

Updated to build with the GCC gfortran compiler on MacOS and Linux.

Currently uses a couple of Numerical-Recipes in FORTRAN routines that are
commercially licensed and not distributed here.

To build use:
```bash
    ln -s gcc.mak Makefile
    make USE_NR=1
```

Notes:
1.  It would be easy to remove the dependency on NR and I encourage 
    someone to do so.
2.  You may **not** add the NR software to this repository and you may only use
    it if you have a valid license to do so

S. Scott Collis

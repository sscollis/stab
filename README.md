# STAB:  Compressible stability solver

This uses a direct approach using either Chebyshev collocation or high-order
finite difference.   This is complemented by the `shoot` solver which uses
an ODE solver with Conte orthogonalization and shooting to satisfy the boundary
conditions.  I suggest using `shoot` to polish the eigensolutions identified
by `stab`.

The `stab` code has been updated to build with the GCC gfortran compiler on
both MacOS and Linux and older builds have not been recently used.

## Building

To build use:
```bash
    ln -s gcc.mak Makefile
    make USE_NR=1
```

## Notes
`stab` currently uses a couple of Numerical-Recipes in FORTRAN
routines that are commercially licensed and therefore not distributed here.
You can easily replace these with public versions and are encouraged to do so.

1.  It would be easy to remove the dependency on NR and I encourage
    someone to do so.
2.  You may **not** add the NR software to this repository and you may only use
    it if you have a valid license to do so.

## Test Cases

Case       |  Description
-----------|---------------------------------------------------------------------------------------
`thesis`   |  Test case from Collis Ph.D. thesis, Chapter 4
`compbl`   |  Similar to `thesis` but using `compbl` instead of `fsc` for mean, slight differences
`TStest`   |  Tollmein-Schlichting test case
`CFtest`   |  Cross-flow vortex test case
`mixl`     |  Mixing layer test case
`test`     |  Test case for experimentation
`vk`       |  Case for Vinod Kumar

## Contact

S. Scott Collis
scollis@gmail.com

# Crossflow test case

Note `evec.ref` is in `big_endian` (written on SGI) so on an Intel machine you need to i
read it using

```bash
GFORTRAN_CONVERT_UNIT='swap' ../getevec
```

The eigenvalue from `sgi` uses $N_y = 64$ reference is:

$$\alpha = -3.7392297537875E+001  -2.9982641664422E-001 $$

whereas from Intel on MacOS is using $N_y=96$ is:

$$\alpha = -3.7392546274473E+001  -2.9977113927779E-001 $$

S. Scott Collis\
sscollis@gmail.com

# Tollmein-Schilichting Test case

This is the TS test case from Collis Ph.D. Thesis, Chapter 4.
It is computing using the `fsc` compressible Falkner-Skan solver 
for the mean flow which matches the thesis results exactly.

## Spatial results

Note that these results were likely computed in quad precision on
a Cray (Cray's rocked).

```bash
FSC: M=0.3, lambda=0, beta_h=0, Tw/T0=1, Re_\delta_1 = 1000, Pr=1

k_x = 2.2804739410500E-001  -6.5163146761218E-003  (yi=1.0, 64)
k_x = 2.2804739411367E-001  -6.5163146912049E-003  (yi=1.0, 96)

\lambda_{TS} = 27.552103
```

## Temporal results

```bash
FSC: M=0.3, lambda=0, beta_h=0, Tw/T0=1, Re_\delta_1 = 1000, Pr=1

omega = 1.1449048146419E-001   2.1951514140605E-003  (yi=1.0, 64)
omega = 1.1449048146373E-001   2.1951514138197E-003  (yi=1.0, 96)

k_x = 0.308620690
```

## Running the test

The script `run.sh` computes the mean flow using `fsc`, runs both
the temporal and spatial programs using `stab` and then polishes
the spatial results using `shoot`.  Note that it assumes locations
for each of these solves to you may need to adjust in practise.

To plot the streamwise velocity for the spatial eigenfunction using `gnuplot`

```bash
gnuplot
set xrange [0:15]
set xlabel "y"
set ylabel "u_r, u_i"
plot "space.1" u 1:4 w l title "u_r", "space.1" u 1:5 w l title "u_i"
```

or use the command file

```bash
gnuplot
load "plot.com"
```

S. Scott Collis

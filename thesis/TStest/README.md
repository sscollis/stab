# Tollmein-Schilichting Test case

This is the TS test case from Collis Ph.D. Thesis, Chapter 4.
It is computing using the `fsc` compressible Falkner-Skan solver 
for the mean flow which matches the thesis results exactly.

## Spatial results

```bash
FSC: M=0.3, lambda=0, beta_h=0, Tw/T0=1, Re_\delta_1 = 1000, Pr=1

k_x = 2.2804739411412E-001  -6.5163146827854E-003  (yi=1.0, 64)
k_x = 2.2804739411185E-001  -6.5163146952733E-003  (yi=1.0, 96)

\omega = 0.08
\lambda_{TS} = 2\pi/k_x = 27.552103068969444
```

## Temporal results

```bash
FSC: M=0.3, lambda=0, beta_h=0, Tw/T0=1, Re_\delta_1 = 1000, Pr=1

omega = 1.1467880189148E-001  2.3844535289045E-003 (yi=1.0, 64)
omega = 1.1467880189410E-001  2.3844535276592E-003 (yi=1.0, 96)

k_x = 0.308620690
\lambda_{TS} = 2\pi/k_x = 20.358924436270254
```

## Running the test

The script `run.sh` computes the mean flow using `fsc`, runs both
the temporal and spatial programs using `stab` and then polishes
the spatial results using `shoot`.  Note that it assumes locations
for each of these solves to you may need to adjust in practise.

```bash
./run.sh
```
## Plotting the results

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
## Outputs

Additional output files include:

File        |    Description
------------|--------------------------------------------
time.1      |    Temporal eigenfunction from `stab`
space.1     |    Spatial eigenfunction from `stab`
efun.out    |    Polished spatial eigenfunction from `shoot`
adj.out     |    Adjoint spatial eigenfunction from `shoot`
output.dat  |    Full `shoot` output with nonparallel corrections
evec.dat    |    The full eigensolution output (from last run of `stab`)
rho.out     |    Densitymean profile and derivatives
u.out       |    Streamwise mean profile and derivatives
w.out       |    Crossflow mean profile and derivatives
t.out       |    Temperature mean profile and derivatives

S. Scott Collis
sscollis@gmail.com

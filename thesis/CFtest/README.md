# Crossflow Test case

This is the CF test case from Collis Ph.D. Thesis, Chapter 4.
It is computing using the `fsc` compressible Falkner-Skan solver 
for the mean flow which matches the thesis results exactly.

## Spatial results

## Temporal results

```bash
M = 0.3, Re = 400, Pr = 1.0

theta = 45
beta_h = 1

T_w = T_0 

k_x = -0.287436451
k_z = 0.35
\omega_{lst} = 0.006533585 i

\Delta t = 0.047643

N_x = 20 N_y = 127
N_x = 80 N_y = 127
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

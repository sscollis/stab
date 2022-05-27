# Tollmein-Schilichting Test case

This is the TS test case *similar* to Collis Ph.D. Thesis, Chapter 4
but instead of using FSC it uses the COMPBL solver.  As such, this gives 
slightly different results.

## Spatial results

```bash
FSC: M=0.3, lambda=0, beta_h=0, Tw/T0=1, Re_\delta_1 = 1000, Pr=1

k_x = 2.2822706636765E-001  -6.1126970532663E-003  (yi=1.0, 64)
k_x = 2.2822706636782E-001  -6.1126970666343E-003  (yi=1.0, 96)

\lambda_{TS} = 27.530412615712063 
```

## Temporal results

```bash
FSC: M=0.3, lambda=0, beta_h=0, Tw/T0=1, Re_\delta_1 = 1000, Pr=1

omega = 1.1449048146417E-001   2.1951514144492E-003  (yi=1.0, 64)
omega = 1.1449048146419E-001   2.1951514140605E-003  (yi=1.0, 96)

k_x = 0.308620690
\lambda_{TS} = 20.358924436270254
```

## Running

To run these simply using the script `run.sh` but note that paths
to `compbl` and `shoot` are hardwired and may need to be altered for
your configuration.

```bash
./run.sh
```

## Outputs

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

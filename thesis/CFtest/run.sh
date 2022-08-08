#!/bin/bash
echo $FSC_DIR
FSC_DIR=${FSC_DIR:-../../../fsc}
$FSC_DIR/fsc < fsc.inp 
#../../../fsc/fsc < fsc.inp
tail -n +2 cprofile.dat > profile.0
tail -n +2 cfirst.dat > first.0
tail -n +2 csecond.dat > second.0
#
# Temporal problem
#
../../stab < temporal.inp
../../getevec -v << EOF
evec.dat
6.3418480187508E-007   6.5335847258858E-003
0 0
EOF
#
# Spatial problem
#
../../stab < spatial.inp
../../getevec -v << EOF
evec.dat
-2.8831962907615E-001  -1.3854663677328E-002
0 0
EOF
#
# Shooting to cleanup spatial 
#
#../../../shoot/shoot < shoot.inp
#
# check quality of solution 
#
ndiff -abserr 1e-8 -relerr 1e-8 space.1 space.ref &&
ndiff -abserr 1e-8 -relerr 1e-8 time.1 time.ref
status=$?
echo
echo 'ndiff completed with status:' $status
exit $status 

#!/bin/bash
../../../fsc/fsc < blasius.inp
tail -n +2 cprofile.dat > profile.0
tail -n +2 cfirst.dat > first.0
tail -n +2 csecond.dat > second.0
#
# Temporal problem
#
../../stab < temporal.inp
../../getevec -v << EOF
evec.dat
1.1467880189410E-001   2.3844535276599E-003
0 0
EOF
#
# Spatial problem
#
../../stab < spatial.inp
../../getevec -v << EOF
evec.dat
2.2804739411180E-001  -6.5163146952626E-003
0 0
EOF
#
# Shooting to cleanup spatial 
#
../../../shoot/shoot < shoot.inp
#
# check quality of solution 
#
echo "============================================="
echo "ndiff of space" && \
ndiff -abserr 1e-8 space.1 space.ref | tee space.log && \
echo "=============================================" && \
echo "ndiff of time" && \
ndiff -abserr 1e-8 time.1 time.ref | tee time.log 
status=$?
echo "============================================="
echo 'ndiff completed with status:' $status
exit $status 

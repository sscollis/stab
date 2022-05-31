#!/bin/bash
../../fsc/fsc < blasius.inp
tail -n +2 cprofile.dat > profile.0
tail -n +2 cfirst.dat > first.0
tail -n +2 csecond.dat > second.0
#
# Temporal problem
#
../stab < temporal.inp
../getevec << EOF
evec.dat
480
0
EOF
#
# Spatial problem
#
../stab < spatial.inp
../getevec << EOF
evec.dat
371
0
EOF
#
# Shooting to cleanup spatial 
#
../../shoot/shoot < shoot.inp
#
exit $?

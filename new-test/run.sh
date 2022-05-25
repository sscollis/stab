#!/bin/bash
../../compbl/compbl < compbl.inp
\mv header.dat header.0
\mv profile.dat profile.0
\mv first.dat first.0
\mv second.dat second.0
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
exit $?

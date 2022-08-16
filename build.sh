#!/bin/bash
#
# Build a fresh version
#
# Revised:  8/15/22
# Author:   S.Scott Collis 
#
set -e
make clean && make $@
exit $? 

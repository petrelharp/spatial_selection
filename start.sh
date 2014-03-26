#!/bin/bash

# USAGE:
#   start.sh 'nsteps=100' 'stepsize=1' 'range=c(201,201)'

ARGS=""
while (( $# > 0 ))
do
    ARGS="$ARGS '"$1"'"
    shift
done

qsub -vARGS="$ARGS" run-sim.pbs 


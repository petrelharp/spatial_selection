#!/bin/bash

# USAGE:
#   restart.sh 7702054-r1-101-sb0.05-sm-0.01-pophistory-run.Rdata 'nsteps=100' 'stepsize=1'

ARGS="'restart=\"$1\"'"
shift
while (( $# > 0 ))
do
    ARGS="$ARGS '"$1"'"
    shift
done
if [ $(echo $ARGS | grep -v burnin) ]
then
    ARGS="$ARGS 'burnin=0'"
fi

qsub -vARGS="$ARGS" run-sim.pbs 

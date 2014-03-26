#!/bin/bash

# USAGE:
#   restart.sh 7702054-r1-101-sb0.05-sm-0.01-pophistory-run.Rdata 'nsteps=100' 'stepsize=1'

qsub -vARGS="'restart=\"$1\"' $*" run-sim.pbs 

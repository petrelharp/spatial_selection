#!/bin/bash

# Simulate time to migration in a new patch in 2D
SB=".01"
N=$1  # do 100 1000

for SM in "-.0001" "-.001" "-.01" "-.1"
do
    for R in 10 20 30 # distance between patches
    do
        THISDIR="migration/sm${SM}/R-${R}/N-${N}/2D"
        mkdir -p $THISDIR
        Rscript generate-patchy-run.R \
            "outdir=\"${THISDIR}\"" \
            "N=${N}" "range=c(101,101)" "sb=${SB}" "sm=${SM}" "nsteps=1000" "stepsize=25" "ntypes=2" \
            "s=matrix(sm,nrow=range[1],ncol=range[2]);s[(abs(abs(row(s)-range[2]/2)-(${R}/2+5))<5)&(col(s)>44)&(col(s)<56)]=sb" \
            'initdist=array(c(rep(N,prod(range)),rep(0,(ntypes-1)*prod(range))),c(range,ntypes));initdist[(row(s)<range[2]/2)&s>0]=floor(.8*N)' &
    done
done

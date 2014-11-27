#!/bin/bash

# Simulate time to migration in a new patch
SB=".01"

for SM in "-.003" "-.03"
do
    for N in 100 1000
    do
        for R in 60 100 120 # distance between patches
        do
            THISDIR="migration/sm${SM}/R-${R}/N-${N}/1D"
            mkdir -p $THISDIR
            Rscript generate-patchy-run.R \
                "outdir=\"${THISDIR}\"" \
                "N=${N}" "range=c(1,501)" "sb=${SB}" "sm=${SM}" "nsteps=1000" "stepsize=25" "ntypes=2" \
                "s=matrix(ifelse(abs(abs(seq.int(1,range[2])-range[2]/2)-(${R}/2+20))<20,sb,sm),nrow=1,ncol=range[2])" \
                'initdist=array(c(rep(N,prod(range)),rep(0,(ntypes-1)*prod(range))),c(range,ntypes));initdist[(seq.int(1,range[2])<range[2]/2)&s>0]=floor(.8*N)' &
        done
    done
done


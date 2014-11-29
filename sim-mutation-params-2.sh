#!/bin/bash

# Simulate time to mutation in a new patch
SB=".01"

for MU in "1e-5" "1e-6"
do
    for SM in "-.001" "-.01" "-.1"
    do
        for N in 50 200 600 800 1000 1200
        do
            THISDIR="mutation/mu${MU}/sm${SM}/N-${N}/1D"
            mkdir -p $THISDIR
            Rscript generate-patchy-run.R \
                "outdir=\"${THISDIR}\"" \
                "N=${N}" "range=c(1,501)" "sb=${SB}" "sm=${SM}" "nsteps=1000" "stepsize=25" "ntypes=2" "mu=1e-5"\
                's=matrix(ifelse((seq.int(1,range[2])>200)&(seq.int(1,range[2])<300),sb,sm),nrow=1,ncol=range[2])' \
                'initdist=array(c(rep(N,prod(range)),rep(0,(ntypes-1)*prod(range))),c(range,ntypes))' &
        done
    done
done

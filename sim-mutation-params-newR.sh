#!/bin/bash

# Simulate time to mutation in a new patch
SB=".01"
SM="-.0001"

for SM in "-.0001" "-.001" "-.01" "-.1"
do
    for N in 100 1000
    do
        for MU in "1e-5" "1e-6"
        do
            THISDIR="mutation/mu${MU}/sm${SM}/N-${N}/1D"
            mkdir -p $THISDIR
            Rscript generate-patchy-run.R \
                "outdir=\"${THISDIR}\"" \
                "N=${N}" "range=c(1,501)" "sb=${SB}" "sm=${SM}" "nsteps=1000" "stepsize=25" "ntypes=2" "mu=${MU}"\
                's=matrix(ifelse((seq.int(1,range[2])>200)&(seq.int(1,range[2])<300),sb,sm),nrow=1,ncol=range[2])' \
                'initdist=array(c(rep(N,prod(range)),rep(0,(ntypes-1)*prod(range))),c(range,ntypes))' &
        done
    done
done


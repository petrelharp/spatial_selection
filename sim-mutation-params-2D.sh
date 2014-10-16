#!/bin/bash

# Simulate time to mutation in a new patch in 2D
SB=".01"

for SM in "-.0001" "-.001" "-.01" "-.1"
do
    for N in 25 100 400 1600
    do
        THISDIR="mutation/sm${SM}/N-${N}/2D"
        mkdir -p $THISDIR
        Rscript generate-patchy-run.R \
            "outdir=\"${THISDIR}\"" \
            "N=${N}" "range=c(101,101)" "sb=${SB}" "sm=${SM}" "nsteps=1000" "stepsize=25" "ntypes=2" "mu=1e-5"\
            's=matrix(sm,nrow=range[1],ncol=range[2]);s[(row(s)>44)&(row(s)<56)&(col(s)>44)&(col(s)<56)]=sb' \
            'initdist=array(c(rep(N,prod(range)),rep(0,(ntypes-1)*prod(range))),c(range,ntypes))' &
    done
done

#!/bin/bash

# Simulate both
SB=".01"
N=$1  # do 100 1000

for SM in "-.0001" "-.001" "-.01" "-.1"
do
    for R in 10 20 40 80 160  # distance between patches
    do
        THISDIR="both/sm${SM}/R-${R}/N-${N}/1D"
        mkdir -p $THISDIR
        Rscript generate-patchy-run.R \
            "outdir=\"${THISDIR}\"" \
            "N=${N}" "range=c(1,501)" "sb=${SB}" "sm=${SM}" "nsteps=1000" "stepsize=25" "ntypes=2" mu="1e-5"\
            "s=matrix(ifelse(abs(abs(seq.int(1,range[2])-range[2]/2)-(${R}/2+20))<20,sb,sm),nrow=1,ncol=range[2])" \
            'initdist=array(c(rep(N,prod(range)),rep(0,(ntypes-1)*prod(range))),c(range,ntypes))' &
    done
done

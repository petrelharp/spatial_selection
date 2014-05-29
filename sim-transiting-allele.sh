#!/bin/bash

# Get a simulation of transits to a new patch
Rscript generate-patchy-run.R 'N=1000' 'range=c(1,301)' 'sb=0.05' 'sm=-0.006' 'nsteps=20' 'stepsize=100' 's=matrix(ifelse((seq.int(range[1],range[2])>185)&(seq.int(range[1],range[2])<205),sb,sm),nrow=1,ncol=range[2])' 'run.id=12345'
Rscript generate-patchy-run.R 'restart="12345-r1-301-sb0.05-sm-0.006-pophistory-run.Rdata"' 'sb=0.05' 'sm=-0.006' 'nsteps=100000' 'stepsize=1' 'range=c(1,301)' 's=matrix(ifelse(((seq.int(range[1],range[2])>185)&(seq.int(range[1],range[2])<205))|seq.int(range[1],range[2])<20,sb,sm),nrow=1,ncol=range[2])' 'run.id=67890'

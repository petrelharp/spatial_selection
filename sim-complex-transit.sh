#!/bin/bash

# ntypes=3
# s=matrix(ifelse( seq.int(1,range[2])<50 | ((seq.int(1,range[2])>160)&(seq.int(1,range[2])<190)) | ((seq.int(1,range[2])>310)&(seq.int(1,range[2])<340)) | seq.int(1,range[2])>450 ,sb,sm) ,nrow=1,ncol=range[2])
# initdist=array(0,dim=c(range,3));
# initdist[1,1,2]=1;
# initdist[1,range[2],3]=1;
# initdist[1,1,1]=initdist[1,range[2],1]=N-1

# Get a simulation of transits to a new patch
Rscript generate-patchy-run.R \
    'N=200' 'range=c(1,501)' 'sb=0.05' 'sm=-0.018' 'nsteps=1000' 'stepsize=100' 'ntypes=3' \
    's=matrix(ifelse(seq.int(1,range[2])<70|((seq.int(1,range[2])>130)&(seq.int(1,range[2])<220))|((seq.int(1,range[2])>280)&(seq.int(1,range[2])<370))|seq.int(1,range[2])>430,sb,sm),nrow=1,ncol=range[2])' \
    'initdist=array(0,c(range,3));initdist[1,1:30,2]=40;initdist[1,range[2]-(0:29),3]=40;initdist[,,1]=N;initdist[1,1:30,1]=initdist[1,range[2]-(0:29),1]=N-40'

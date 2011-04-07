
# makefile for making plots

proportion-by-sigma.eps: standing-variation-plots.R
    nice -19 R --vanilla < $< 



proportion-by-sigma.eps proportion-by-sd.eps charlen-by-sd-limit.eps charlen-by-sd-sigma.eps: Spatial_adaptation/standing-variation-plots.R Spatial_adaptation/standing-variation-fns.R
	nice -19 R --vanilla < $< 


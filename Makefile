
proportion-by-sigma.eps proportion-by-sd.eps charlen-by-sd-limit.eps charlen-by-sd-sigma.eps: Spatial_adaptation/standing-variation-plots.R
	nice -19 R --vanilla < $< 


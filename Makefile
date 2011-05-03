
proportion-by-sigma.eps charlen-by-sd.eps: Spatial_adaptation/standing-variation-plots.R
	nice -19 R --vanilla < $< 


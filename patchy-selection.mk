
plots = 

$(plots): Spatial_adaptation/patchy-selection-plots.R Spatial_adaptation/patchy-selection-fns.R
	nice -19 R --vanilla < $< 


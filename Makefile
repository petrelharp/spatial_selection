
plots = proportion-by-sigma.eps proportion-by-sd.eps charlen-by-sd-limit.eps charlen-by-sd-sigma.eps

patchyplots = phase-diagram-log.pdf phase-diagram.pdf prob-establishment.pdf

$(plots): Spatial_adaptation/standing-variation-plots.R Spatial_adaptation/standing-variation-fns.R
	nice -19 R --vanilla < $< 

$(patchyplots): patchy-selection-plots.R patchy-selection-fns.R plotting-fns.R
	nice -19 R --vanilla < patchy-selection-plots.R

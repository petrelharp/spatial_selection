
panmixia = min-s-versus-r.pdf f-versus-gamma.pdf gamma-contour-1e-8.pdf

examples.tex: Spatial_adaptation/compute-things.R
	nice -19 R --vanilla < $< 

$(plots): Spatial_adaptation/panmixia.R
	nice -19 R --vanilla < $<


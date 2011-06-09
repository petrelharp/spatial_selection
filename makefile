
panmixia = min-s-versus-r.png f-versus-gamma.png gamma-contour-1e-8.png

examples.tex: Spatial_adaptation/compute-things.R
	nice -19 R --vanilla < $< 

$(panmixia): Spatial_adaptation/panmixia.R
	nice -19 R --vanilla < $<


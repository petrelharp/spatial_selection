
panmixia = Spatial_adaptation/min-s-versus-r.pdf Spatial_adaptation/f-versus-gamma.pdf Spatial_adaptation/gamma-contour-1e-8.pdf

examples.tex: Spatial_adaptation/compute-things.R
	nice -19 R --vanilla < $< 

$(panmixia): Spatial_adaptation/panmixia.R
	nice -19 R --vanilla < $<


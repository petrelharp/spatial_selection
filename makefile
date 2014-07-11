

PATCHY_FIGS = $(shell grep includegr patchy-selection-paper.tex | sed -e 's/.*{\(.*\)}.*/\1/' )
PATCHY_EPS = $(patsubst %,%.eps,$(PATCHY_FIGS))

.PHONY : figs

figs : $(PATCHY_EPS)

panmixia = min-s-versus-r.pdf f-versus-gamma.pdf gamma-contour-1e-8.pdf

examples.tex: Spatial_adaptation/compute-things.R
	nice -19 R --vanilla < $< 

$(panmixia): Spatial_adaptation/panmixia.R
	nice -19 R --vanilla < $<

%.eps : %.pdf
	inkscape --without-gui --export-eps=$@ $<


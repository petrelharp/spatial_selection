
plots = helianthus-ex-table.tex helianthus-standing-ex-table.tex

$(plots): patchy-selection-plots.R patchy-selection-fns.R
	nice -19 R --vanilla < $< 


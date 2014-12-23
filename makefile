PATCHY_TARGETS = patchy-selection-paper-with-figs.pdf patchy-selection-review-responses.pdf patchy-selection-paper-no-figs.pdf patchy-selection-paper-arxiv.pdf patchy-tab-S1.pdf patchy-tab-S2.pdf patchy-supp-info.pdf

sfigs = patchy-fig-S3.pdf patchy-fig-S4.pdf patchy-fig-S5.pdf patchy-fig-S6.pdf patchy-fig-S7.pdf patchy-fig-S8.pdf patchy-fig-S9.pdf patchy-fig-S10.pdf

PATCHY_FIGS = $(shell grep "^ *[^%].*includegr" patchy-selection-paper.tex | sed -e 's/.*{\([^}]*\).*/\1/' )
PATCHY_EPS = $(patsubst %,%.eps,$(PATCHY_FIGS))
SFIG_EPS = $(patsubst %.pdf,%.eps,$(sfigs))

.PHONY : figs patchy clean

all : patchy figs
	# tar -cvzhf patchy-suppmat.tar.gz patchy-supp-info.pdf S*_Table.pdf S*_Figure.pdf

patchy : $(PATCHY_TARGETS)

clean :
	-rm -f *.aux
	-rm -f *.log
	-rm -f *.bbl
	-rm -f *.blg
	-rm -f *.fff
	-rm -f *.lof
	-rm -f *.lot
	-rm -f *.out
	-rm -f *.ttt

superclean : clean
	rm -f $(PATCHY_TARGETS)

figs : $(PATCHY_EPS) $(SFIG_EPS)
	-ln -f -s sim-occupation-freqs.eps Fig1.eps
	-ln -f -s sim-snapshots.eps Fig2.eps
	-ln -f -s branching-concept.eps  Fig3.eps
	-ln -f -s phase-diagram-log.eps Fig4.eps
	-ln -f -s sim-transit.eps Fig5.eps
	-ln -s -f prob-mutation-compared.eps Fig6.eps
	-ln -f -s Lava_flow_mice_prob_parallel.eps Fig7.eps
	-ln -f -s prob-establishment.eps Fig8.eps
	-ln -f -s mutation-times-predicted.pdf S1_Figure.pdf
	-ln -f -s migration-time-predicted.pdf S2_Figure.pdf
	-ln -f -s patchy-fig-S3.pdf S3_Figure.pdf
	-ln -f -s patchy-fig-S4.pdf S4_Figure.pdf
	-ln -f -s patchy-fig-S5.pdf S5_Figure.pdf
	-ln -f -s patchy-fig-S6.pdf S6_Figure.pdf
	-ln -f -s patchy-fig-S7.pdf S7_Figure.pdf
	-ln -f -s patchy-fig-S8.pdf S8_Figure.pdf
	-ln -f -s patchy-fig-S9.pdf S9_Figure.pdf
	-ln -f -s patchy-fig-S10.pdf S10_Figure.pdf
	-ln -f -s patchy-tab-S1.pdf S1_Table.pdf
	-ln -f -s patchy-tab-S2.pdf S2_Table.pdf

patchy-selection-paper-submission-no-figs.pdf : patchy-selection-paper.tex standing_patches_refs.bib patchy-review-responses.tex
	-rm -f *.{aux,bbl}
	pdflatex -jobname patchy-selection-paper-submission-no-figs '\def\submissionnofigs{a} \input{patchy-selection-paper}'
	bibtex patchy-selection-paper-submission-no-figs
	pdflatex -jobname patchy-selection-paper-submission-no-figs '\def\submissionnofigs{b} \input{patchy-selection-paper}'
	pdflatex -jobname patchy-selection-paper-submission-no-figs '\def\submissionnofigs{c} \input{patchy-selection-paper}'
	pdflatex -jobname patchy-selection-paper-submission-no-figs '\def\submissionnofigs{d} \input{patchy-selection-paper}'

patchy-selection-paper-submission-figs.pdf : patchy-selection-paper.tex standing_patches_refs.bib patchy-review-responses.tex
	-rm -f *.{aux,bbl}
	pdflatex -jobname patchy-selection-paper-submission-figs '\def\submissionfigs{a} \input{patchy-selection-paper}'
	bibtex patchy-selection-paper-submission-figs
	pdflatex -jobname patchy-selection-paper-submission-figs '\def\submissionfigs{a} \input{patchy-selection-paper}'
	pdflatex -jobname patchy-selection-paper-submission-figs '\def\submissionfigs{a} \input{patchy-selection-paper}'
	pdflatex -jobname patchy-selection-paper-submission-figs '\def\submissionfigs{a} \input{patchy-selection-paper}'

patchy-selection-paper-arxiv.pdf : patchy-selection-paper.tex standing_patches_refs.bib patchy-review-responses.tex
	-rm -f *.{aux,bbl}
	pdflatex -jobname patchy-selection-paper-arxiv '\def\arxiv{a} \input{patchy-selection-paper}'
	bibtex patchy-selection-paper-arxiv
	pdflatex -jobname patchy-selection-paper-arxiv '\def\arxiv{a} \input{patchy-selection-paper}'
	pdflatex -jobname patchy-selection-paper-arxiv '\def\arxiv{a} \input{patchy-selection-paper}'
	pdflatex -jobname patchy-selection-paper-arxiv '\def\arxiv{a} \input{patchy-selection-paper}'

# Submit: paper with figures in reasonable order, and page and line numbers referred to in review-reponses below.
patchy-selection-paper-with-figs.pdf : patchy-selection-paper-submission-figs.pdf
	pdfjam --outfile $@ $< 1-44

# Submit: responses to reviews
patchy-selection-review-responses.pdf : patchy-selection-paper-submission-figs.pdf
	pdfjam --outfile $@ $< 45-

# Submit: the paper without figures but lists of figure legends
patchy-selection-paper-no-figs.pdf : patchy-selection-paper-submission-no-figs.pdf
	pdfjam --outfile $@ $< 1-28

# Submit: the "appendix" text
patchy-supp-info.pdf : patchy-selection-paper-submission-no-figs.pdf
	pdfjam --outfile $@ $< 29-30

# Submit: table S1
patchy-tab-S1.pdf : patchy-selection-paper-submission-no-figs.pdf
	pdfjam --outfile $@ $< 39

# Submit: table S2
patchy-tab-S2.pdf : patchy-selection-paper-submission-no-figs.pdf
	pdfjam --outfile $@ $< 40

patchy-fig-S3.pdf : example-mutation-sims/18885-r1-501-sb0_01-sm-0_1-N1200-pophistory-run.pdf example-mutation-sims/56325-r1-501-sb0_01-sm-0_1-N600-pophistory-run.pdf example-mutation-sims/28432-r1-501-sb0_01-sm-0_1-N50-pophistory-run.pdf
	pdfjam --outfile $@ $? --noautoscale false --nup 1x3

patchy-fig-S4.pdf : example-mutation-sims/43099-r1-501-sb0_01-sm-0_01-N1200-pophistory-run.pdf example-mutation-sims/69787-r1-501-sb0_01-sm-0_01-N600-pophistory-run.pdf example-mutation-sims/11821-r1-501-sb0_01-sm-0_01-N50-pophistory-run.pdf
	pdfjam --outfile $@ $? --noautoscale false --nup 1x3

patchy-fig-S5.pdf : example-mutation-sims/59611-r1-501-sb0_01-sm-0_001-N1200-pophistory-run.pdf example-mutation-sims/93713-r1-501-sb0_01-sm-0_001-N600-pophistory-run.pdf example-mutation-sims/29850-r1-501-sb0_01-sm-0_001-N50-pophistory-run.pdf
	pdfjam --outfile $@ $? --noautoscale false --nup 1x3

patchy-fig-S6.pdf : example-mutation-sims/15449-r1-501-sb0_01-sm-1e-04-N1600-pophistory-run.pdf example-mutation-sims/5582-r1-501-sb0_01-sm-1e-04-N400-pophistory-run.pdf example-mutation-sims/24639-r1-501-sb0_01-sm-1e-04-N25-pophistory-run.pdf
	pdfjam --outfile $@ $? --noautoscale false --nup 1x3

patchy-fig-S7.pdf : example-migration-sims/89826-r1-501-sb0_01-sm-0_1-N4000-pophistory-run.pdf example-migration-sims/76178-r1-501-sb0_01-sm-0_1-N1000-pophistory-run.pdf example-migration-sims/97545-r1-501-sb0_01-sm-0_01-N100-pophistory-run.pdf
	pdfjam --outfile $@ $? --noautoscale false --nup 1x3

patchy-fig-S8.pdf : example-migration-sims/42080-r1-501-sb0_01-sm-0_01-N4000-pophistory-run.pdf example-migration-sims/82698-r1-501-sb0_01-sm-0_01-N1000-pophistory-run.pdf example-migration-sims/97545-r1-501-sb0_01-sm-0_01-N100-pophistory-run.pdf
	pdfjam --outfile $@ $? --noautoscale false --nup 1x3

patchy-fig-S9.pdf : example-migration-sims/88079-r1-501-sb0_01-sm-0_001-N4000-pophistory-run.pdf example-migration-sims/3464-r1-501-sb0_01-sm-0_001-N1000-pophistory-run.pdf example-migration-sims/37774-r1-501-sb0_01-sm-0_001-N100-pophistory-run.pdf
	pdfjam --outfile $@ $? --noautoscale false --nup 1x3

patchy-fig-S10.pdf : example-migration-sims/26044-r1-501-sb0_01-sm-1e-04-N4000-pophistory-run.pdf example-migration-sims/85750-r1-501-sb0_01-sm-1e-04-N1000-pophistory-run.pdf example-migration-sims/4111-r1-501-sb0_01-sm-1e-04-N100-pophistory-run.pdf
	pdfjam --outfile $@ $? --noautoscale false --nup 1x3

panmixia = min-s-versus-r.pdf f-versus-gamma.pdf gamma-contour-1e-8.pdf

examples.tex: Spatial_adaptation/compute-things.R
	nice -19 R --vanilla < $< 

$(panmixia): Spatial_adaptation/panmixia.R
	nice -19 R --vanilla < $<

%.eps : %.pdf
	inkscape --without-gui --export-eps=$@ $<


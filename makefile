

.PHONY : figs patchy clean

all : patchy figs

patchy : patchy-selection-paper-with-figs.pdf patchy-selection-review-responses.pdf patchy-selection-paper-no-figs.pdf patchy-selection-paper-arxiv.pdf

sfigs = patchy-fig-S3.pdf patchy-fig-S4.pdf patchy-fig-S5.pdf patchy-fig-S6.pdf patchy-fig-S7.pdf patchy-fig-S8.pdf patchy-fig-S9.pdf patchy-fig-S10.pdf

PATCHY_FIGS = $(shell grep "^ *[^%].*includegr" patchy-selection-paper.tex | sed -e 's/.*{\([^}]*\).*/\1/' )
PATCHY_EPS = $(patsubst %,%.eps,$(PATCHY_FIGS))
SFIG_EPS = $(patsubst %.pdf,%.eps,$(sfigs))

clean :
	-rm -f *.{aux,log,bbl,blg,fff,lof,lot,out,ttt}

figs : $(PATCHY_EPS) $(SFIG_EPS)

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

patchy-selection-paper-with-figs.pdf : patchy-selection-paper-submission-figs.pdf
	pdfjam --outfile $@ $< 1-42

patchy-selection-review-responses.pdf : patchy-selection-paper-submission-figs.pdf
	pdfjam --outfile $@ $< 43-

patchy-selection-paper-no-figs.pdf : patchy-selection-paper-submission-version.pdf
	pdfjam --outfile $@ $< 1-31

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


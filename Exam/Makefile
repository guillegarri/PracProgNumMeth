

.PHONY: all
all: regular irregular comparison project.pdf

project.pdf: project.tex regular.pdf irregular.pdf times.pdf
	pdflatex project
	pdflatex project

regular.pdf: regular.gpi regular
	gnuplot $<

irregular.pdf: irregular.gpi irregular
	gnuplot $<

times.pdf: time.gpi comparison
	gnuplot $<

.PHONY: regular
regular:
	$(MAKE) -C Regular_grid/

.PHONY: irregular
irregular:
	$(MAKE) -C Irregular_grid/

.PHONY: comparison
comparison:
	$(MAKE) -C Comparison/

.PHONY: cleanproject
cleanproject:
	rm -f *.pdf
	rm -f *.log
	rm -f *.aux


.PHONY: clean
clean: cleanproject
	cd Regular_grid/ && $(MAKE) clean
	cd Irregular_grid/ && $(MAKE) clean
	cd Comparison/ && $(MAKE) clean

CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONY: all

all: report.pdf plot.svg

report.pdf: report.tex plot.tex
	latex report
	latex report
	dvipdf report

plot.tex: texplot.gpi data.out
	gnuplot $<

plot.svg: svgplot.gpi data.out
	gnuplot $<

data.out:arctan
	./$< > $@

arctan:arctan.c

.PHONY: clean
clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.out
	rm -f *.svg
	rm -f plot.tex

CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
PARAMS = -3 3 0.2

.PHONY: all clean view

all: report.pdf

view: report.pdf
	xdg-open report.pdf

report.pdf: report.tex plot.tex
	latex report
	latex report
	dvipdf report

plot.tex:plot.gpi data.out
	gnuplot $<

data.out:errorfunction
	./$< $(PARAMS) > $@

errorfunction:errorfunction.c

clean:
	find . -type f -executable -delete
	rm *.out
	rm *.pdf
	rm plot.tex

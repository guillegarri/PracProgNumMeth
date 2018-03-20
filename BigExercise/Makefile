CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

plot.tex: texplot.gpi data.out
	gnuplot $<

plot.svg: svgplot.gpi data.out
	gnuplot $<

data.out:arctan
	./$< > $@

arctan:arctan.c

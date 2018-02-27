CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
MADEFILES = logistic logistic.svg logistic.dat GR uplot.svg orbit.dat circ.dat elip.dat GR.dat elip.svg circ.svg GR.svg

.PHONY: clean orbitplots all

all: logistic.svg uplot.svg elip.svg circ.svg GR.svg

logistic.svg:logistic.gpi logistic.dat
	gnuplot $<

logistic.dat:logistic
	./$< > $@

uplot.svg: uplot.gpi orbit.dat
	gnuplot $<

circ.svg: circ.gpi circ.dat
	gnuplot $<

elip.svg: elip.gpi elip.dat
	gnuplot $<

GR.svg: GR.gpi GR.dat
		gnuplot $<

circ.dat elip.dat GR.dat:orbit.dat
	sed -n 1,3145p orbit.dat > circ.dat
	sed -n 3148,6292p orbit.dat > elip.dat
	sed -n 6295,9439p orbit.dat > GR.dat

orbit.dat:GR
	./$< > $@

clean:
	rm -f $(MADEFILES)

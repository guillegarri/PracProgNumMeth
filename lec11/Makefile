CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: leastsquares.svg RosenPath.svg

leastsquares.svg: fit.out leastsquares.gpi
	gnuplot leastsquares.gpi

RosenPath.svg: RosenPath.out RosenPath.gpi
	gnuplot RosenPath.gpi

fit.out params.out: leastsquares
	./$< 1>fit.out 2>params.out

RosenMin.out RosenPath.out:rosenbrock
	./$< 1>RosenMin.out 2>RosenPath.out

rosenbrock: rosenbrock.c

leastsquares: leastsquares.c

clean:
	find . -type f -executable -delete
	rm *.out
	rm *.svg

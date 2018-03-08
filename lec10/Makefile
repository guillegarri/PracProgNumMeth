CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

rosenbrock.out: rosenbrock
	./$< > $@

rosenbrock: rosenbrock.c

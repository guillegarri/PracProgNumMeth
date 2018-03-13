CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all: wavefunc.svg rosenbrock.out

wavefunc.svg: epsilon.out plot.gpi
	gnuplot plot.gpi

epsilon.out log.out: hydrogen
	./$< 1> epsilon.out 2> log.out

hydrogen: hydrogen.c

rosenbrock.out: rosenbrock
	./$< > $@

rosenbrock: rosenbrock.c

clean:
		find . -type f -executable -delete
		rm *.out

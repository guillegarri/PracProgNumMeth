CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

svgs = rosenbrock_int.svg rosenbrock_real.svg himmelblau_int.svg himmelblau_real.svg 2dGauss_int.svg 2dGauss_real.svg linear_int.svg linear_real.svg gridpoints.svg

.PHONY: all
all: alldata.out $(svgs)

rosenbrock_int.svg: introsenbrock.gpi grid_rosenbrock_int.out
	gnuplot $<

rosenbrock_real.svg: realrosenbrock.gpi grid_rosenbrock.out
	gnuplot $<

himmelblau_int.svg: inthimmelblau.gpi grid_himmelblau_int.out
	gnuplot $<

himmelblau_real.svg: realhimmelblau.gpi grid_himmelblau.out
	gnuplot $<

2dGauss_int.svg: intdgauss.gpi grid_dgauss_int.out
	gnuplot $<

2dGauss_real.svg: realdgauss.gpi grid_dgauss.out
	gnuplot $<

linear_int.svg: intlinear.gpi grid_linear_int.out
	gnuplot $<

linear_real.svg: reallinear.gpi grid_linear.out
	gnuplot $<

gridpoints.svg: gridpoints.gpi XYvals.out
	gnuplot $<

alldata.out grid_rosenbrock.out grid_himmelblau.out grid_dgauss.out grid_rosenbrock_int.out grid_himmelblau_int.out grid_dgauss_int.out: main
	./main 1> alldata.out 2>log.out

main: main.o bilinear.o



clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.out
	rm -f *.svg

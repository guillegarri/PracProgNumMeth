CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

.PHONY: all
all: plot.svg plot2d.svg

plot2d.svg: plot2d.gpi results2d.out func2d.out
	gnuplot $<

plot.svg: plot.gpi results.out
	gnuplot $<

results2d.out func2d.out: main2d
	./$< > results2d.out 2> func2d.out

results.out: main
	./$< > $@ 2>log.out

main2d: main2d.o ann2d.o

main: main.o ann.o

clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.out
	rm -f *.svg

CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

plot.svg: plot.gpi results.out
	gnuplot $<

results.out: main
	./$< > $@ 2>log

main: main.o ann.o

clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.out

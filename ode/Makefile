CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

sine.svg: sine.gpi sine.out mathsine.out
	gnuplot $<

results.out sine.out: main
	./$< > results.out

mathsine.out: mathsine
	./$< > $@

mathsine: mathsine.c

main: main.o rungekutta.o

clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.out

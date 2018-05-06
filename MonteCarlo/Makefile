CFLAGS += -Wall -O -std=gnu1x
LDLIBS += -lm

error.svg: plot.gpi error.out
	gnuplot $<

results.out error.out: main
	./$< 1> results.out 2> error.out

main: main.o montecarlo.o

clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.out

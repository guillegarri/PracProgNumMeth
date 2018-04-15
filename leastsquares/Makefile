CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

plot1.svg: plot1.gpi data.out
	gnuplot $<

fit.out data.out: main
	./$< 1> fit.out 2> data.out

main: main.o gramschmidtfuncs.o

clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.out

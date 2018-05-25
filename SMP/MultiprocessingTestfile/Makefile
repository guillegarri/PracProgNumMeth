CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags` -fopenmp
LDLIBS += `gsl-config --libs` -lgomp

results.out: main
	./$< > $@

main: main.o integration.o

clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.out

CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

all:qdat.dat

intres.txt: myint
	./$< > $@

qdat.dat: quantum
	./$< 1> $@ 2>log.dat

clean:
	rm -f *.o *.dat
	find . -type f -executable -delete

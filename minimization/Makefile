CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

.PHONY: all
all: broyden.out newton.out

broyden.out: testbroyden
	./$< > $@

newton.out: main
	./$< > $@

testbroyden: testbroyden.o broyden.o

main: main.o gramschmidtfuncs.o newton.o

clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.out
